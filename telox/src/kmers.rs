/*
 * TeloX - Telomere Motif Extraction Tool
 *
 * Copyright (c) 2025 Yumi Sims Wellcome Sanger Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

use needletail::{parse_fastx_file, Sequence};
use std::collections::HashMap;
use std::path::Path;
use anyhow::{anyhow, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use csv::{ReaderBuilder, WriterBuilder};
use rayon::prelude::*;
use hashbrown::HashMap as Hm;
use ahash::RandomState;
use std::process::Command;

#[derive(Debug)]
pub struct KmerPair {
    pub forward: u32,
    pub rc: u32,
}

#[derive(Debug)]
pub struct ContinuityMetrics {
    pub window_coverage: f64,  // percentage of window covered by motifs
    pub motif_density: f64,    // motifs per base pair
    pub gap_count: usize,      // number of gaps between motifs
    pub max_gap_size: usize,   // largest gap between motifs
    pub avg_gap_size: f64,     // average gap size
}

#[derive(Debug)]
pub struct StrandBiasAnalysis {
    pub kmer: String,
    pub forward_count: u32,
    pub rc_count: u32,
    pub total_count: u32,
    pub bias_ratio: f64,
    pub bias_direction: String,  // "forward", "reverse", or "balanced"
    pub significance: String,    // "strong", "moderate", or "weak"
    pub longest_stretch: usize,  // longest continuous stretch in the sequence(s)
    pub longest_stretch_with_indels: usize,  // longest stretch allowing indels
    pub indel_tolerance_used: bool,  // whether indel tolerance was beneficial
}

fn is_homopolymer(seq: &str) -> bool {
    if seq.is_empty() {
        return false;
    }
    let first = seq.chars().next().unwrap();
    seq.chars().all(|c| c == first)
}

fn is_dinucleotide_repeat(seq: &str) -> bool {
    is_dinucleotide_repeat_bytes(seq.as_bytes())
}



fn is_simple_repeat(seq: &str) -> bool {
    is_simple_repeat_bytes(seq.as_bytes())
}

fn g_content(seq: &str) -> f64 {
    let g_count = seq.chars().filter(|&c| c == 'G' || c == 'g').count();
    g_count as f64 / seq.len() as f64
}

fn c_content(seq: &str) -> f64 {
    let c_count = seq.chars().filter(|&c| c == 'C' || c == 'c').count();
    c_count as f64 / seq.len() as f64
}

/// Generate reverse complement of a DNA sequence
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'N' | 'n' => 'N',
            _ => 'N',
        })
        .collect()
}

/// For each k-mer, find the longest continuous stretch (number of consecutive, non-overlapping occurrences)
/// in the sequence, checking both forward and reverse complement orientations, using a fast single-pass approach.
pub fn longest_continuous_stretch_for_kmers(
    sequence: &str,
    kmers: &[String],
) -> HashMap<String, usize> {
    let seq_bytes = sequence.as_bytes();
    if kmers.is_empty() || sequence.is_empty() {
        return HashMap::new();
    }
    let k = kmers[0].len();
    if k == 0 || k > seq_bytes.len() {
        return kmers.iter().map(|kmer| (kmer.clone(), 0)).collect();
    }

    // Map canonical k-mer to (set of all forms, current_run, max_run)
    let mut kmer_forms: HashMap<Vec<u8>, String> = HashMap::new();
    let mut run_map: HashMap<String, (usize, usize)> = HashMap::new();
    for kmer in kmers {
        let rc = reverse_complement(kmer);
        kmer_forms.insert(kmer.as_bytes().to_vec(), kmer.clone());
        kmer_forms.insert(rc.as_bytes().to_vec(), kmer.clone());
        run_map.insert(kmer.clone(), (0, 0));
    }

    let mut i = 0;
    let mut last_match: Option<String> = None;
    while i + k <= seq_bytes.len() {
        let window = &seq_bytes[i..i + k];
        if let Some(canonical) = kmer_forms.get(window) {
            // Continue run for this k-mer
            if last_match.as_ref() == Some(canonical) {
                let entry = run_map.get_mut(canonical).unwrap();
                entry.0 += 1;
                if entry.0 > entry.1 {
                    entry.1 = entry.0;
                }
            } else {
                // Reset all current runs except this one
                for (k, v) in run_map.iter_mut() {
                    if k == canonical {
                        v.0 = 1;
                        if v.0 > v.1 {
                            v.1 = v.0;
                        }
                    } else {
                        v.0 = 0;
                    }
                }
                last_match = Some(canonical.clone());
            }
            i += k;
        } else {
            // Reset all current runs
            for v in run_map.values_mut() {
                v.0 = 0;
            }
            last_match = None;
            i += 1;
        }
    }
    run_map.into_iter().map(|(k, (_cur, max))| (k, max)).collect()
}

/// Advanced longest stretch analysis allowing for small indels and gaps.
/// This function is more biologically relevant for telomeric sequences which may contain
/// sequencing errors, assembly gaps, or natural variations.
pub fn longest_continuous_stretch_with_indels(
    sequence: &str,
    kmers: &[String],
    max_indels: usize,      // Maximum number of indels allowed per stretch
    max_gap_size: usize,    // Maximum gap size before resetting (default: 5)
) -> HashMap<String, usize> {
    let seq_bytes = sequence.as_bytes();
    if kmers.is_empty() || sequence.is_empty() {
        return HashMap::new();
    }
    let k = kmers[0].len();
    if k == 0 || k > seq_bytes.len() {
        return kmers.iter().map(|kmer| (kmer.clone(), 0)).collect();
    }

    // Map canonical k-mer to all its forms (forward and reverse complement)
    let mut kmer_forms: HashMap<Vec<u8>, String> = HashMap::new();
    for kmer in kmers {
        let rc = reverse_complement(kmer);
        kmer_forms.insert(kmer.as_bytes().to_vec(), kmer.clone());
        kmer_forms.insert(rc.as_bytes().to_vec(), kmer.clone());
    }

    // Track state for each k-mer: (current_stretch, max_stretch, current_indels, last_match_pos)
    let mut kmer_states: HashMap<String, (usize, usize, usize, Option<usize>)> = HashMap::new();
    for kmer in kmers {
        kmer_states.insert(kmer.clone(), (0, 0, 0, None));
    }

    let mut i = 0;
    while i + k <= seq_bytes.len() {
        let window = &seq_bytes[i..i + k];
        
        if let Some(canonical) = kmer_forms.get(window) {
            // Found a k-mer match
            let state = kmer_states.get_mut(canonical).unwrap();
            
            match state.3 {
                None => {
                    // First match for this k-mer
                    state.0 = 1;
                    state.1 = state.1.max(1);
                    state.2 = 0;
                    state.3 = Some(i);
                }
                Some(last_pos) => {
                    // Calculate gap since last match
                    let expected_next = last_pos + k;
                    let gap_size = if i >= expected_next {
                        i - expected_next
                    } else {
                        0  // Overlapping matches, treat as continuous
                    };
                    
                    if gap_size == 0 {
                        // Perfect continuation
                        state.0 += 1;
                        state.1 = state.1.max(state.0);
                    } else if gap_size <= max_gap_size {
                        // Small gap - consume indel budget based on gap size
                        let indel_cost = match gap_size {
                            1..=2 => 1,
                            3..=5 => 2,
                            _ => gap_size / 2,  // Larger gaps cost more
                        };
                        
                        if state.2 + indel_cost <= max_indels {
                            // Can afford this gap
                            state.0 += 1;
                            state.1 = state.1.max(state.0);
                            state.2 += indel_cost;
                        } else {
                            // Too many indels, reset this k-mer's stretch
                            state.0 = 1;
                            state.2 = 0;
                        }
                    } else {
                        // Gap too large, reset
                        state.0 = 1;
                        state.2 = 0;
                    }
                    
                    state.3 = Some(i);
                }
            }
            
            // Reset other k-mers if they had gaps
            for (other_kmer, other_state) in kmer_states.iter_mut() {
                if other_kmer != canonical {
                    if let Some(other_last_pos) = other_state.3 {
                        let other_expected = other_last_pos + k;
                        let other_gap = if i >= other_expected {
                            i - other_expected
                        } else {
                            0
                        };
                        
                        // If gap is too large for the other k-mer, reset it
                        if other_gap > max_gap_size {
                            other_state.0 = 0;
                            other_state.2 = 0;
                            other_state.3 = None;
                        }
                    }
                }
            }
            
            i += k;  // Advance by k-mer length
        } else {
            // No match at this position, advance by 1
            // Update gap tracking for all k-mers
            for (_kmer, state) in kmer_states.iter_mut() {
                if let Some(last_pos) = state.3 {
                    let expected_next = last_pos + k;
                    let gap_size = if i >= expected_next {
                        i - expected_next + 1  // +1 for current mismatch
                    } else {
                        1
                    };
                    
                    // If gap becomes too large, reset
                    if gap_size > max_gap_size {
                        state.0 = 0;
                        state.2 = 0;
                        state.3 = None;
                    }
                }
            }
            
            i += 1;
        }
    }

    // Extract final max stretches
    kmer_states.into_iter()
        .map(|(kmer, (_current, max, _indels, _last_pos))| (kmer, max))
        .collect()
}

fn is_homopolymer_bytes(seq: &[u8]) -> bool {
    if seq.is_empty() {
        return false;
    }
    let first = seq[0];
    seq.iter().all(|&c| c == first)
}

fn g_content_bytes(seq: &[u8]) -> f64 {
    let g_count = seq.iter().filter(|&&c| c == b'G' || c == b'g').count();
    g_count as f64 / seq.len() as f64
}

fn c_content_bytes(seq: &[u8]) -> f64 {
    let c_count = seq.iter().filter(|&&c| c == b'C' || c == b'c').count();
    c_count as f64 / seq.len() as f64
}

fn is_dinucleotide_repeat_bytes(seq: &[u8]) -> bool {
    let len = seq.len();
    if len < 4 || len % 2 != 0 {
        return false;
    }
    let first = seq[0];
    let second = seq[1];
    if first == second {
        return false; // not a dinucleotide repeat if both are the same
    }
    for i in (0..len).step_by(2) {
        if seq[i] != first || seq[i + 1] != second {
            return false;
        }
    }
    true
}



fn is_simple_repeat_bytes(seq: &[u8]) -> bool {
    let len = seq.len();
    if len < 4 {
        return false;
    }
    
    // Check for repeats of length 1-4
    for repeat_len in 1..=4 {
        if len % repeat_len == 0 && len >= repeat_len * 2 {
            let pattern = &seq[0..repeat_len];
            let mut is_repeat = true;
            for i in (repeat_len..len).step_by(repeat_len) {
                if &seq[i..i + repeat_len] != pattern {
                    is_repeat = false;
                    break;
                }
            }
            if is_repeat {
                // Additional check: avoid flagging complex sequences as simple repeats
                let unique_bases: std::collections::HashSet<u8> = pattern.iter().cloned().collect();
                if unique_bases.len() <= 2 {
                    return true;
                }
            }
        }
    }
    false
}

/// Check if a k-mer is composed of tandem repeats of a shorter unit (up to k/2 length)
/// This helps identify cases like AACCGAACCG (which is 2x AACCG) that should be filtered out
fn is_tandem_repeat_kmer(seq: &str) -> bool {
    let len = seq.len();
    if len < 6 {  // Only check k-mers of 6+ bp
        return false;
    }
    
    // Check for tandem repeats of length 2 to k/2
    for repeat_len in 2..=(len/2) {
        if len % repeat_len == 0 && len >= repeat_len * 2 {
            let pattern = &seq[0..repeat_len];
            let mut is_repeat = true;
            for i in (repeat_len..len).step_by(repeat_len) {
                if &seq[i..i + repeat_len] != pattern {
                    is_repeat = false;
                    break;
                }
            }
            if is_repeat {
                return true;
            }
        }
    }
    false
}

pub fn count_kmers_in_fasta(
    fasta_path: impl AsRef<Path>,
    k: usize,
) -> Result<HashMap<String, KmerPair>> {
    if k == 0 {
        return Err(anyhow!("k must be greater than 0"));
    }
    if k > 32 {
        return Err(anyhow!("k must be <= 32 for fixed-size array optimization"));
    }

    let mut counts: Hm<[u8; 32], KmerPair, RandomState> = Hm::with_hasher(RandomState::new());
    let mut reader = parse_fastx_file(fasta_path)?;
    while let Some(record) = reader.next() {
        let seqrec = record?;
        seqrec.normalize(true);
        let seq_len = seqrec.seq().len();
        let min_pos = if seq_len > 5000 { seq_len - 5000 } else { 0 };
        for (pos, kmer) in seqrec.kmers(k as u8).enumerate() {
            if pos < min_pos {
                continue;
            }
            let mut forward = [0u8; 32];
            let mut rc = [0u8; 32];
            let kmer_vec = kmer.to_vec();
            let rc_vec = kmer.reverse_complement().to_vec();
            forward[..k].copy_from_slice(&kmer_vec);
            rc[..k].copy_from_slice(&rc_vec);
            if forward[..k].contains(&b'N') || rc[..k].contains(&b'N') {
                continue;
            }
            if is_homopolymer_bytes(&forward[..k]) {
                continue;
            }
            if is_dinucleotide_repeat_bytes(&forward[..k]) {
                continue;
            }
            if is_simple_repeat_bytes(&forward[..k]) {
                continue;
            }
            // Filter out tandem repeat k-mers (e.g., AACCGAACCG which is 2x AACCG)
            if let Ok(kmer_str) = std::str::from_utf8(&forward[..k]) {
                if is_tandem_repeat_kmer(kmer_str) {
                    continue;
                }
            }
            if forward[..k].iter().any(|&c| c.is_ascii_whitespace()) || rc[..k].iter().any(|&c| c.is_ascii_whitespace()) {
                continue;
            }
            if !(g_content_bytes(&forward[..k]) > 0.28 ||
                 c_content_bytes(&forward[..k]) > 0.28 ||
                 g_content_bytes(&rc[..k]) > 0.28 ||
                 c_content_bytes(&rc[..k]) > 0.28) {
                continue;
            }
            let (canonical, is_forward) = if forward[..k] < rc[..k] {
                (forward, true)
            } else if rc[..k] < forward[..k] {
                (rc, false)
            } else {
                (forward, true) // palindromic, treat as forward
            };
            let entry = counts.entry(canonical).or_insert(KmerPair { forward: 0, rc: 0 });
            if forward[..k] == rc[..k] {
                entry.forward += 1;
                entry.rc += 1;
            } else if is_forward {
                entry.forward += 1;
            } else {
                entry.rc += 1;
            }
        }
    }
    // Convert [u8; 32] keys to String for output
    let final_counts: HashMap<String, KmerPair> = counts.into_iter()
        .filter_map(|(kmer, pair)| {
            match std::str::from_utf8(&kmer[..k]) {
                Ok(s) => Some((s.to_string(), pair)),
                Err(_) => None,
            }
        })
        .collect();
    Ok(final_counts)
}

pub fn print_kmer_table(counts: &HashMap<String, KmerPair>, top_n: usize) {
    let mut sorted: Vec<_> = counts.iter().collect();
    sorted.sort_by(|a, b| {
        let total_a = a.1.forward + a.1.rc;
        let total_b = b.1.forward + b.1.rc;
        total_b.cmp(&total_a)  // Sort by total count descending
    });

    println!("\n{:<10} {:<15} {:<15} {:<15} {:<15}", 
             "Rank", "K-mer", "Total Count", "Forward", "RC");
    println!("{}", "-".repeat(70));

    for (i, (kmer, counts)) in sorted.iter().take(top_n).enumerate() {
        let total = counts.forward + counts.rc;
        println!("{:<10} {:<15} {:<15} {:<15} {:<15}",
                 i+1, kmer, total, counts.forward, counts.rc);
    }
}

pub fn analyze_strand_bias(
    counts: &HashMap<String, KmerPair>,
    longest_stretch_map: Option<&HashMap<String, usize>>,
) -> Vec<StrandBiasAnalysis> {
    analyze_strand_bias_with_indels(counts, longest_stretch_map, None)
}

pub fn analyze_strand_bias_with_indels(
    counts: &HashMap<String, KmerPair>,
    longest_stretch_map: Option<&HashMap<String, usize>>,
    longest_stretch_indels_map: Option<&HashMap<String, usize>>,
) -> Vec<StrandBiasAnalysis> {
    let mut bias_analyses = Vec::new();
    
    for (kmer, pair) in counts {
        let total = pair.forward + pair.rc;
        let bias_ratio = if pair.rc > 0 {
            pair.forward as f64 / pair.rc as f64
        } else {
            f64::INFINITY // All on forward strand
        };
        
        let bias_direction = if bias_ratio > 2.0 {
            "forward".to_string()
        } else if bias_ratio < 0.5 {
            "reverse".to_string()
        } else {
            "balanced".to_string()
        };
        
        let significance = if bias_ratio > 1.5 || bias_ratio < 0.66 {
            "strong".to_string()
        } else {
            "weak".to_string()
        };
        
        let longest_stretch = longest_stretch_map
            .and_then(|m| m.get(kmer))
            .copied()
            .unwrap_or(0);
            
        let longest_stretch_with_indels = longest_stretch_indels_map
            .and_then(|m| m.get(kmer))
            .copied()
            .unwrap_or(longest_stretch);  // Default to exact stretch if not provided
            
        let indel_tolerance_used = longest_stretch_with_indels > longest_stretch;
        
        bias_analyses.push(StrandBiasAnalysis {
            kmer: kmer.clone(),
            forward_count: pair.forward,
            rc_count: pair.rc,
            total_count: total,
            bias_ratio,
            bias_direction,
            significance,
            longest_stretch,
            longest_stretch_with_indels,
            indel_tolerance_used,
        });
    }
    
    // Sort by bias ratio (most biased first)
    bias_analyses.sort_by(|a, b| {
        let a_abs = if a.bias_ratio == f64::INFINITY { 1000.0 } else { a.bias_ratio };
        let b_abs = if b.bias_ratio == f64::INFINITY { 1000.0 } else { b.bias_ratio };
        b_abs.partial_cmp(&a_abs).unwrap_or(std::cmp::Ordering::Equal)
    });
    
    bias_analyses
}

pub fn print_strand_bias_table(analyses: &[StrandBiasAnalysis], top_n: usize) {
    println!("\n{:<15} {:<10} {:<10} {:<10} {:<12} {:<10} {:<10} {:<8} {:<8} {:<8}", 
             "K-mer", "Forward", "RC", "Total", "Bias Ratio", "Direction", "Significance", "Stretch", "StretchI", "Improved");
    println!("{}", "-".repeat(115));

    for analysis in analyses.iter().take(top_n) {
        let ratio_str = if analysis.bias_ratio == f64::INFINITY {
            "∞".to_string()
        } else {
            format!("{:.2}", analysis.bias_ratio)
        };
        
        let improved = if analysis.indel_tolerance_used { "Yes" } else { "No" };
        
        println!("{:<15} {:<10} {:<10} {:<10} {:<12} {:<10} {:<10} {:<8} {:<8} {:<8}",
                 analysis.kmer,
                 analysis.forward_count,
                 analysis.rc_count,
                 analysis.total_count,
                 ratio_str,
                 analysis.bias_direction,
                 analysis.significance,
                 analysis.longest_stretch,
                 analysis.longest_stretch_with_indels,
                 improved);
    }
}

pub fn save_strand_bias_table(analyses: &[StrandBiasAnalysis], output_path: &str) -> Result<()> {
    let mut content = String::new();
    content.push_str("Kmer\tForward\tRC\tTotal\tBiasRatio\tDirection\tSignificance\tLongestStretch\tLongestStretchIndels\tIndelImproved\n");
    
    for analysis in analyses {
        let ratio_str = if analysis.bias_ratio == f64::INFINITY {
            "inf".to_string()
        } else {
            format!("{:.3}", analysis.bias_ratio)
        };
        
        let improved = if analysis.indel_tolerance_used { "true" } else { "false" };
        
        content.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            analysis.kmer,
            analysis.forward_count,
            analysis.rc_count,
            analysis.total_count,
            ratio_str,
            analysis.bias_direction,
            analysis.significance,
            analysis.longest_stretch,
            analysis.longest_stretch_with_indels,
            improved
        ));
    }
    
    std::fs::write(output_path, content)?;
    Ok(())
}

pub fn filter_by_strand_bias(analyses: &[StrandBiasAnalysis], min_bias_ratio: f64) -> Vec<&StrandBiasAnalysis> {
    analyses.iter()
        .filter(|analysis| {
            if analysis.bias_ratio == f64::INFINITY {
                true // Include k-mers found only on forward strand
            } else {
                analysis.bias_ratio >= min_bias_ratio || analysis.bias_ratio <= (1.0 / min_bias_ratio)
            }
        })
        .collect()
}

pub fn get_strand_bias_summary(analyses: &[StrandBiasAnalysis]) -> (usize, usize, usize, usize) {
    let total = analyses.len();
    let forward_biased = analyses.iter().filter(|a| a.bias_direction == "forward").count();
    let reverse_biased = analyses.iter().filter(|a| a.bias_direction == "reverse").count();
    let balanced = analyses.iter().filter(|a| a.bias_direction == "balanced").count();
    
    (total, forward_biased, reverse_biased, balanced)
}

/// Filter multiple strand bias TSV files by LongestStretch >= 2 and Significance != "weak".
/// Writes all passing rows to a single output file, with header only once.

pub fn filter_strand_bias_tsvs(input_files: &[&str], output_file: &str) -> Result<()> {
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(BufWriter::new(File::create(output_file)?));
    let mut header_written = false;
    for &input_path in input_files {
        let file = File::open(input_path)?;
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(BufReader::new(file));
        let headers = reader.headers()?.clone();
        if !header_written {
            writer.write_record(&headers)?;
            header_written = true;
        }
        for result in reader.records() {
            let record = result?;
            let stretch_idx = headers.iter().position(|h| h == "LongestStretch").ok_or_else(|| anyhow!("No LongestStretch column in {}", input_path))?;
            let sig_idx = headers.iter().position(|h| h == "Significance").ok_or_else(|| anyhow!("No Significance column in {}", input_path))?;
            let stretch: usize = record.get(stretch_idx).unwrap_or("0").parse().unwrap_or(0);
            let significance = record.get(sig_idx).unwrap_or("").trim().to_lowercase();
            if stretch < 2 { continue; }
            if significance == "weak" { continue; }
            writer.write_record(&record)?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Given a list of k-mers, return a deduplicated list where each group of rotationally
/// equivalent k-mers is represented by the lexicographically smallest rotation.
pub fn consolidate_rotational_kmers(kmers: &[String]) -> Vec<String> {
    use std::collections::HashSet;

    fn min_rotation(s: &str) -> String {
        let k = s.len();
        let mut min = s.to_string();
        let doubled = s.repeat(2);
        for i in 1..k {
            let rot = &doubled[i..i + k];
            if rot < &min {
                min = rot.to_string();
            }
        }
        min
    }

    let mut seen = HashSet::new();
    let mut result = Vec::new();
    for kmer in kmers {
        let canonical = min_rotation(kmer);
        if seen.insert(canonical.clone()) {
            result.push(canonical);
        }
    }
    result
}

/// Filter bias analyses to keep only those with longest_stretch >= 2 and significance != "weak".
pub fn filter_bias_analyses<'a>(analyses: &'a [StrandBiasAnalysis]) -> Vec<&'a StrandBiasAnalysis> {
    analyses
        .iter()
        .filter(|a| a.longest_stretch >= 2 && a.significance != "weak")
        .collect()
}

/// Run KMC to count k-mers from a FASTA file
pub fn run_kmc(input_fasta: &str, k: usize, db_prefix: &str) -> std::io::Result<()> {
    let status = Command::new("kmc")
        .args(&[
            format!("-k{}", k),
            "-ci1".to_string(),
            "-cs1000000".to_string(),
            "-fm".to_string(),
            input_fasta.to_string(),
            db_prefix.to_string(),
            ".".to_string(),
        ])
        .status()?;
    if !status.success() {
        panic!("KMC failed with exit code: {:?}", status.code());
    }
    Ok(())
}

/// Run kmc_dump to export KMC results to a text file
pub fn run_kmc_dump(db_prefix: &str, output_txt: &str) -> std::io::Result<()> {
    let status = Command::new("kmc_dump")
        .args(&[db_prefix, output_txt])
        .status()?;
    if !status.success() {
        panic!("kmc_dump failed with exit code: {:?}", status.code());
    }
    Ok(())
}

/// Read dumped KMC k-mer counts from a text file
pub fn read_kmc_counts(path: &str) -> std::io::Result<Vec<(String, u64)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut kmers = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split_whitespace();
        if let (Some(kmer), Some(count)) = (parts.next(), parts.next()) {
            if let Ok(count) = count.parse() {
                kmers.push((kmer.to_string(), count));
            }
        }
    }
    Ok(kmers)
}

/// Full pipeline: run KMC, dump, and load k-mer counts
pub fn kmc_pipeline(input_fasta: &str, k: usize, db_prefix: &str, output_txt: &str) -> std::io::Result<Vec<(String, u64)>> {
    run_kmc(input_fasta, k, db_prefix)?;
    run_kmc_dump(db_prefix, output_txt)?;
    read_kmc_counts(output_txt)
}

/// Filter and analyze k-mers loaded from KMC output, applying biological filters and bias analysis.
pub fn analyze_kmc_kmers(
    kmers: Vec<(String, u64)>,
    k: usize,
) -> Vec<StrandBiasAnalysis> {
    use std::collections::HashMap;
    let mut canonical_counts: HashMap<String, KmerPair> = HashMap::new();

    for (kmer, count) in kmers {
        if kmer.len() != k {
            continue;
        }
        let rc = reverse_complement(&kmer);
        let (canonical, is_forward) = if kmer < rc {
            (kmer.clone(), true)
        } else {
            (rc.clone(), false)
        };

        // Apply filters to canonical k-mer
        if canonical.contains('N') || reverse_complement(&canonical).contains('N') {
            continue;
        }
        if is_homopolymer(&canonical) {
            continue;
        }
        if is_dinucleotide_repeat_bytes(canonical.as_bytes()) {
            continue;
        }
        if is_simple_repeat_bytes(canonical.as_bytes()) {
            continue;
        }
        // Filter out tandem repeat k-mers
        if is_tandem_repeat_kmer(&canonical) {
            continue;
        }
        if canonical.chars().any(|c| c.is_whitespace()) || reverse_complement(&canonical).chars().any(|c| c.is_whitespace()) {
            continue;
        }
        if !(g_content(&canonical) > 0.28 || c_content(&canonical) > 0.28 ||
             g_content(&reverse_complement(&canonical)) > 0.28 ||
             c_content(&reverse_complement(&canonical)) > 0.28) {
            continue;
        }

        let entry = canonical_counts.entry(canonical).or_insert(KmerPair { forward: 0, rc: 0 });
        if is_forward {
            entry.forward += count as u32;
        } else {
            entry.rc += count as u32;
        }
    }

    analyze_strand_bias(&canonical_counts, None)
}

/// Optimized k-mer counting for the last N bp of each scaffold using parallel processing
/// Only processes scaffolds larger than 1MB (1,000,000 bp)
pub fn count_kmers_last_n_bp_parallel(
    fasta_path: impl AsRef<Path>,
    k: usize,
    n: usize,
) -> Result<HashMap<String, KmerPair>> {
    if k == 0 {
        return Err(anyhow!("k must be greater than 0"));
    }
    if k > 32 {
        return Err(anyhow!("k must be <= 32 for fixed-size array optimization"));
    }

    // Read all sequences first, filtering by size
    let mut reader = parse_fastx_file(fasta_path)?;
    let mut sequences = Vec::new();
    let mut filtered_count = 0;
    let mut processed_count = 0;
    
    while let Some(record) = reader.next() {
        let seqrec = record?;
        seqrec.normalize(true);
        let seq = seqrec.seq().to_vec();
        
        // Only process scaffolds larger than 1MB (1,000,000 bp)
        if seq.len() < 1_000_000 {
            filtered_count += 1;
            continue;
        }
        
        processed_count += 1;
        sequences.push(seq);
    }
    
    eprintln!("K-mer counting: processing {} scaffolds >= 1MB, filtered out {} smaller scaffolds", processed_count, filtered_count);

    // Process sequences in parallel
    let results: Vec<HashMap<[u8; 32], KmerPair, RandomState>> = sequences.par_iter().map(|seq| {
        let mut local_counts: HashMap<[u8; 32], KmerPair, RandomState> = HashMap::with_hasher(RandomState::new());
        
        // Get the last N bp slice
        let seq_len = seq.len();
        let start_pos = if seq_len > n { seq_len - n } else { 0 };
        let seq_slice = &seq[start_pos..];
        
        // Count k-mers in the slice
        for window in seq_slice.windows(k) {
            let mut forward = [0u8; 32];
            let mut rc = [0u8; 32];
            
            // Copy k-mer to forward array
            forward[..k].copy_from_slice(window);
            
            // Generate reverse complement
            for (i, &base) in window.iter().enumerate() {
                rc[k - 1 - i] = match base {
                    b'A' | b'a' => b'T',
                    b'T' | b't' => b'A',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    b'N' | b'n' => b'N',
                    _ => b'N',
                };
            }
            
            // Apply filters
            if forward[..k].contains(&b'N') || rc[..k].contains(&b'N') {
                continue;
            }
            if is_homopolymer_bytes(&forward[..k]) {
                continue;
            }
            if is_dinucleotide_repeat_bytes(&forward[..k]) {
                continue;
            }
            if is_simple_repeat_bytes(&forward[..k]) {
                continue;
            }
            // Filter out tandem repeat k-mers (e.g., AACCGAACCG which is 2x AACCG)
            if let Ok(kmer_str) = std::str::from_utf8(&forward[..k]) {
                if is_tandem_repeat_kmer(kmer_str) {
                    continue;
                }
            }
            if forward[..k].iter().any(|&c| c.is_ascii_whitespace()) || rc[..k].iter().any(|&c| c.is_ascii_whitespace()) {
                continue;
            }
            if !(g_content_bytes(&forward[..k]) > 0.28 ||
                 c_content_bytes(&forward[..k]) > 0.28 ||
                 g_content_bytes(&rc[..k]) > 0.28 ||
                 c_content_bytes(&rc[..k]) > 0.28) {
                continue;
            }
            
            // Determine canonical form
            let (canonical, is_forward) = if forward[..k] < rc[..k] {
                (forward, true)
            } else if rc[..k] < forward[..k] {
                (rc, false)
            } else {
                (forward, true) // palindromic
            };
            
            let entry = local_counts.entry(canonical).or_insert(KmerPair { forward: 0, rc: 0 });
            if forward[..k] == rc[..k] {
                entry.forward += 1;
                entry.rc += 1;
            } else if is_forward {
                entry.forward += 1;
            } else {
                entry.rc += 1;
            }
        }
        local_counts
    }).collect();

    // Merge all local results
    let mut final_counts: HashMap<[u8; 32], KmerPair, RandomState> = HashMap::with_hasher(RandomState::new());
    for local in results {
        for (kmer, pair) in local {
            let entry = final_counts.entry(kmer).or_insert(KmerPair { forward: 0, rc: 0 });
            entry.forward += pair.forward;
            entry.rc += pair.rc;
        }
    }

    // Convert to String keys for output
    let string_counts: HashMap<String, KmerPair> = final_counts.into_iter()
        .filter_map(|(kmer, pair)| {
            match std::str::from_utf8(&kmer[..k]) {
                Ok(s) => Some((s.to_string(), pair)),
                Err(_) => None,
            }
        })
        .collect();
    
    Ok(string_counts)
}

/// Extract the last N bp of each scaffold/sequence from a FASTA file and write to a new FASTA file.
/// Only processes scaffolds larger than 1MB and merges them into a single sequence separated by spaces.
pub fn extract_last_n_bp_to_fasta(input_fasta: &str, output_fasta: &str, n: usize) -> std::io::Result<()> {
    let mut reader = parse_fastx_file(input_fasta)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, format!("{}", e)))?;
    let mut writer = std::fs::File::create(output_fasta)?;
    
    let mut merged_sequence = Vec::new();
    let mut scaffold_count = 0;
    let mut filtered_count = 0;
    
    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, format!("{}", e)))?;
        let seq = seqrec.seq();
        let len = seq.len();
        
        // Only process scaffolds larger than 1MB (1,000,000 bp)
        if len < 1_000_000 {
            filtered_count += 1;
            continue;
        }
        
        scaffold_count += 1;
        let start = if len > n { len - n } else { 0 };
        let last_n = &seq[start..];
        
        // Add space separator if this is not the first sequence
        if !merged_sequence.is_empty() {
            merged_sequence.push(b' ');
        }
        
        // Add the sequence
        merged_sequence.extend_from_slice(last_n);
    }
    
    // Write the merged sequence as a single FASTA entry
    writeln!(writer, ">merged_large_scaffolds_last{}bp_count{}", n, scaffold_count)?;
    
    // Write sequence in 60bp lines
    for chunk in merged_sequence.chunks(60) {
        writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap_or(""))?;
    }
    
    eprintln!("Processed {} scaffolds >= 1MB, filtered out {} smaller scaffolds", scaffold_count, filtered_count);
    eprintln!("Total merged sequence length: {} bp", merged_sequence.len());
    
    Ok(())
}

/// Improved: Parse KMC -b output and match Rust-native k-mer/RC counting logic
pub fn parse_kmc_b_output_and_analyze(
    path: &str,
    k: usize,
    longest_stretch_map: Option<&std::collections::HashMap<String, usize>>,
) -> std::io::Result<Vec<StrandBiasAnalysis>> {
    use std::collections::HashMap;
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    // First, build a map of all k-mers and their (count, rc_count)
    let mut raw_counts: HashMap<String, (u32, u32)> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split_whitespace();
        let kmer = match parts.next() {
            Some(s) => s.to_string(),
            None => continue,
        };
        let count = match parts.next() {
            Some(s) => s.parse::<u32>().unwrap_or(0),
            None => continue,
        };
        let rc_count = match parts.next() {
            Some(s) => s.parse::<u32>().unwrap_or(0),
            None => 0,
        };
        if kmer.len() != k {
            continue;
        }
        raw_counts.insert(kmer, (count, rc_count));
    }
    // Now, for each canonical k-mer, sum counts from both k-mer and its RC
    let mut canonical_counts: HashMap<String, KmerPair> = HashMap::new();
    let mut seen = std::collections::HashSet::new();
    for (kmer, &(count, rc_count)) in &raw_counts {
        let rc = reverse_complement(kmer);
        let canonical = if kmer < &rc { kmer } else { &rc };
        if seen.contains(canonical) {
            continue;
        }
        seen.insert(canonical.clone());
        // Get counts for both kmer and rc
        let (fwd1, rc1) = raw_counts.get(kmer).copied().unwrap_or((0, 0));
        let (fwd2, rc2) = raw_counts.get(&rc).copied().unwrap_or((0, 0));
        // Assign forward/rc as in Rust-native: forward = fwd1 + rc2, rc = rc1 + fwd2
        let forward = fwd1 + rc2;
        let rc_val = rc1 + fwd2;
        // Apply same biological filters as in Rust-native pipeline
        if canonical.contains('N') || reverse_complement(canonical).contains('N') {
            continue;
        }
        if is_homopolymer(canonical) {
            continue;
        }
        if is_dinucleotide_repeat_bytes(canonical.as_bytes()) {
            continue;
        }
        if is_simple_repeat_bytes(canonical.as_bytes()) {
            continue;
        }
        // Filter out tandem repeat k-mers
        if is_tandem_repeat_kmer(&canonical) {
            continue;
        }
        if canonical.chars().any(|c| c.is_whitespace()) || reverse_complement(canonical).chars().any(|c| c.is_whitespace()) {
            continue;
        }
        if !(g_content(canonical) > 0.28 || c_content(canonical) > 0.28 ||
             g_content(&reverse_complement(canonical)) > 0.28 ||
             c_content(&reverse_complement(canonical)) > 0.28) {
            continue;
        }
        canonical_counts.insert(canonical.to_string(), KmerPair { forward, rc: rc_val });
    }
    Ok(analyze_strand_bias(&canonical_counts, longest_stretch_map))
}

pub fn gather_and_rank_filtered_results(filtered_files: &[&str], output_file: &str) -> Result<()> {
    #[derive(Debug)]
    struct RankedKmer {
        kmer: String,
        forward: u32,
        rc: u32,
        total: u32,
        bias_ratio: String,
        direction: String,
        significance: String,
        longest_stretch: usize,
    }

    let mut all_kmers = Vec::new();

    for &file in filtered_files {
        let f = File::open(file)?;
        let reader = BufReader::new(f);
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if i == 0 || line.trim().is_empty() { continue; } // skip header/empty
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 8 { continue; }
            all_kmers.push(RankedKmer {
                kmer: cols[0].to_string(),
                forward: cols[1].parse().unwrap_or(0),
                rc: cols[2].parse().unwrap_or(0),
                total: cols[3].parse().unwrap_or(0),
                bias_ratio: cols[4].to_string(),
                direction: cols[5].to_string(),
                significance: cols[6].to_string(),
                longest_stretch: cols[7].parse().unwrap_or(0),
            });
        }
    }

    // Sort by longest_stretch DESC, then total DESC
    all_kmers.sort_by(|a, b| {
        b.longest_stretch.cmp(&a.longest_stretch)
            .then(b.total.cmp(&a.total))
    });

    // Write to output
    let mut out = File::create(output_file)?;
    writeln!(out, "Kmer\tForward\tRC\tTotal\tBiasRatio\tDirection\tSignificance\tLongestStretch")?;
    for k in all_kmers {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            k.kmer, k.forward, k.rc, k.total, k.bias_ratio, k.direction, k.significance, k.longest_stretch
        )?;
    }

    Ok(())
}

/// Given a rank.tsv file, return a Vec<String> of unique motifs consolidated by their minimal rotation (lex smallest rotation), like TELO_MOTIF_DB.
pub fn consolidate_ranked_motifs_by_rotation(rank_tsv: &str) -> anyhow::Result<Vec<String>> {
    use std::collections::HashMap;
    
    #[derive(Debug)]
    struct MotifData {
        kmer: String,
        forward: u32,
        reverse: u32,
        longest_stretch: usize,
    }
    
    let file = std::fs::File::open(rank_tsv)?;
    let reader = std::io::BufReader::new(file);
    let mut motif_data = Vec::new();
    
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 || line.trim().is_empty() { continue; } // skip header/empty
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 { continue; }
        
        let kmer = cols[0].to_string();
        let forward: u32 = cols[1].parse().unwrap_or(0);
        let reverse: u32 = cols[2].parse().unwrap_or(0);
        let longest_stretch: usize = cols[7].parse().unwrap_or(0);
        
        motif_data.push(MotifData {
            kmer,
            forward,
            reverse,
            longest_stretch,
        });
    }
    
    // Helper function to find minimal rotation
    fn min_rotation(s: &str) -> String {
        let k = s.len();
        let mut min = s.to_string();
        let doubled = s.repeat(2);
        for i in 1..k {
            let rot = &doubled[i..i + k];
            if rot < &min {
                min = rot.to_string();
            }
        }
        min
    }
    
    // Group motifs by their canonical (minimal) rotation
    let mut rotation_groups: HashMap<String, Vec<MotifData>> = HashMap::new();
    for data in motif_data {
        let canonical = min_rotation(&data.kmer);
        rotation_groups.entry(canonical).or_insert_with(Vec::new).push(data);
    }
    
    // Create consolidated groups with combined statistics
    let mut consolidated_groups = Vec::new();
    for (canonical, variants) in rotation_groups {
        let all_variants: Vec<String> = variants.iter().map(|v| v.kmer.clone()).collect();
        let forward_total: u32 = variants.iter().map(|v| v.forward).sum();
        let reverse_total: u32 = variants.iter().map(|v| v.reverse).sum();
        let max_longest_stretch: usize = variants.iter().map(|v| v.longest_stretch).max().unwrap_or(0);
        
        consolidated_groups.push((canonical, all_variants, forward_total, reverse_total, max_longest_stretch));
    }
    
    // Sort by max longest stretch (descending), then by total frequency (descending)
    consolidated_groups.sort_by(|a, b| b.4.cmp(&a.4).then((b.2 + b.3).cmp(&(a.2 + a.3))));
    
    // Print consolidated rotational groups table
    println!("\nConsolidated Rotational Groups:");
    println!("{}", "=".repeat(90));
    println!("{:<15} {:<30} {:<15} {:<15} {:<15}", 
             "Canonical_kmer", "kmer_group", "forward_total", "reverse_total", "longeststretch");
    println!("{}", "-".repeat(90));
    
    for (canonical, variants, forward_total, reverse_total, max_stretch) in &consolidated_groups {
        let variants_str = variants.join(",");
        println!("{:<15} {:<30} {:<15} {:<15} {:<15}", 
                 canonical, variants_str, forward_total, reverse_total, max_stretch);
    }
    
    // Identify the most likely telomere motif
    let most_likely_motif = if let Some((canonical, variants, forward_total, reverse_total, max_stretch)) = consolidated_groups.first() {
        let total_frequency = forward_total + reverse_total;
        
        // Create telomere motif result
        let telomere_result = serde_json::json!({
            "most_likely_telomere_motif": {
                "canonical_sequence": canonical,
                "rotational_variants": variants,
                "statistics": {
                    "forward_count": forward_total,
                    "reverse_count": reverse_total,
                    "total_frequency": total_frequency,
                    "longest_continuous_stretch": max_stretch
                },
                "analysis_rationale": format!(
                    "Selected as most likely telomere motif based on highest longest stretch ({}) and total frequency ({})",
                    max_stretch, total_frequency
                )
            },
            "analysis_metadata": {
                "total_motif_groups_analyzed": consolidated_groups.len(),
                "selection_criteria": ["longest_continuous_stretch", "total_frequency"],
                "data_source": "k-mer_discovery_pipeline"
            }
        });
        
        // Write to JSON file
        match std::fs::write("telomere_motif_final.json", serde_json::to_string_pretty(&telomere_result).unwrap()) {
            Ok(_) => {
                println!("\nMost likely telomere motif identified:");
                println!("Canonical sequence: {}", canonical);
                println!("Rotational variants: {}", variants.join(", "));
                println!("Total frequency: {} (Forward: {}, Reverse: {})", total_frequency, forward_total, reverse_total);
                println!("Longest stretch: {}", max_stretch);
                println!("Result written to: telomere_motif_final.json");
            },
            Err(e) => {
                eprintln!("Error writing telomere_motif_final.json: {}", e);
            }
        }
        
        Some(canonical.clone())
    } else {
        println!("No telomere motifs found in consolidated results");
        None
    };

    // Return unique canonical representatives for the motif database
    let result: Vec<String> = consolidated_groups.into_iter()
        .map(|(canonical, _, _, _, _)| canonical)
        .collect();
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_indel_tolerant_stretch() {
        // Test sequence with a telomere-like pattern with gaps
        let sequence = "TTAGGGTTAGGGTTACGGTTAGGGTTAGGG"; // TTAGGG with 2bp gap
        let kmers = vec!["TTAGGG".to_string()];
        
        // Exact matching should give shorter stretch
        let exact_result = longest_continuous_stretch_for_kmers(sequence, &kmers);
        
        // Indel-tolerant should give longer stretch
        let indel_result = longest_continuous_stretch_with_indels(sequence, &kmers, 5, 5);
        
        println!("Test sequence: {}", sequence);
        println!("Exact stretch: {:?}", exact_result);
        println!("Indel-tolerant stretch: {:?}", indel_result);
        
        // The indel-tolerant version should handle the 2bp gap better
        assert!(indel_result.get("TTAGGG").unwrap_or(&0) >= exact_result.get("TTAGGG").unwrap_or(&0));
    }
    
    #[test]
    fn test_strict_vs_indel_tolerant() {
        // Sequence with perfect telomere repeats plus some with gaps
        let sequence = "TTAGGGTTAGGGTTAGGGNNNAGGGTTAGGGTTAGGG";
        let kmers = vec!["TTAGGG".to_string()];
        
        let exact = longest_continuous_stretch_for_kmers(sequence, &kmers);
        let indel = longest_continuous_stretch_with_indels(sequence, &kmers, 3, 3);
        
        println!("Sequence with gaps: {}", sequence);
        println!("Exact: {:?}", exact);
        println!("Indel-tolerant: {:?}", indel);
        
        // Indel-tolerant should be >= exact
        assert!(indel.get("TTAGGG").unwrap_or(&0) >= exact.get("TTAGGG").unwrap_or(&0));
    }
}
