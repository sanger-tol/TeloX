/*
 * TeloX - Telomere Motif Extraction Tool
 *
 * Copyright (c) 2024 [Your Name]
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
}

fn is_homopolymer(seq: &str) -> bool {
    if seq.is_empty() {
        return false;
    }
    let first = seq.chars().next().unwrap();
    seq.chars().all(|c| c == first)
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
    use std::collections::HashSet;
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
        
        bias_analyses.push(StrandBiasAnalysis {
            kmer: kmer.clone(),
            forward_count: pair.forward,
            rc_count: pair.rc,
            total_count: total,
            bias_ratio,
            bias_direction,
            significance,
            longest_stretch,
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
    println!("\n{:<15} {:<10} {:<10} {:<10} {:<12} {:<10} {:<10}", 
             "K-mer", "Forward", "RC", "Total", "Bias Ratio", "Direction", "Significance");
    println!("{}", "-".repeat(85));

    for analysis in analyses.iter().take(top_n) {
        let ratio_str = if analysis.bias_ratio == f64::INFINITY {
            "∞".to_string()
        } else {
            format!("{:.2}", analysis.bias_ratio)
        };
        
        println!("{:<15} {:<10} {:<10} {:<10} {:<12} {:<10} {:<10}",
                 analysis.kmer,
                 analysis.forward_count,
                 analysis.rc_count,
                 analysis.total_count,
                 ratio_str,
                 analysis.bias_direction,
                 analysis.significance);
    }
}

pub fn save_strand_bias_table(analyses: &[StrandBiasAnalysis], output_path: &str) -> Result<()> {
    let mut content = String::new();
    content.push_str("Kmer\tForward\tRC\tTotal\tBiasRatio\tDirection\tSignificance\n");
    
    for analysis in analyses {
        let ratio_str = if analysis.bias_ratio == f64::INFINITY {
            "inf".to_string()
        } else {
            format!("{:.3}", analysis.bias_ratio)
        };
        
        content.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            analysis.kmer,
            analysis.forward_count,
            analysis.rc_count,
            analysis.total_count,
            ratio_str,
            analysis.bias_direction,
            analysis.significance
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

/// Optimized k-mer counting for the last 5000bp of each scaffold using parallel processing
pub fn count_kmers_last_5000bp_parallel(
    fasta_path: impl AsRef<Path>,
    k: usize,
) -> Result<HashMap<String, KmerPair>> {
    if k == 0 {
        return Err(anyhow!("k must be greater than 0"));
    }
    if k > 32 {
        return Err(anyhow!("k must be <= 32 for fixed-size array optimization"));
    }

    // Read all sequences first
    let mut reader = parse_fastx_file(fasta_path)?;
    let mut sequences = Vec::new();
    while let Some(record) = reader.next() {
        let seqrec = record?;
        seqrec.normalize(true);
        sequences.push(seqrec.seq().to_vec());
    }

    // Process sequences in parallel
    let results: Vec<HashMap<[u8; 32], KmerPair, RandomState>> = sequences.par_iter().map(|seq| {
        let mut local_counts: HashMap<[u8; 32], KmerPair, RandomState> = HashMap::with_hasher(RandomState::new());
        
        // Get the last 5000bp slice
        let seq_len = seq.len();
        let start_pos = if seq_len > 5000 { seq_len - 5000 } else { 0 };
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
pub fn extract_last_n_bp_to_fasta(input_fasta: &str, output_fasta: &str, n: usize) -> std::io::Result<()> {
    let mut reader = parse_fastx_file(input_fasta)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, format!("{}", e)))?;
    let mut writer = std::fs::File::create(output_fasta)?;
    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, format!("{}", e)))?;
        let name = std::str::from_utf8(seqrec.id()).unwrap_or("scaffold");
        let seq = seqrec.seq();
        let len = seq.len();
        let start = if len > n { len - n } else { 0 };
        let last_n = &seq[start..];
        // Write FASTA header
        writeln!(writer, ">{}", name)?;
        // Write sequence in 60bp lines
        for chunk in last_n.chunks(60) {
            writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap_or(""))?;
        }
    }
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
    use std::collections::HashSet;
    let file = std::fs::File::open(rank_tsv)?;
    let reader = std::io::BufReader::new(file);
    let mut motifs = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 || line.trim().is_empty() { continue; } // skip header/empty
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.is_empty() { continue; }
        motifs.push(cols[0].to_string());
    }
    // Use min_rotation logic from consolidate_rotational_kmers
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
    for motif in motifs {
        let canonical = min_rotation(&motif);
        if seen.insert(canonical.clone()) {
            result.push(canonical);
        }
    }
    Ok(result)
}
