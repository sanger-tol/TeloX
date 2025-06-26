use needletail::{parse_fastx_file, Sequence};
use std::collections::HashMap;
use std::path::Path;
use anyhow::{anyhow, Result};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use csv::{ReaderBuilder, WriterBuilder};
use rayon::prelude::*;

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

pub fn count_kmers_in_fasta(
    fasta_path: impl AsRef<Path>,
    k: usize,
) -> Result<HashMap<String, KmerPair>> {
    if k == 0 {
        return Err(anyhow!("k must be greater than 0"));
    }

    use std::collections::HashMap as StdHashMap;
    let mut counts: StdHashMap<Vec<u8>, KmerPair> = StdHashMap::new();
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
            let forward = kmer.to_vec();
            let rc = kmer.reverse_complement().to_vec();
            if forward.contains(&b'N') || rc.contains(&b'N') {
                continue;
            }
            if is_homopolymer_bytes(&forward) {
                continue;
            }
            if forward.iter().any(|&c| c.is_ascii_whitespace()) || rc.iter().any(|&c| c.is_ascii_whitespace()) {
                continue;
            }
            if !(g_content_bytes(&forward) > 0.28 ||
                 c_content_bytes(&forward) > 0.28 ||
                 g_content_bytes(&rc) > 0.28 ||
                 c_content_bytes(&rc) > 0.28) {
                continue;
            }
            let (canonical, is_forward) = if forward < rc {
                (forward.clone(), true)
            } else if rc < forward {
                (rc.clone(), false)
            } else {
                (forward.clone(), true) // palindromic, treat as forward
            };
            let entry = counts.entry(canonical).or_insert(KmerPair { forward: 0, rc: 0 });
            if forward == rc {
                entry.forward += 1;
                entry.rc += 1;
            } else if is_forward {
                entry.forward += 1;
            } else {
                entry.rc += 1;
            }
        }
    }
    // Convert Vec<u8> keys to String for output
    let final_counts: HashMap<String, KmerPair> = counts.into_iter()
        .filter_map(|(kmer, pair)| {
            match String::from_utf8(kmer) {
                Ok(s) => Some((s, pair)),
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

pub fn save_kmer_table(counts: &HashMap<String, KmerPair>, output_path: &str) -> Result<()> {
    let mut entries: Vec<_> = counts.iter().collect();

    // Sort by total count (forward + rc), descending
    entries.sort_by(|a, b| {
        let total_a = a.1.forward + a.1.rc;
        let total_b = b.1.forward + b.1.rc;
        total_b.cmp(&total_a)
    });

    let content = entries.iter()
        .map(|(kmer, counts)| {
            format!(
                "{}\t{}\t{}\t{}",
                kmer,
                counts.forward + counts.rc,
                counts.forward,
                counts.rc
            )
        })
        .collect::<Vec<_>>()
        .join("\n");

    std::fs::write(output_path, content)?;
    Ok(())
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
        
        let significance = if bias_ratio > 5.0 || bias_ratio < 0.2 {
            "strong".to_string()
        } else if bias_ratio > 2.0 || bias_ratio < 0.5 {
            "moderate".to_string()
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
