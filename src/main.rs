/*
 * TeloX - Telomere Motif Extraction Tool
 *
 * Copyright (c) 2025 Yumi Sims, Wellcome Sanger Institute
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

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, Write, BufReader};
use std::path::Path;
use std::env;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use flate2::read::GzDecoder;
mod kmers;
mod annotation;
use kmers::{
    count_kmers_last_n_bp_parallel, analyze_strand_bias, analyze_strand_bias_with_indels, 
    consolidate_ranked_motifs_by_rotation, extract_last_n_bp_to_fasta, run_kmc, run_kmc_dump, 
    read_kmc_counts, analyze_kmc_kmers, longest_continuous_stretch_for_kmers, 
    longest_continuous_stretch_with_indels, calculate_kmer_window_density
};
use annotation::{annotate_genome_with_motifs, create_annotation_config};

// Constants
const TELO_PENALTY: i64 = 1;
const TELO_MAX_DROP: i64 = 2000;
const TELO_MIN_SCORE: i64 = 300;
static TELO_MOTIF_DB: &[&str] = &[
    "AAAATTGTCCGTCC",
    "AAACCACCCT",
    "AAACCC",
    "AAACCCC",
    "AAACCCT",
    "AAAGAACCT",
    "AAATGTGGAGG",
    "AACAGACCCG",
    "AACCATCCCT",
    "AACCC",
    "AACCCAGACCC",
    "AACCCAGACCT",
    "AACCCAGACGC",
    "AACCCCAACCT",
    "AACCCGAACCT",
    "AACCCT",
    "AACCCTG",
    "AACCCTGACGC",
    "AACCT",
    "AAGGAC",
    "ACCCAG",
    "ACCTG",
    "ACGGCAGCG",
];


// Nucleotide table (similar to seq_nt6_table)
const SEQ_NT6_TABLE: [u8; 256] = [
    0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
];

struct Sequence {
    name: String,
    seq: String,
}

struct SequenceDict {
    sequences: Vec<Sequence>,
}

#[derive(Serialize, Deserialize)]
struct MotifFrequency {
    motif: String,
    frequency: usize,
    percentage: f64,
}

#[derive(Serialize, Deserialize)]
struct TelomereCandidate {
    status: String,
    primary_motif: Option<String>,
    frequency: Option<usize>,
    percentage: Option<f64>,
    source: Option<String>,
    total_occurrences: Option<usize>,
    unique_motifs: Option<usize>,
    top_motifs: Vec<MotifFrequency>,
    message: Option<String>,
}

impl SequenceDict {
    fn from_fasta<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(&path)?;
        let path_str = path.as_ref().to_string_lossy();
        
        // Create appropriate reader based on file extension
        let reader: Box<dyn BufRead> = if path_str.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        
        let mut sequences = Vec::new();
        let mut current_name = String::new();
        let mut current_seq = String::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                if !current_name.is_empty() {
                    sequences.push(Sequence {
                        name: current_name.clone(),
                        seq: current_seq.clone(),
                    });
                    current_seq.clear();
                }
                current_name = line[1..].trim().to_string();
            } else {
                current_seq.push_str(&line.trim().to_uppercase());
            }
        }

        if !current_name.is_empty() {
            sequences.push(Sequence {
                name: current_name,
                seq: current_seq,
            });
        }

        Ok(Self { sequences })
    }
}

fn check_motif(motif: &str) -> bool {
    for c in motif.chars() {
        let val = SEQ_NT6_TABLE[c as usize];
        if val < 1 || val > 4 {
            return false;
        }
    }
    true
}

fn list_telo_motifs(out: &mut impl Write) -> io::Result<()> {
    for (i, motif) in TELO_MOTIF_DB.iter().enumerate() {
        writeln!(out, "[{:2}] {}", i + 1, motif)?;
    }
    Ok(())
}

fn telo_finder_core(
    sequence: &str,
    mtab: &HashSet<u64>,
    mlen: usize,
    penalty: i64,
    max_drop: i64,
    min_score: i64,
) -> u64 {
    let slen = sequence.len();
    let mut sum_telo = 0u64;
    let mask = (1u64 << (2 * mlen)) - 1;

    // Check 5' end
    let mut score = 0i64;
    let mut max_score = 0i64;
    let mut max_i = -1isize;
    let mut x = 0u64;
    let mut l = 0usize;

    for (i, c) in sequence.chars().enumerate() {
        let nt = SEQ_NT6_TABLE[c as usize];
        let hit = if (1..=4).contains(&nt) {
            x = ((x << 2) | (nt - 1) as u64) & mask;
            l += 1;
            l >= mlen && mtab.contains(&x)
        } else {
            l = 0;
            x = 0;
            false
        };

        if i >= mlen {
            score += if hit { 1 } else { -penalty };
        }
        if score > max_score {
            max_score = score;
            max_i = i as isize;
        } else if max_score - score > max_drop as i64 {
            break;
        }
    }

    let mut st = 0;
    if max_score >= min_score as i64 {
        sum_telo += ((max_i + 1) as u64) << 32;
        st = (max_i + 1) as usize;
    }

    // Check 3' end
    score = 0;
    max_score = 0;
    max_i = -1;
    x = 0;
    l = 0;

    for (i, c) in sequence.chars().rev().take(slen - st).enumerate() {
        let nt = SEQ_NT6_TABLE[c as usize];
        let hit = if (1..=4).contains(&nt) {
            x = ((x << 2) | (4 - nt) as u64) & mask;
            l += 1;
            l >= mlen && mtab.contains(&x)
        } else {
            l = 0;
            x = 0;
            false
        };

        if i >= mlen {
            score += if hit { 1 } else { -penalty };
        }
        if score > max_score {
            max_score = score;
            max_i = i as isize;
        } else if max_score - score > max_drop as i64 {
            break;
        }
    }

    if max_score >= min_score as i64 {
        sum_telo += (slen - (max_i as usize)) as u64;
    }

    sum_telo
}

pub fn telo_finder<P: AsRef<Path>>(
    fasta_path: P,
    custom_motif: Option<&str>,
    mut output: Option<&mut dyn Write>,
    motif_list: Option<&[String]>,
) -> io::Result<Vec<(bool, bool)>> {
    let sdict = SequenceDict::from_fasta(fasta_path)?;
    let nseq = sdict.sequences.len();
    let mut telo_ends = vec![(false, false); nseq];

    let motifs: Vec<String> = if let Some(motif_list) = motif_list {
        motif_list.to_vec()
    } else if let Some(motif) = custom_motif {
        if !check_motif(motif) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Invalid motif characters",
            ));
        }
        vec![motif.to_string()]
    } else {
        TELO_MOTIF_DB.iter().map(|s| s.to_string()).collect()
    };

    for motif in motifs {
        let mlen = motif.len();
        let mut mtab = HashSet::new();

        // Generate all rotations of the motif
        for i in 0..mlen {
            let mut x = 0u64;
            for j in 0..mlen {
                let c = SEQ_NT6_TABLE[motif.as_bytes()[(i + j) % mlen] as usize];
                assert!((1..=4).contains(&c));
                x = (x << 2) | (c - 1) as u64;
            }
            mtab.insert(x);
        }

        for (i, seq) in sdict.sequences.iter().enumerate() {
            let telo = telo_finder_core(
                &seq.seq,
                &mtab,
                mlen,
                TELO_PENALTY,
                TELO_MAX_DROP,
                TELO_MIN_SCORE,
            );

            if telo == 0 {
                continue;
            }

            if (telo >> 32) > 0 {
                telo_ends[i].0 = true;
                if let Some(out) = output.as_mut() {
                    writeln!(
                        out,
                        "{}\t0\t{}\t{}",
                        seq.name,
                        (telo >> 32) as u32,
                        motif
                    )?;
                } else {
                    eprintln!(
                        "[INFO] found telo motif {} in sequence {} 5'-end up to position {}",
                        motif,
                        seq.name,
                        (telo >> 32) as u32
                    );
                }
            }

            if (telo as u32) > 0 {
                telo_ends[i].1 = true;
                if let Some(out) = output.as_mut() {
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}",
                        seq.name,
                        seq.seq.len() - (telo as u32) as usize,
                        seq.seq.len(),
                        motif
                    )?;
                } else {
                    eprintln!(
                        "[INFO] found telo motif {} in sequence {} 3'-end from position {}",
                        motif,
                        seq.name,
                        seq.seq.len() - (telo as u32) as usize
                    );
                }
            }
        }
    }

    Ok(telo_ends)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_check_motif() {
        assert!(check_motif("AACCCT"));
        assert!(!check_motif("AANCC")); // Contains N which is invalid
    }

    #[test]
    fn test_telo_finder_core() {
        let motif = "AACCCT";
        let mlen = motif.len();
        let mut mtab = HashSet::new();

        // Generate all rotations
        for i in 0..mlen {
            let mut x = 0u64;
            for j in 0..mlen {
                let c = SEQ_NT6_TABLE[motif.as_bytes()[(i + j) % mlen] as usize];
                x = (x << 2) | (c - 1) as u64;
            }
            mtab.insert(x);
        }

        let seq = "AACCCTAACCCTAACCCTGGG";
        let result = telo_finder_core(seq, &mtab, mlen, 1, 2000, 300);
        assert_ne!(result, 0);
    }
}

/// Analyze annotation file and report the most frequent telomere motifs
fn analyze_annotation_file(anno_file: &str) -> Result<Vec<(String, usize)>> {
    let file = File::open(anno_file)?;
    let reader = io::BufReader::new(file);
    let mut motif_counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 4 {
            let motif = parts[3].trim();
            if !motif.is_empty() {
                *motif_counts.entry(motif.to_string()).or_insert(0) += 1;
            }
        }
    }
    
    // Sort by frequency (descending)
    let mut sorted_motifs: Vec<(String, usize)> = motif_counts.into_iter().collect();
    sorted_motifs.sort_by(|a, b| b.1.cmp(&a.1));
    
    Ok(sorted_motifs)
}

/// Extract the most frequent motif from annotation file
fn extract_primary_motif_from_anno(anno_file: &str) -> Result<String> {
    let file = File::open(anno_file)?;
    let reader = io::BufReader::new(file);
    let mut motif_counts = std::collections::HashMap::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 4 {
            let motif = parts[3].to_string();
            *motif_counts.entry(motif).or_insert(0) += 1;
        }
    }
    
    if motif_counts.is_empty() {
        return Err(anyhow::anyhow!("No motifs found in annotation file"));
    }
    
    // Find most frequent motif
    let primary_motif = motif_counts
        .iter()
        .max_by_key(|(_, count)| *count)
        .map(|(motif, _)| motif.clone())
        .ok_or_else(|| anyhow::anyhow!("Could not determine primary motif"))?;
    
    println!("Motif frequency analysis from {}:", anno_file);
    let mut sorted_motifs: Vec<_> = motif_counts.iter().collect();
    sorted_motifs.sort_by(|a, b| b.1.cmp(a.1));
    
    for (motif, count) in sorted_motifs.iter().take(5) {
        let percentage = (**count as f64 / motif_counts.values().sum::<usize>() as f64) * 100.0;
        println!("  {}: {} occurrences ({:.1}%)", motif, count, percentage);
    }
    
    Ok(primary_motif)
}

/// Report the most frequent telomere motifs from annotation files
fn report_most_frequent_motifs(initial_anno_file: &str, final_anno_file: &str) -> Result<()> {
    println!("\n=== TELOMERE MOTIF FREQUENCY ANALYSIS ===");
    
    // Check which annotation file to analyze
    let (anno_file, source) = if std::path::Path::new(initial_anno_file).exists() {
        (initial_anno_file, "predefined database")
    } else if std::path::Path::new(final_anno_file).exists() {
        (final_anno_file, "discovered database")
    } else {
        println!("No annotation files found. No telomere motifs detected.");
        // Write empty result to telomere_candidate.txt
        let candidate = TelomereCandidate {
            status: "no_motifs".to_string(),
            primary_motif: None,
            frequency: None,
            percentage: None,
            source: None,
            total_occurrences: Some(0),
            unique_motifs: Some(0),
            top_motifs: vec![],
            message: Some("No telomere motifs detected in the genome.".to_string()),
        };
        let json_output = serde_json::to_string_pretty(&candidate)?;
        std::fs::write("telomere_candidate.json", json_output)?;
        return Ok(());
    };
    
    println!("Analyzing annotation file: {} (from {})", anno_file, source);
    
    let motif_frequencies = analyze_annotation_file(anno_file)?;
    
    if motif_frequencies.is_empty() {
        println!("No telomere motifs found in the annotation file.");
        // Write empty result to telomere_candidate.txt
        let candidate = TelomereCandidate {
            status: "no_motifs".to_string(),
            primary_motif: None,
            frequency: None,
            percentage: None,
            source: Some(source.to_string()),
            total_occurrences: Some(0),
            unique_motifs: Some(0),
            top_motifs: vec![],
            message: Some("No telomere motifs found in the annotation file.".to_string()),
        };
        let json_output = serde_json::to_string_pretty(&candidate)?;
        std::fs::write("telomere_candidate.json", json_output)?;
        return Ok(());
    }
    
    println!("\nMost frequent telomere motifs:");
    println!("{:<20} {:<10} {:<10}", "Motif", "Frequency", "Percentage");
    println!("{}", "-".repeat(40));
    
    let total_occurrences: usize = motif_frequencies.iter().map(|(_, count)| count).sum();
    
    for (motif, count) in motif_frequencies.iter().take(10) {
        let percentage = (*count as f64 / total_occurrences as f64) * 100.0;
        println!("{:<20} {:<10} {:<10.1}%", motif, count, percentage);
    }
    
    if motif_frequencies.len() > 10 {
        println!("... and {} more motifs", motif_frequencies.len() - 10);
    }
    
    // Report the top motif as the most likely telomere motif
    if let Some((top_motif, top_count)) = motif_frequencies.first() {
        let percentage = (*top_count as f64 / total_occurrences as f64) * 100.0;
        println!("\n=== PRIMARY TELOMERE MOTIF ===");
        println!("Most frequent motif: {}", top_motif);
        println!("Occurrences: {} ({:.1}% of all detected motifs)", top_count, percentage);
        println!("Source: {}", source);
        
        // Create top motifs list
        let top_motifs: Vec<MotifFrequency> = motif_frequencies.iter().take(10).map(|(motif, count)| {
            MotifFrequency {
                motif: motif.clone(),
                frequency: *count,
                percentage: (*count as f64 / total_occurrences as f64) * 100.0,
            }
        }).collect();
        
        // Write the primary telomere motif to telomere_candidate.txt as JSON
        let candidate = TelomereCandidate {
            status: "success".to_string(),
            primary_motif: Some(top_motif.clone()),
            frequency: Some(*top_count),
            percentage: Some(percentage),
            source: Some(source.to_string()),
            total_occurrences: Some(total_occurrences),
            unique_motifs: Some(motif_frequencies.len()),
            top_motifs,
            message: None,
        };
        
        let json_output = serde_json::to_string_pretty(&candidate)?;
        std::fs::write("telomere_candidate.json", json_output)?;
        
        println!("Primary telomere motif written to: telomere_candidate.json");
    }
    
    // Save detailed results to a file
    let output_file = "telomere_motif_summary.txt";
    let mut summary_file = File::create(output_file)?;
    writeln!(summary_file, "TELOMERE MOTIF FREQUENCY SUMMARY")?;
    writeln!(summary_file, "=================================")?;
    writeln!(summary_file, "Source: {}", source)?;
    writeln!(summary_file, "Annotation file: {}", anno_file)?;
    writeln!(summary_file, "Total motif occurrences: {}", total_occurrences)?;
    writeln!(summary_file, "Unique motifs: {}", motif_frequencies.len())?;
    writeln!(summary_file)?;
    writeln!(summary_file, "{:<20} {:<10} {:<10}", "Motif", "Frequency", "Percentage")?;
    writeln!(summary_file, "{}", "-".repeat(40))?;
    
    for (motif, count) in &motif_frequencies {
        let percentage = (*count as f64 / total_occurrences as f64) * 100.0;
        writeln!(summary_file, "{:<20} {:<10} {:<10.1}%", motif, count, percentage)?;
    }
    
    println!("Detailed summary saved to: {}", output_file);
    
    Ok(())
}

fn parse_args(args: &[String]) -> (String, usize, usize, bool, bool, f64, bool) {
    let mut fasta_file = String::new();
    let mut max_indels = 5;
    let mut max_gap_size = 5;
    let mut strict_mode = false;
    let mut annotate_genome = false;
    let mut annotation_threshold = 0.5;
    let mut use_window_density = false;
    
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--max-indels" => {
                if i + 1 < args.len() {
                    max_indels = args[i + 1].parse().unwrap_or(5);
                    i += 2;
                } else {
                    eprintln!("Error: --max-indels requires a value");
                    std::process::exit(1);
                }
            }
            "--max-gap-size" => {
                if i + 1 < args.len() {
                    max_gap_size = args[i + 1].parse().unwrap_or(5);
                    i += 2;
                } else {
                    eprintln!("Error: --max-gap-size requires a value");
                    std::process::exit(1);
                }
            }
            "--strict" => {
                strict_mode = true;
                i += 1;
            }
            "--annotate" => {
                annotate_genome = true;
                i += 1;
            }
            "--window-density" => {
                use_window_density = true;
                i += 1;
            }
            "--annotation-threshold" => {
                if i + 1 < args.len() {
                    annotation_threshold = args[i + 1].parse().unwrap_or(0.5);
                    i += 2;
                } else {
                    eprintln!("Error: --annotation-threshold requires a value");
                    std::process::exit(1);
                }
            }
            _ => {
                if !args[i].starts_with("--") && fasta_file.is_empty() {
                    fasta_file = args[i].clone();
                }
                i += 1;
            }
        }
    }
    
    (fasta_file, max_indels, max_gap_size, strict_mode, annotate_genome, annotation_threshold, use_window_density)
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <fasta_file> [options]", args[0]);
        eprintln!("       {} kmc <input_fasta> <k> <db_prefix> <output_txt>", args[0]);
        eprintln!("       {} <fasta_file> extract_lastN <N> <output_fasta>", args[0]);
        eprintln!();
        eprintln!("Note: Supports both .fasta and .fasta.gz files");
        eprintln!();
        eprintln!("Options:");
        eprintln!("  --max-indels <N>      Maximum indels allowed in stretch analysis (default: 5)");
        eprintln!("  --max-gap-size <N>    Maximum gap size before resetting stretch (default: 5)");
        eprintln!("  --strict              Use exact matching only (no indel tolerance)");
        eprintln!("  --annotate            Generate genome annotation with BED file output");
        eprintln!("  --annotation-threshold <F>  Score threshold for annotations (default: 0.5)");
        eprintln!("  --window-density      Use 1000bp window density instead of longest stretch");
        eprintln!();
        eprintln!("Examples:");
        eprintln!("  {} genome.fasta                           # Standard analysis with indel tolerance", args[0]);
        eprintln!("  {} genome.fasta --max-indels 3 --max-gap-size 3  # More strict indel tolerance", args[0]);
        eprintln!("  {} genome.fasta --strict                  # Exact matching only", args[0]);
        eprintln!("  {} genome.fasta --annotate                # Generate BED file annotations", args[0]);
        eprintln!("  {} genome.fasta.gz --window-density       # Use window-based k-mer density (gzipped)", args[0]);
        eprintln!("  {} genome.fasta --annotate --annotation-threshold 0.7  # Higher threshold", args[0]);
        eprintln!("  {} genome.fasta extract_lastN 10000 last10000.fasta", args[0]);
        eprintln!("  {} kmc last5000.fasta 7 kmc_db kmc_dump.txt", args[0]);
        std::process::exit(1);
    }

    // KMC canonicalization pipeline
    if args.len() == 6 && args[1] == "kmc" {
        let input_fasta = &args[2];
        let k: usize = args[3].parse().expect("k must be an integer");
        let db_prefix = &args[4];
        let output_txt = &args[5];
        // Run KMC and dump
        kmers::run_kmc(input_fasta, k, db_prefix)?;
        kmers::run_kmc_dump(db_prefix, output_txt)?;
        // Read and canonicalize
        let kmc_counts = kmers::read_kmc_counts(output_txt)?;
        let analyses = kmers::analyze_kmc_kmers(kmc_counts, k);
        // Save canonicalized k-mer table
        let mut file = std::fs::File::create("canonical_kmers.tsv")?;
        writeln!(file, "Kmer\tForward\tRC\tTotal")?;
        for a in &analyses {
            writeln!(file, "{}\t{}\t{}\t{}", a.kmer, a.forward_count, a.rc_count, a.total_count)?;
        }
        println!("Saved canonicalized k-mer table to canonical_kmers.tsv ({} k-mers)", analyses.len());
        return Ok(());
    }

    // Handle special modes first
    if args.contains(&"extract_lastN".to_string()) {
        if let Some(pos) = args.iter().position(|x| x == "extract_lastN") {
            if pos + 2 < args.len() {
                let fasta_path = &args[1];
                let n: usize = args[pos + 1].parse().expect("N must be an integer");
                let output_fasta = &args[pos + 2];
                extract_last_n_bp_to_fasta(fasta_path, output_fasta, n)?;
                println!("Extracted last {} bp of each scaffold ≥1MB to {}", n, output_fasta);
                return Ok(());
            }
        }
    }

    // Parse command line arguments
    let (fasta_path, max_indels, max_gap_size, strict_mode, annotate_genome, annotation_threshold, use_window_density) = parse_args(&args);
    
    if fasta_path.is_empty() {
        eprintln!("Error: FASTA file not specified");
        std::process::exit(1);
    }

    // Default analysis: extract last 5000bp and run telomere analysis
    const DEFAULT_BP_SIZE: usize = 5000;
    
    println!("=== TeloX Telomere Analysis ===");
    println!("Input: {}", fasta_path);
    println!("Analyzing last {} bp of scaffolds ≥1MB", DEFAULT_BP_SIZE);
    if strict_mode {
        println!("Mode: Exact matching (no indel tolerance)");
    } else {
        println!("Mode: Indel-tolerant (max_indels: {}, max_gap_size: {})", max_indels, max_gap_size);
    }
    println!();
    
    // Step 1: First run telo_finder on TELO_MOTIF_DB
    println!("Step 1: Running telo_finder with predefined TELO_MOTIF_DB...");
    let mut initial_anno_file = std::fs::File::create("initial_anno.txt")?;
    let initial_telo_results = telo_finder(&fasta_path, None, Some(&mut initial_anno_file), None)?;
    
    // Check if any telomere motifs were found
    let found_any_telomeres = initial_telo_results.iter().any(|(five_prime, three_prime)| *five_prime || *three_prime);
    
    if found_any_telomeres {
        println!("Found telomere motifs with predefined database. Results written to initial_anno.txt");
        
        // Create JSON output for telo_finder results
        // Extract primary motif and calculate statistics
        let primary_motif = match extract_primary_motif_from_anno("initial_anno.txt") {
            Ok(motif) => motif,
            Err(e) => {
                eprintln!("Warning: Could not extract primary motif from telo_finder results: {}", e);
                eprintln!("Attempting to create JSON with first available motif...");
                // Fallback: try to get any motif from the file
                let mut fallback_motif: Option<String> = None;
                if let Ok(contents) = std::fs::read_to_string("initial_anno.txt") {
                    for line in contents.lines() {
                        if !line.trim().is_empty() && !line.starts_with('#') {
                            let parts: Vec<&str> = line.split('\t').collect();
                            if parts.len() >= 4 {
                                fallback_motif = Some(parts[3].to_string());
                                break;
                            }
                        }
                    }
                }
                fallback_motif.ok_or_else(|| anyhow::anyhow!("Could not determine primary motif from telo_finder results"))?
            }
        };
        
        println!("Creating JSON output for telo_finder results...");
        
        // Parse initial_anno.txt to count occurrences and calculate longest stretch
        let mut total_count = 0;
        let mut forward_count = 0;
        let mut reverse_count = 0;
        let mut longest_stretch = 0usize;
        
        if let Ok(contents) = std::fs::read_to_string("initial_anno.txt") {
            for line in contents.lines() {
                if !line.trim().is_empty() && !line.starts_with('#') {
                    let parts: Vec<&str> = line.split('\t').collect();
                    if parts.len() >= 4 && parts[3] == primary_motif {
                        total_count += 1;
                        
                        // Determine strand and calculate stretch from positions
                        // Format: chromosome/seq_name, start, end, motif
                        // telo_finder writes: 5' end as (seq_name, 0, end_pos, motif) and 3' end as (seq_name, start_pos, seq_len, motif)
                        if parts.len() >= 4 {
                            if let (Ok(start), Ok(end)) = (parts[1].parse::<usize>(), parts[2].parse::<usize>()) {
                                let stretch = end - start;
                                if stretch > longest_stretch {
                                    longest_stretch = stretch;
                                }
                                
                                // Infer strand from position:
                                // - start == 0 means 5' end (forward strand, +)
                                // - start > 0 means 3' end (reverse strand, -)
                                if start == 0 {
                                    forward_count += 1;
                                } else {
                                    reverse_count += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Convert longest_stretch from bp to number of motif repeats (approximate)
        let motif_len = primary_motif.len();
        let longest_stretch_repeats = if motif_len > 0 {
            longest_stretch / motif_len
        } else {
            0
        };
        
        // Create JSON structure matching Stage 2 format
        let telomere_result = serde_json::json!({
            "canonical_motif": primary_motif,
            "rotational_variants": [primary_motif.clone()],
            "forward_count": forward_count,
            "reverse_count": reverse_count,
            "total_count": total_count,
            "longest_stretch": longest_stretch_repeats,
            "longest_stretch_bp": longest_stretch,
            "data_source": "telo_finder_predefined_database"
        });
        
        // Write to JSON file (always create, even if some stats are missing)
        match std::fs::write("telomere_motif_final.json", serde_json::to_string_pretty(&telomere_result).unwrap()) {
            Ok(_) => {
                println!("✓ Primary telomere motif identified by telo_finder:");
                println!("  Canonical sequence: {}", primary_motif);
                println!("  Total occurrences: {} (Forward: {}, Reverse: {})", total_count, forward_count, reverse_count);
                println!("  Longest stretch: {} bp (approximately {} repeats)", longest_stretch, longest_stretch_repeats);
                println!("  Result written to: telomere_motif_final.json");
            },
            Err(e) => {
                eprintln!("❌ Error writing telomere_motif_final.json: {}", e);
                return Err(anyhow::anyhow!("Failed to write telomere_motif_final.json: {}", e));
            }
        }
        
        // Step 4: Genome annotation (if requested)
        if annotate_genome {
            println!("\nStep 4: Generating genome annotations with identified telomere motifs...");
            
            // Extract the most frequent motif from initial_anno.txt
            match extract_primary_motif_from_anno("initial_anno.txt") {
                Ok(primary_motif) => {
                    println!("Primary motif identified from telo_finder results: {}", primary_motif);
                    
                    let primary_motif_vec = vec![primary_motif.clone()];
                    let annotation_config = create_annotation_config(
                        Some(annotation_threshold),
                        Some(max_gap_size * 10), // Use larger gap for annotation merging
                        Some(true), // Enable merging
                    );
                    
                    let output_prefix = "telox_genome";
                    match annotate_genome_with_motifs(&fasta_path, &primary_motif_vec, output_prefix, Some(annotation_config)) {
                        Ok(_annotations) => {
                            println!("Genome annotation completed successfully with motif: {}", primary_motif);
                        },
                        Err(e) => eprintln!("Warning: Genome annotation failed: {}", e),
                    }
                },
                Err(e) => {
                    eprintln!("Warning: Could not extract primary motif from initial_anno.txt: {}", e);
                    println!("Skipping genome annotation.");
                }
            }
        }
        
        println!("Analysis complete.");
        return Ok(());
    }
    
    println!("No telomere motifs found with predefined database. Proceeding to k-mer analysis...");
    
    // Step 2: If no telomeres found, run k-mer analysis to generate new motifs
    println!("Step 2: Running k-mer analysis to discover potential telomere motifs...");
    
    for k in 5..=15 {
        println!("Processing {}-mers...", k);
        
        // Use optimized k-mer counting for last N bp of each scaffold ≥1MB
        let counts = kmers::count_kmers_last_n_bp_parallel(&fasta_path, k, DEFAULT_BP_SIZE)
            .with_context(|| format!("Failed to count {}-mers in last {}bp of {}", k, DEFAULT_BP_SIZE, &fasta_path))?;

        // Print k-mer count table for debugging
        println!("Top 20 {}-mers by count:", k);
        kmers::print_kmer_table(&counts, 20);

        // Calculate longest stretch for each k-mer using only the last N bp of each scaffold ≥1MB
        let sdict = SequenceDict::from_fasta(&fasta_path)?;
        let kmer_list: Vec<String> = counts.keys().cloned().collect();
        let mut longest_stretch_map: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
        let mut longest_stretch_indels_map: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
        let mut processed_scaffolds = 0;
        let mut filtered_scaffolds = 0;
        
        // Use configured indel tolerance parameters
        
        for seq in &sdict.sequences {
            // Only process scaffolds larger than 1MB (1,000,000 bp)
            if seq.seq.len() < 1_000_000 {
                filtered_scaffolds += 1;
                continue;
            }
            
            processed_scaffolds += 1;
            let region = if seq.seq.len() <= DEFAULT_BP_SIZE {
                &seq.seq[..]
            } else {
                &seq.seq[seq.seq.len() - DEFAULT_BP_SIZE..]
            };
            
            // Calculate stretch using selected method
            let stretch_map = if use_window_density {
                // Use 1000bp window density approach
                kmers::calculate_kmer_window_density(region, &kmer_list, 1000)
            } else {
                // Use traditional longest stretch
                kmers::longest_continuous_stretch_for_kmers(region, &kmer_list)
            };
            
            // Calculate indel-tolerant stretch (only if not in strict mode and not using window density)
            let stretch_indels_map = if strict_mode || use_window_density {
                stretch_map.clone()  // Use exact stretch for both
            } else {
                kmers::longest_continuous_stretch_with_indels(region, &kmer_list, max_indels, max_gap_size)
            };
            
            for (kmer, stretch) in stretch_map {
                let entry = longest_stretch_map.entry(kmer).or_insert(0);
                if stretch > *entry {
                    *entry = stretch;
                }
            }
            
            for (kmer, stretch) in stretch_indels_map {
                let entry = longest_stretch_indels_map.entry(kmer).or_insert(0);
                if stretch > *entry {
                    *entry = stretch;
                }
            }
        }
        
        if k == 6 { // Only print this once to avoid spam
            eprintln!("Longest stretch calculation: processed {} scaffolds >= 1MB, filtered out {} smaller scaffolds", processed_scaffolds, filtered_scaffolds);
            eprintln!("Using last {} bp of each large scaffold for analysis", DEFAULT_BP_SIZE);
        }
        
        // Print top longest stretch results for debugging
        let mut stretch_debug: Vec<_> = longest_stretch_map.iter().collect();
        stretch_debug.sort_by(|a, b| b.1.cmp(a.1));
        let method_name = if use_window_density { "window density (1000bp)" } else { "longest stretch" };
        println!("Top 10 {}-mers by {} method:", k, method_name);
        for (i, (kmer, stretch)) in stretch_debug.iter().take(10).enumerate() {
            let count = counts.get(*kmer).map(|p| p.forward + p.rc).unwrap_or(0);
            println!("  {}: {} (stretch: {}, count: {})", i+1, kmer, stretch, count);
        }

        // Strand bias analysis with both exact and indel-tolerant stretch metrics
        let bias_analyses = kmers::analyze_strand_bias_with_indels(&counts, Some(&longest_stretch_map), Some(&longest_stretch_indels_map));
        
        // Save strand bias results with new indel-tolerant metrics
        let bias_output_filename = format!("strand_bias_{}mer.tsv", k);
        kmers::save_strand_bias_table(&bias_analyses, &bias_output_filename)?;
        println!("Saved strand bias results to {}", bias_output_filename);
        
        // Print summary of improvements from indel tolerance
        let improved_count = bias_analyses.iter().filter(|a| a.indel_tolerance_used).count();
        if improved_count > 0 {
            println!("  -> {} k-mers benefited from indel tolerance ({}%)", 
                     improved_count, 
                     (improved_count * 100) / bias_analyses.len());
        }
    }

    // Instead of reading filtered files, gather all filtered analyses in memory and rank them
    println!("Step 2 complete. Now ranking and consolidating k-mers...");
    let mut all_filtered_analyses = Vec::new();
    for k in 5..=15 {
        let bias_output_filename = format!("strand_bias_{}mer.tsv", k);
        println!("Reading {}...", bias_output_filename);
        
        let file = match std::fs::File::open(&bias_output_filename) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Warning: Could not open {}: {}", bias_output_filename, e);
                continue;
            }
        };
        
        let reader = std::io::BufReader::new(file);
        let mut line_count = 0;
        let mut added_count = 0;
        
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            line_count += 1;
            if i == 0 || line.trim().is_empty() { continue; }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 10 { 
                eprintln!("Warning: Line {} in {} has only {} columns, expected 10", i+1, bias_output_filename, cols.len());
                continue; 
            }
            let total: u32 = cols[3].parse().unwrap_or(0);
            let longest_stretch_exact: usize = cols[7].parse().unwrap_or(0);
            let longest_stretch_indels: usize = cols[8].parse().unwrap_or(0);
            let significance = cols[6];
            
            // Use indel-tolerant stretch for filtering and ranking
            let ranking_stretch = longest_stretch_indels;
            
            if ranking_stretch <= 3 || significance == "weak" || total <= 10 { continue; }
            all_filtered_analyses.push((
                cols[0].to_string(), // kmer
                cols[1].parse().unwrap_or(0), // forward
                cols[2].parse().unwrap_or(0), // rc
                cols[3].parse().unwrap_or(0), // total
                cols[4].to_string(), // bias_ratio
                cols[5].to_string(), // direction
                cols[6].to_string(), // significance
                ranking_stretch,     // use indel-tolerant stretch for ranking
            ));
            added_count += 1;
        }
        println!("  Read {} lines, added {} filtered k-mers", line_count, added_count);
    }
    // Sort by longest_stretch DESC, then total DESC
    println!("Sorting {} total k-mers...", all_filtered_analyses.len());
    all_filtered_analyses.sort_by(|a, b| b.7.cmp(&a.7).then(b.3.cmp(&a.3)));
    
    // Write to rank.tsv
    println!("Writing ranked results to rank.tsv...");
    let mut out = std::fs::File::create("rank.tsv")?;
    let header_method = if use_window_density { "WindowDensity1000bp" } else { "LongestStretchIndels" };
    writeln!(out, "Kmer\tForward\tRC\tTotal\tBiasRatio\tDirection\tSignificance\t{}", header_method)?;
    for k in all_filtered_analyses {
        writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", k.0, k.1, k.2, k.3, k.4, k.5, k.6, k.7)?;
    }
    println!("Ranked k-mers written to rank.tsv");

    // Step 3: Generate new motif array from k-mer analysis and run telo_finder
    println!("Step 3: Generating new telomere motif database from k-mer analysis...");
    let discovered_motifs = kmers::consolidate_ranked_motifs_by_rotation("rank.tsv")?;
    println!("Running telo_finder with {} discovered motifs...", discovered_motifs.len());
    
    // Run telo_finder with the new motifs array and print results to anno.txt
    let mut anno_file = std::fs::File::create("anno.txt")?;
    let final_telo_results = telo_finder(&fasta_path, None, Some(&mut anno_file), Some(&discovered_motifs))?;
    
    // Check if any telomere motifs were found with the new database
    let found_any_telomeres_final = final_telo_results.iter().any(|(five_prime, three_prime)| *five_prime || *three_prime);
    
    if found_any_telomeres_final {
        println!("Found telomere motifs with discovered database. Results written to anno.txt");
        
        // Create JSON output for k-mer discovery telo_finder results
        match extract_primary_motif_from_anno("anno.txt") {
            Ok(primary_motif) => {
                println!("Creating JSON output for k-mer discovery telo_finder results...");
                
                // Count total occurrences of the primary motif
                let mut total_count = 0;
                let mut forward_count = 0;
                let mut reverse_count = 0;
                
                if let Ok(contents) = std::fs::read_to_string("anno.txt") {
                    for line in contents.lines() {
                        if !line.trim().is_empty() && !line.starts_with('#') {
                            let parts: Vec<&str> = line.split('\t').collect();
                            if parts.len() >= 4 && parts[3] == primary_motif {
                                total_count += 1;
                                if parts.len() >= 5 && parts[4] == "+" {
                                    forward_count += 1;
                                } else if parts.len() >= 5 && parts[4] == "-" {
                                    reverse_count += 1;
                                }
                            }
                        }
                    }
                }
                
                // Get longest stretch from the consolidated ranking if available
                let mut longest_stretch = 0;
                if let Ok(contents) = std::fs::read_to_string("rank.tsv") {
                    for line in contents.lines() {
                        if !line.trim().is_empty() && !line.starts_with('#') {
                            let parts: Vec<&str> = line.split('\t').collect();
                            if parts.len() >= 2 && parts[0] == primary_motif {
                                if let Ok(stretch) = parts[1].parse::<usize>() {
                                    longest_stretch = stretch;
                                }
                                break;
                            }
                        }
                    }
                }
                
                // Create JSON structure
                let telomere_result = serde_json::json!({
                    "canonical_motif": primary_motif,
                    "rotational_variants": [primary_motif],
                    "forward_count": forward_count,
                    "reverse_count": reverse_count,
                    "total_count": total_count,
                    "longest_stretch": longest_stretch,
                    "data_source": "telo_finder_kmer_discovery"
                });
                
                // Write to JSON file (overwrite if exists)
                match std::fs::write("telomere_motif_final.json", serde_json::to_string_pretty(&telomere_result).unwrap()) {
                    Ok(_) => {
                        println!("Primary telomere motif identified by k-mer discovery telo_finder:");
                        println!("Canonical sequence: {}", primary_motif);
                        println!("Total occurrences: {} (Forward: {}, Reverse: {})", total_count, forward_count, reverse_count);
                        println!("Longest stretch: {}", longest_stretch);
                        println!("Result written to: telomere_motif_final.json");
                    },
                    Err(e) => {
                        eprintln!("Error writing telomere_motif_final.json: {}", e);
                    }
                }
            },
            Err(e) => {
                eprintln!("Warning: Could not extract primary motif from k-mer discovery telo_finder results: {}", e);
            }
        }
    } else {
        println!("No telomere motifs found even with discovered database.");
    }
    
    println!("Analysis complete.");

    // Step 4: Genome annotation (if requested)
    if annotate_genome {
        println!("\nStep 4: Generating genome annotations...");
        
        // Use only the most likely telomere motif (first in consolidated ranking)
        if let Some(primary_motif) = discovered_motifs.first() {
            println!("Annotating genome with primary telomere motif: {}", primary_motif);
            
            let primary_motif_vec = vec![primary_motif.clone()];
            let annotation_config = create_annotation_config(
                Some(annotation_threshold),
                Some(max_gap_size * 10), // Use larger gap for annotation merging
                Some(true), // Enable merging
            );
            
            let output_prefix = "telox_genome";
            match annotate_genome_with_motifs(&fasta_path, &primary_motif_vec, output_prefix, Some(annotation_config)) {
                Ok(_annotations) => {
                    println!("Genome annotation completed successfully with motif: {}", primary_motif);
                },
                Err(e) => eprintln!("Warning: Genome annotation failed: {}", e),
            }
        } else {
            println!("No telomere motifs available for annotation.");
        }
    }

    Ok(())
}
