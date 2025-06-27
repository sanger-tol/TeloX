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
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::env;
use anyhow::{Context, Result};
mod kmers;
use kmers::{count_kmers_in_fasta, print_kmer_table, analyze_strand_bias, print_strand_bias_table, save_strand_bias_table, get_strand_bias_summary, longest_continuous_stretch_for_kmers, filter_strand_bias_tsvs, consolidate_rotational_kmers, extract_last_n_bp_to_fasta};

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

impl SequenceDict {
    fn from_fasta<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        let reader = io::BufReader::new(file);
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
        let mut candidate_file = File::create("telomere_candidate.txt")?;
        writeln!(candidate_file, "No telomere motifs detected in the genome.")?;
        return Ok(());
    };
    
    println!("Analyzing annotation file: {} (from {})", anno_file, source);
    
    let motif_frequencies = analyze_annotation_file(anno_file)?;
    
    if motif_frequencies.is_empty() {
        println!("No telomere motifs found in the annotation file.");
        // Write empty result to telomere_candidate.txt
        let mut candidate_file = File::create("telomere_candidate.txt")?;
        writeln!(candidate_file, "No telomere motifs found in the annotation file.")?;
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
        
        // Write the primary telomere motif to telomere_candidate.txt
        let mut candidate_file = File::create("telomere_candidate.txt")?;
        writeln!(candidate_file, "PRIMARY TELOMERE MOTIF CANDIDATE")?;
        writeln!(candidate_file, "=================================")?;
        writeln!(candidate_file, "Motif: {}", top_motif)?;
        writeln!(candidate_file, "Frequency: {} occurrences", top_count)?;
        writeln!(candidate_file, "Percentage: {:.1}% of all detected motifs", percentage)?;
        writeln!(candidate_file, "Source: {}", source)?;
        writeln!(candidate_file, "Total motif occurrences: {}", total_occurrences)?;
        writeln!(candidate_file, "Unique motifs detected: {}", motif_frequencies.len())?;
        writeln!(candidate_file)?;
        writeln!(candidate_file, "Top 10 most frequent motifs:")?;
        writeln!(candidate_file, "{:<20} {:<10} {:<10}", "Motif", "Frequency", "Percentage")?;
        writeln!(candidate_file, "{}", "-".repeat(40))?;
        
        for (motif, count) in motif_frequencies.iter().take(10) {
            let percentage = (*count as f64 / total_occurrences as f64) * 100.0;
            writeln!(candidate_file, "{:<20} {:<10} {:<10.1}%", motif, count, percentage)?;
        }
        
        println!("Primary telomere motif written to: telomere_candidate.txt");
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

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <fasta_file> [extract_lastN <N> <output_fasta>]", args[0]);
        eprintln!("       {} kmc <input_fasta> <k> <db_prefix> <output_txt>", args[0]);
        eprintln!("Example: {} genome.fasta", args[0]);
        eprintln!("         {} genome.fasta extract_lastN 5000 last5000.fasta", args[0]);
        eprintln!("         {} kmc last5000.fasta 7 kmc_db kmc_dump.txt", args[0]);
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

    let fasta_path = args[1].clone();

    // Optionally extract last N bp to a new FASTA file
    if args.len() == 5 && args[2] == "extract_lastN" {
        let n: usize = args[3].parse().expect("N must be an integer");
        let output_fasta = &args[4];
        extract_last_n_bp_to_fasta(&fasta_path, output_fasta, n)?;
        println!("Extracted last {} bp of each scaffold to {}", n, output_fasta);
        return Ok(());
    }

    // Step 1: First run telo_finder on TELO_MOTIF_DB
    println!("Step 1: Running telo_finder with predefined TELO_MOTIF_DB...");
    let mut initial_anno_file = std::fs::File::create("initial_anno.txt")?;
    let initial_telo_results = telo_finder(&fasta_path, None, Some(&mut initial_anno_file), None)?;
    
    // Check if any telomere motifs were found
    let found_any_telomeres = initial_telo_results.iter().any(|(five_prime, three_prime)| *five_prime || *three_prime);
    
    if found_any_telomeres {
        println!("Found telomere motifs with predefined database. Results written to initial_anno.txt");
        println!("Analysis complete.");
        return Ok(());
    }
    
    println!("No telomere motifs found with predefined database. Proceeding to k-mer analysis...");
    
    // Step 2: If no telomeres found, run k-mer analysis to generate new motifs
    println!("Step 2: Running k-mer analysis to discover potential telomere motifs...");
    
    for k in 5..=12 {
        println!("Processing {}-mers...", k);
        
        // Use optimized k-mer counting for last 5000bp of each scaffold
        let counts = kmers::count_kmers_last_5000bp_parallel(&fasta_path, k)
            .with_context(|| format!("Failed to count {}-mers in last 5000bp of {}", k, &fasta_path))?;

        // Calculate longest stretch for each k-mer using only the last 5000bp of each scaffold
        let sdict = SequenceDict::from_fasta(&fasta_path)?;
        let kmer_list: Vec<String> = counts.keys().cloned().collect();
        let mut longest_stretch_map: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
        for seq in &sdict.sequences {
            let region = if seq.seq.len() <= 5000 {
                &seq.seq[..]
            } else {
                &seq.seq[seq.seq.len() - 5000..]
            };
            let stretch_map = kmers::longest_continuous_stretch_for_kmers(region, &kmer_list);
            for (kmer, stretch) in stretch_map {
                let entry = longest_stretch_map.entry(kmer).or_insert(0);
                if stretch > *entry {
                    *entry = stretch;
                }
            }
        }
        
        // Strand bias analysis
        let mut bias_analyses = kmers::analyze_strand_bias(&counts, Some(&longest_stretch_map));
        
        // Filter out k-mers with longest stretch < 2 and weak significance
        let filtered_analyses: Vec<&kmers::StrandBiasAnalysis> = bias_analyses.iter()
            .filter(|analysis| {
                analysis.longest_stretch >= 2 && analysis.significance != "weak"
            })
            .collect();
        
        // Save strand bias results
        let bias_output_filename = format!("strand_bias_{}mer.tsv", k);
        let mut content = String::new();
        content.push_str("Kmer\tForward\tRC\tTotal\tBiasRatio\tDirection\tSignificance\tLongestStretch\n");
        for analysis in &bias_analyses {
            let ratio_str = if analysis.bias_ratio == f64::INFINITY {
                "inf".to_string()
            } else {
                format!("{:.3}", analysis.bias_ratio)
            };
            content.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                analysis.kmer,
                analysis.forward_count,
                analysis.rc_count,
                analysis.total_count,
                ratio_str,
                analysis.bias_direction,
                analysis.significance,
                analysis.longest_stretch
            ));
        }
        std::fs::write(&bias_output_filename, content)?;
        println!("Saved strand bias results to {}", bias_output_filename);
    }

    // Instead of reading filtered files, gather all filtered analyses in memory and rank them
    let mut all_filtered_analyses = Vec::new();
    for k in 5..=12 {
        let bias_output_filename = format!("strand_bias_{}mer.tsv", k);
        let file = std::fs::File::open(&bias_output_filename)?;
        let reader = std::io::BufReader::new(file);
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if i == 0 || line.trim().is_empty() { continue; }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 8 { continue; }
            let longest_stretch: usize = cols[7].parse().unwrap_or(0);
            let significance = cols[6];
            if longest_stretch < 2 || significance == "weak" { continue; }
            all_filtered_analyses.push((
                cols[0].to_string(), // kmer
                cols[1].parse().unwrap_or(0), // forward
                cols[2].parse().unwrap_or(0), // rc
                cols[3].parse().unwrap_or(0), // total
                cols[4].to_string(), // bias_ratio
                cols[5].to_string(), // direction
                cols[6].to_string(), // significance
                longest_stretch,
            ));
        }
    }
    // Sort by longest_stretch DESC, then total DESC
    all_filtered_analyses.sort_by(|a, b| b.7.cmp(&a.7).then(b.3.cmp(&a.3)));
    // Write to rank.tsv
    let mut out = std::fs::File::create("rank.tsv")?;
    writeln!(out, "Kmer\tForward\tRC\tTotal\tBiasRatio\tDirection\tSignificance\tLongestStretch")?;
    for k in all_filtered_analyses {
        writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", k.0, k.1, k.2, k.3, k.4, k.5, k.6, k.7)?;
    }
    println!("Ranked k-mers written to rank.tsv");

    // Step 3: Generate new motif array from k-mer analysis and run telo_finder
    println!("Step 3: Generating new telomere motif database from k-mer analysis...");
    let motifs = kmers::consolidate_ranked_motifs_by_rotation("rank.tsv")?;
    println!("Running telo_finder with {} discovered motifs...", motifs.len());
    
    // Run telo_finder with the new motifs array and print results to anno.txt
    let mut anno_file = std::fs::File::create("anno.txt")?;
    let final_telo_results = telo_finder(&fasta_path, None, Some(&mut anno_file), Some(&motifs))?;
    
    // Check if any telomere motifs were found with the new database
    let found_any_telomeres_final = final_telo_results.iter().any(|(five_prime, three_prime)| *five_prime || *three_prime);
    
    if found_any_telomeres_final {
        println!("Found telomere motifs with discovered database. Results written to anno.txt");
    } else {
        println!("No telomere motifs found even with discovered database.");
    }
    
    println!("Analysis complete.");

    // Report most frequent telomere motifs
    report_most_frequent_motifs("initial_anno.txt", "anno.txt")?;

    Ok(())
}
