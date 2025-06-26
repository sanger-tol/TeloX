use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::env;
use anyhow::{Context, Result};
mod kmers;
use kmers::{count_kmers_in_fasta, print_kmer_table, analyze_strand_bias, print_strand_bias_table, save_strand_bias_table, get_strand_bias_summary, longest_continuous_stretch_for_kmers, filter_strand_bias_tsvs, consolidate_rotational_kmers};

// Constants
const TELO_PENALTY: i64 = 1;
const TELO_MAX_DROP: i64 = 2000;
const TELO_MIN_SCORE: i64 = 300;
/*static TELO_MOTIF_DB: &[&str] = &[
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
];*/

static TELO_MOTIF_DB: &[&str] = &[
"CCTAA",
"AATTC",
"TTAGG",
"GCCTAA",
"TTAGGC",
"CCACAA",
"TTGTGG",
"CCCTAA",
"TTAGGG",
"TGCAA",
"TTGCA",
"CCCCAAAA",
"TTTTGGGG",
"CCCTAAAA",
"TTTTAGGG",
"CCCTA",
"TAGGG",
"CAATCGTCC",
"GGACGATTG",
"ACACCAGT",
"ACTGGTGT",
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
) -> io::Result<Vec<(bool, bool)>> {
    let sdict = SequenceDict::from_fasta(fasta_path)?;
    let nseq = sdict.sequences.len();
    let mut telo_ends = vec![(false, false); nseq];

    let motifs = if let Some(motif) = custom_motif {
        if !check_motif(motif) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Invalid motif characters",
            ));
        }
        vec![motif]
    } else {
        TELO_MOTIF_DB.to_vec()
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

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: {} <fasta_file>", args[0]);
        eprintln!("Example: {} genome.fasta", args[0]);
        std::process::exit(1);
    }

    let fasta_path = args[1].clone();

    for k in 5..=12 {
        // Use Rust k-mer counting for last 5000bp of each scaffold
        let counts = kmers::count_kmers_in_fasta(&fasta_path, k)
            .with_context(|| format!("Failed to count {}-mers in {}", k, &fasta_path))?;

        // Calculate longest stretch for each k-mer using the original FASTA
        let sdict = SequenceDict::from_fasta(&fasta_path)?;
        let kmer_list: Vec<String> = counts.keys().cloned().collect();
        let mut longest_stretch_map: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
        for seq in &sdict.sequences {
            let seq_slice = if seq.seq.len() > 5000 {
                &seq.seq[seq.seq.len() - 5000..]
            } else {
                &seq.seq
            };
            let stretch_map = kmers::longest_continuous_stretch_for_kmers(seq_slice, &kmer_list);
            for (kmer, stretch) in stretch_map {
                let entry = longest_stretch_map.entry(kmer).or_insert(0);
                if stretch > *entry {
                    *entry = stretch;
                }
            }
        }
        // Strand bias analysis
        let mut bias_analyses = kmers::analyze_strand_bias(&counts, Some(&longest_stretch_map));
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
    }

    Ok(())
}
