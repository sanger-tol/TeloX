use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::collections::HashMap;

const SEQ_NT4_TABLE: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = 1;
    table[b'C' as usize] = 2;
    table[b'G' as usize] = 3;
    table[b'T' as usize] = 4;
    table[b'a' as usize] = 1;
    table[b'c' as usize] = 2;
    table[b'g' as usize] = 3;
    table[b't' as usize] = 4;
    table
};

fn check_motif(motif: &str) -> bool {
    motif.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'))
}

fn encode_kmer(motif: &str) -> u64 {
    motif.bytes().fold(0, |acc, b| {
        (acc << 2) | SEQ_NT4_TABLE[b as usize] as u64
    })
}

fn encode_revcomp_kmer(motif: &str) -> u64 {
    motif
        .bytes()
        .rev()
        .fold(0, |acc, b| {
            let nt = SEQ_NT4_TABLE[b as usize];
            let rc = match nt {
                1 => 4, // A -> T
                2 => 3, // C -> G
                3 => 2, // G -> C
                4 => 1, // T -> A
                _ => 0,
            };
            (acc << 2) | rc as u64
        })
}

fn telo_finder_core(sequence: &str, motif: &str) -> u64 {
    let slen = sequence.len();
    if slen < motif.len() || motif.len() == 0 || !check_motif(motif) {
        return 0;
    }

    let mlen = motif.len();
    let mask: u64 = (1 << (2 * mlen)) - 1;
    let mut sum_telo = 0;

    let motif_val = encode_kmer(motif);
    let motif_rev = encode_revcomp_kmer(motif);

    // --- 5' end ---
    let mut max_score = 0;
    let mut max_i = None;
    let mut score = 0;
    let mut x = 0;

    for i in 0..slen {
        let nt = SEQ_NT4_TABLE[sequence.as_bytes()[i] as usize];
        if nt == 0 {
            x = 0;
            score = 0;
            continue;
        }

        x = ((x << 2) | nt as u64) & mask;

        if i + 1 >= mlen {
            if x == motif_val {
                score += 1;
            } else {
                score = 0;
            }
            if score > max_score {
                max_score = score;
                max_i = Some(i);
            }
        }
    }

    if let Some(i) = max_i {
        let end_pos = i + 1;
        sum_telo += (end_pos as u64) << 32;
    }

    // --- 3' end ---
    max_score = 0;
    max_i = None;
    score = 0;
    x = 0;

    for i in (0..slen).rev() {
        let nt = SEQ_NT4_TABLE[sequence.as_bytes()[i] as usize];
        if nt == 0 {
            x = 0;
            score = 0;
            continue;
        }

        // Reverse complement logic
        let rc = match nt {
            1 => 4,
            2 => 3,
            3 => 2,
            4 => 1,
            _ => 0,
        };

        x = ((x << 2) | rc as u64) & mask;

        if slen - i >= mlen {
            if x == motif_rev {
                score += 1;
            } else {
                score = 0;
            }
            if score > max_score {
                max_score = score;
                max_i = Some(i);
            }
        }
    }

    if let Some(i) = max_i {
        sum_telo += (slen - i) as u64;
    }

    sum_telo
}

pub fn telo_finder<P: AsRef<std::path::Path>>(
    fasta_file: P,
    motif: Option<&str>,
    mut output: Option<&mut dyn Write>,
) -> std::io::Result<()> {
    let motif = motif.unwrap_or("TTAGGG").to_string();

    let mut fasta = BufReader::new(File::open(fasta_file)?);
    let mut line = String::new();
    let mut name = String::new();
    let mut seq = String::new();

    while fasta.read_line(&mut line)? > 0 {
        let trimmed = line.trim();
        if trimmed.starts_with('>') {
            if !seq.is_empty() {
                let telo = telo_finder_core(&seq, &motif);
                if telo > 0 {
                    if let Some(w) = output.as_deref_mut() {
                        writeln!(w, "{}\t{}\t{}\t{}", name, telo >> 32, telo & 0xffffffff, motif)?;
                    } else {
                        println!("{}\t{}\t{}\t{}", name, telo >> 32, telo & 0xffffffff, motif);
                    }
                }
                seq.clear();
            }
            name = trimmed[1..].to_string();
        } else {
            seq.push_str(trimmed);
        }
        line.clear();
    }

    if !seq.is_empty() {
        let telo = telo_finder_core(&seq, &motif);
        if telo > 0 {
            if let Some(w) = output.as_deref_mut() {
                writeln!(w, "{}\t{}\t{}\t{}", name, telo >> 32, telo & 0xffffffff, motif)?;
            } else {
                println!("{}\t{}\t{}\t{}", name, telo >> 32, telo & 0xffffffff, motif);
            }
        }
    }

    Ok(())
}


fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <fasta_file> [motif]", args[0]);
        std::process::exit(1);
    }

    let motif = if args.len() >= 3 { Some(&args[2][..]) } else { None };
    telo_finder(&args[1], motif, None)
}
