/*
 * TeloX - Telomere Motif Extraction Tool
 * Genome Annotation Module
 *
 * Copyright (c) 2025 Yumi Sims, Chenxi Zhou, Will Eagles, Wellcome Sanger Institute
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

use needletail::parse_fastx_file;
use std::path::Path;
use anyhow::Result;
use std::io::Write;
use std::collections::{HashMap, HashSet};

/// Represents a telomere annotation in the genome
#[derive(Debug, Clone)]
pub struct TelomereAnnotation {
    pub chromosome: String,
    pub start: usize,
    pub end: usize,
    pub motif: String,
    pub strand: char,
    pub score: f64,
    pub region_type: String, // "telomere", "subtelomere", "telomere_cluster"
}

/// Configuration for genome annotation
#[derive(Debug, Clone)]
pub struct AnnotationConfig {
    pub min_motif_length: usize,
    pub max_gap_size: usize,
    pub min_cluster_size: usize,
    pub score_threshold: f64,
    pub include_reverse_complement: bool,
    pub merge_overlapping: bool,
    pub end_region_size: f64, // Fraction of sequence considered "end region" for scoring bonus
}

impl Default for AnnotationConfig {
    fn default() -> Self {
        Self {
            min_motif_length: 3,
            max_gap_size: 50,
            min_cluster_size: 2,
            score_threshold: 0.5,
            include_reverse_complement: true,
            merge_overlapping: true,
            end_region_size: 0.1, // 10% of sequence ends
        }
    }
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

/// Scan genome for telomere motif occurrences and return annotations
pub fn annotate_telomere_motifs(
    fasta_path: impl AsRef<Path>,
    motifs: &[String],
    config: &AnnotationConfig,
) -> Result<Vec<TelomereAnnotation>> {
    let mut annotations = Vec::new();
    let mut reader = parse_fastx_file(fasta_path)?;
    
    println!("Scanning genome for telomere motifs...");
    println!("Memory optimization: processing sequences individually");
    let mut sequence_count = 0;
    let mut total_sequence_length = 0;
    
    while let Some(record) = reader.next() {
        let seqrec = record?;
        let seq_name = std::str::from_utf8(seqrec.id())
            .unwrap_or("unknown")
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();
        
        // Process sequence in chunks to reduce memory usage
        let seq_bytes = seqrec.seq();
        let seq_len = seq_bytes.len();
        total_sequence_length += seq_len;
        
        sequence_count += 1;
        if sequence_count % 100 == 0 {
            println!("  Processed {} sequences ({:.1}M bp total)...", 
                     sequence_count, total_sequence_length as f64 / 1_000_000.0);
        }
        
        // Skip very small sequences to save memory
        if seq_len < 100 {
            continue;
        }
        
        // Convert to string only when needed, process in chunks for large sequences
        if seq_len > 10_000_000 {  // 10MB threshold
            println!("  Large sequence {}: {} bp - processing in chunks", seq_name, seq_len);
            // Process in 1MB chunks with overlap
            let chunk_size = 1_000_000;
            let overlap = 1000;
            
            for chunk_start in (0..seq_len).step_by(chunk_size - overlap) {
                let chunk_end = (chunk_start + chunk_size).min(seq_len);
                let chunk_bytes = &seq_bytes[chunk_start..chunk_end];
                
                if let Ok(chunk_sequence) = std::str::from_utf8(chunk_bytes) {
                    let chunk_sequence = chunk_sequence.to_uppercase();
                    
                    for motif in motifs {
                        let mut chunk_annotations = find_motif_occurrences(
                            &seq_name,
                            &chunk_sequence,
                            motif,
                            config,
                        );
                        
                        // Adjust coordinates to global sequence coordinates
                        for ann in &mut chunk_annotations {
                            ann.start += chunk_start;
                            ann.end += chunk_start;
                        }
                        
                        annotations.extend(chunk_annotations);
                    }
                }
                
                if chunk_end >= seq_len {
                    break;
                }
            }
        } else {
            // Process normally for smaller sequences
            if let Ok(sequence) = std::str::from_utf8(&seq_bytes) {
                let sequence = sequence.to_uppercase();
                
                for motif in motifs {
                    let motif_annotations = find_motif_occurrences(
                        &seq_name,
                        &sequence,
                        motif,
                        config,
                    );
                    annotations.extend(motif_annotations);
                }
            }
        }
    }
    
    println!("Processed {} sequences total", sequence_count);
    
    // Post-process annotations
    if config.merge_overlapping {
        println!("Merging overlapping annotations...");
        annotations = merge_overlapping_annotations(annotations, config.max_gap_size);
    }
    
    // Sort by chromosome and position
    annotations.sort_by(|a, b| {
        a.chromosome.cmp(&b.chromosome)
            .then(a.start.cmp(&b.start))
    });
    
    Ok(annotations)
}

/// Find all occurrences of a specific motif in a sequence
fn find_motif_occurrences(
    seq_name: &str,
    sequence: &str,
    motif: &str,
    config: &AnnotationConfig,
) -> Vec<TelomereAnnotation> {
    let mut annotations = Vec::new();
    let seq_bytes = sequence.as_bytes();
    let motif_bytes = motif.as_bytes();
    let motif_len = motif.len();
    
    if motif_len < config.min_motif_length || sequence.len() < motif_len {
        return annotations;
    }
    
    // Search for forward motif
    let mut i = 0;
    while i + motif_len <= seq_bytes.len() {
        if &seq_bytes[i..i + motif_len] == motif_bytes {
            // Calculate score based on local context
            let score = calculate_motif_score(seq_bytes, i, motif_len, config);
            
            if score >= config.score_threshold {
                annotations.push(TelomereAnnotation {
                    chromosome: seq_name.to_string(),
                    start: i,
                    end: i + motif_len,
                    motif: motif.to_string(),
                    strand: '+',
                    score,
                    region_type: "telomere".to_string(),
                });
            }
            i += motif_len; // Non-overlapping search
        } else {
            i += 1;
        }
    }
    
    // Search for reverse complement if enabled
    if config.include_reverse_complement {
        let rc_motif = reverse_complement(motif);
        let rc_bytes = rc_motif.as_bytes();
        
        let mut i = 0;
        while i + motif_len <= seq_bytes.len() {
            if &seq_bytes[i..i + motif_len] == rc_bytes {
                let score = calculate_motif_score(seq_bytes, i, motif_len, config);
                
                if score >= config.score_threshold {
                    annotations.push(TelomereAnnotation {
                        chromosome: seq_name.to_string(),
                        start: i,
                        end: i + motif_len,
                        motif: rc_motif.clone(),
                        strand: '-',
                        score,
                        region_type: "telomere".to_string(),
                    });
                }
                i += motif_len;
            } else {
                i += 1;
            }
        }
    }
    
    annotations
}

/// Calculate a quality score for a motif occurrence based on local context
fn calculate_motif_score(
    sequence: &[u8],
    pos: usize,
    motif_len: usize,
    config: &AnnotationConfig,
) -> f64 {
    let seq_len = sequence.len();
    let window_size = motif_len * 5; // Look at 5x motif length around the match
    
    let start = pos.saturating_sub(window_size);
    let end = (pos + motif_len + window_size).min(seq_len);
    
    if end <= start {
        return 0.5; // Default score
    }
    
    let window = &sequence[start..end];
    let motif = &sequence[pos..pos + motif_len];
    
    // Count motif occurrences in the window
    let mut motif_count = 0;
    let mut i = 0;
    while i + motif_len <= window.len() {
        if &window[i..i + motif_len] == motif {
            motif_count += 1;
            i += motif_len;
        } else {
            i += 1;
        }
    }
    
    // Score based on local density
    let max_possible_motifs = window.len() / motif_len;
    let density = if max_possible_motifs > 0 {
        motif_count as f64 / max_possible_motifs as f64
    } else {
        0.0
    };
    
    // Bonus for being near sequence ends (telomeres are typically at ends)
    let distance_from_start = pos as f64 / seq_len as f64;
    let distance_from_end = (seq_len - pos - motif_len) as f64 / seq_len as f64;
    let end_bonus = if distance_from_start < config.end_region_size || 
                       distance_from_end < config.end_region_size {
        0.3
    } else {
        0.0
    };
    
    // Base score from density + end bonus
    let base_score = density * 0.7 + 0.3; // Ensure minimum base score
    (base_score + end_bonus).min(1.0)
}

/// Merge overlapping or nearby annotations
fn merge_overlapping_annotations(
    mut annotations: Vec<TelomereAnnotation>,
    max_gap: usize,
) -> Vec<TelomereAnnotation> {
    if annotations.is_empty() {
        return annotations;
    }
    
    // Sort by chromosome and position
    annotations.sort_by(|a, b| {
        a.chromosome.cmp(&b.chromosome)
            .then(a.start.cmp(&b.start))
    });
    
    let mut merged = Vec::new();
    let mut current = annotations[0].clone();
    
    for annotation in annotations.into_iter().skip(1) {
        // Check if we should merge with current
        if annotation.chromosome == current.chromosome &&
           annotation.strand == current.strand &&
           annotation.start <= current.end + max_gap {
            
            // Merge annotations
            current.end = current.end.max(annotation.end);
            current.score = (current.score + annotation.score) / 2.0;
            
            // Update motif to show it's a cluster
            if current.motif != annotation.motif {
                current.motif = format!("{}|{}", current.motif, annotation.motif);
            }
            current.region_type = "telomere_cluster".to_string();
        } else {
            // Start new annotation
            merged.push(current);
            current = annotation;
        }
    }
    merged.push(current);
    
    merged
}

/// Write annotations to BED format file
pub fn write_bed_file(
    annotations: &[TelomereAnnotation],
    output_path: &str,
    include_header: bool,
) -> Result<()> {
    let mut file = std::fs::File::create(output_path)?;
    
    if include_header {
        writeln!(file, "track name=\"TeloX_Telomeres\" description=\"Telomere motifs identified by TeloX\" itemRgb=\"On\"")?;
    }
    
    for annotation in annotations {
        // Simple BED format: chr start end motif strand color
        let color = match annotation.strand {
            '+' => "255,0,0",    // Red for forward strand
            '-' => "0,0,255",    // Blue for reverse strand
            _ => "128,128,128",  // Gray for unknown
        };
        
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}",
            annotation.chromosome,
            annotation.start,
            annotation.end,
            annotation.motif,
            annotation.strand,
            color
        )?;
    }
    
    Ok(())
}

/// Generate comprehensive telomere annotation report
pub fn generate_annotation_report(
    annotations: &[TelomereAnnotation],
    output_path: &str,
) -> Result<()> {
    let mut file = std::fs::File::create(output_path)?;
    
    writeln!(file, "# TeloX Telomere Annotation Report")?;
    writeln!(file, "# Generated: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC"))?;
    writeln!(file, "")?;
    
    // Summary statistics
    let total_annotations = annotations.len();
    let chromosomes: HashSet<_> = annotations.iter()
        .map(|a| &a.chromosome)
        .collect();
    let unique_motifs: HashSet<_> = annotations.iter()
        .map(|a| &a.motif)
        .collect();
    
    writeln!(file, "## Summary Statistics")?;
    writeln!(file, "Total telomere annotations: {}", total_annotations)?;
    writeln!(file, "Chromosomes with telomeres: {}", chromosomes.len())?;
    writeln!(file, "Unique motifs found: {}", unique_motifs.len())?;
    writeln!(file, "")?;
    
    // Motif frequency table
    let mut motif_counts: HashMap<String, usize> = HashMap::new();
    for annotation in annotations {
        *motif_counts.entry(annotation.motif.clone()).or_insert(0) += 1;
    }
    
    let mut sorted_motifs: Vec<_> = motif_counts.iter().collect();
    sorted_motifs.sort_by(|a, b| b.1.cmp(a.1));
    
    writeln!(file, "## Motif Frequency")?;
    writeln!(file, "Motif\tCount\tPercentage")?;
    for (motif, count) in sorted_motifs {
        let percentage = (*count as f64 / total_annotations as f64) * 100.0;
        writeln!(file, "{}\t{}\t{:.1}%", motif, count, percentage)?;
    }
    writeln!(file, "")?;
    
    // Chromosome distribution
    let mut chr_counts: HashMap<String, usize> = HashMap::new();
    for annotation in annotations {
        *chr_counts.entry(annotation.chromosome.clone()).or_insert(0) += 1;
    }
    
    let mut sorted_chrs: Vec<_> = chr_counts.iter().collect();
    sorted_chrs.sort_by(|a, b| b.1.cmp(a.1));
    
    writeln!(file, "## Chromosome Distribution")?;
    writeln!(file, "Chromosome\tTelomere_Count")?;
    for (chr, count) in sorted_chrs {
        writeln!(file, "{}\t{}", chr, count)?;
    }
    writeln!(file, "")?;
    
    // Region type distribution
    let mut type_counts: HashMap<String, usize> = HashMap::new();
    for annotation in annotations {
        *type_counts.entry(annotation.region_type.clone()).or_insert(0) += 1;
    }
    
    writeln!(file, "## Region Type Distribution")?;
    writeln!(file, "Region_Type\tCount")?;
    for (region_type, count) in type_counts {
        writeln!(file, "{}\t{}", region_type, count)?;
    }
    
    Ok(())
}

/// Write detailed annotations to TSV format
pub fn write_annotation_tsv(annotations: &[TelomereAnnotation], output_path: &str) -> Result<()> {
    let mut file = std::fs::File::create(output_path)?;
    
    writeln!(file, "chromosome\tstart\tend\tmotif\tstrand\tscore\tregion_type\tlength")?;
    
    for annotation in annotations {
        let length = annotation.end - annotation.start;
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
            annotation.chromosome,
            annotation.start,
            annotation.end,
            annotation.motif,
            annotation.strand,
            annotation.score,
            annotation.region_type,
            length
        )?;
    }
    
    Ok(())
}

/// Main function to annotate genome with discovered telomere motifs
pub fn annotate_genome_with_motifs(
    fasta_path: &str,
    motifs: &[String],
    output_prefix: &str,
    config: Option<AnnotationConfig>,
) -> Result<Vec<TelomereAnnotation>> {
    let config = config.unwrap_or_default();
    
    println!("\n=== GENOME ANNOTATION WITH TELOMERE MOTIFS ===");
    println!("Annotating genome with {} telomere motifs...", motifs.len());
    println!("Motifs: {:?}", motifs);
    println!();
    println!("Configuration:");
    println!("  Min motif length: {}", config.min_motif_length);
    println!("  Max gap size: {}", config.max_gap_size);
    println!("  Min cluster size: {}", config.min_cluster_size);
    println!("  Score threshold: {:.2}", config.score_threshold);
    println!("  Include reverse complement: {}", config.include_reverse_complement);
    println!("  Merge overlapping: {}", config.merge_overlapping);
    println!("  End region size: {:.1}%", config.end_region_size * 100.0);
    println!();
    
    // Annotate genome
    let annotations = annotate_telomere_motifs(fasta_path, motifs, &config)?;
    
    println!("Found {} telomere annotations", annotations.len());
    
    if annotations.is_empty() {
        println!("No telomere annotations found with current parameters.");
        println!("Debugging information:");
        println!("  - Motifs searched: {:?}", motifs);
        println!("  - Score threshold: {:.2}", config.score_threshold);
        println!("  - Min motif length: {}", config.min_motif_length);
        println!("  - Include reverse complement: {}", config.include_reverse_complement);
        println!();
        println!("Consider:");
        println!("  - Lowering score threshold: --annotation-threshold 0.3");
        println!("  - Checking if motifs are correct");
        println!("  - Verifying FASTA file format");
        return Ok(Vec::new());
    }
    
    // Write BED file
    let bed_file = format!("{}_telomeres.bed", output_prefix);
    write_bed_file(&annotations, &bed_file, true)?;
    println!("BED file written to: {}", bed_file);
    
    // Write detailed report
    let report_file = format!("{}_annotation_report.txt", output_prefix);
    generate_annotation_report(&annotations, &report_file)?;
    println!("Annotation report written to: {}", report_file);
    
    // Write detailed TSV
    let tsv_file = format!("{}_annotations.tsv", output_prefix);
    write_annotation_tsv(&annotations, &tsv_file)?;
    println!("Detailed annotations written to: {}", tsv_file);
    
    // Summary statistics
    let chromosomes: HashSet<_> = annotations.iter().map(|a| &a.chromosome).collect();
    let unique_motifs: HashSet<_> = annotations.iter().map(|a| &a.motif).collect();
    
    println!();
    println!("=== ANNOTATION SUMMARY ===");
    println!("Total annotations: {}", annotations.len());
    println!("Chromosomes annotated: {}", chromosomes.len());
    println!("Unique motifs found: {}", unique_motifs.len());
    
    // Show top chromosomes
    let mut chr_counts: HashMap<String, usize> = HashMap::new();
    for annotation in &annotations {
        *chr_counts.entry(annotation.chromosome.clone()).or_insert(0) += 1;
    }
    let mut sorted_chrs: Vec<_> = chr_counts.iter().collect();
    sorted_chrs.sort_by(|a, b| b.1.cmp(a.1));
    
    println!("Top chromosomes by annotation count:");
    for (chr, count) in sorted_chrs.iter().take(5) {
        println!("  {}: {} annotations", chr, count);
    }
    
    Ok(annotations)
}

/// Create annotation configuration with custom parameters
pub fn create_annotation_config(
    score_threshold: Option<f64>,
    max_gap_size: Option<usize>,
    merge_overlapping: Option<bool>,
) -> AnnotationConfig {
    let mut config = AnnotationConfig::default();
    
    if let Some(threshold) = score_threshold {
        config.score_threshold = threshold;
    }
    if let Some(gap_size) = max_gap_size {
        config.max_gap_size = gap_size;
    }
    if let Some(merge) = merge_overlapping {
        config.merge_overlapping = merge;
    }
    
    config
}
