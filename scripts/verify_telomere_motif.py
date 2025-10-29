#!/usr/bin/env python3
"""
TeloX Telomere Motif Verification Script
Comprehensive validation plots to verify identified telomere motifs

Usage:
    python3 verify_telomere_motif.py telox_genome_telomeres.bed
    python3 verify_telomere_motif.py telox_genome_annotations.tsv --format tsv
"""

import sys
import argparse
import json
from collections import defaultdict, Counter
import numpy as np

def load_bed_file(bed_file):
    """Load telomere annotations from BED file"""
    annotations = []
    
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('track') or not line:
                continue
                
            cols = line.split('\t')
            if len(cols) >= 5:
                annotations.append({
                    'chr': cols[0],
                    'start': int(cols[1]),
                    'end': int(cols[2]),
                    'motif': cols[3],
                    'strand': cols[4],
                    'length': int(cols[2]) - int(cols[1])
                })
    
    return annotations

def load_motif_info():
    """Load primary motif information from TeloX JSON output"""
    try:
        with open('telomere_motif_final.json', 'r') as f:
            data = json.load(f)
            return data.get('most_likely_telomere_motif', {})
    except FileNotFoundError:
        print("Warning: telomere_motif_final.json not found")
        return {}

def generate_rotational_variants(motif):
    """Generate all rotational variants of a motif"""
    variants = []
    for i in range(len(motif)):
        variant = motif[i:] + motif[:i]
        variants.append(variant)
    
    # Add reverse complement variants
    rc_motif = reverse_complement(motif)
    for i in range(len(rc_motif)):
        variant = rc_motif[i:] + rc_motif[:i]
        variants.append(variant)
    
    return list(set(variants))  # Remove duplicates

def reverse_complement(seq):
    """Generate reverse complement"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def calculate_chromosome_lengths(annotations):
    """Estimate chromosome lengths from annotations"""
    chr_lengths = {}
    for ann in annotations:
        chr_name = ann['chr']
        if chr_name not in chr_lengths:
            chr_lengths[chr_name] = ann['end']
        else:
            chr_lengths[chr_name] = max(chr_lengths[chr_name], ann['end'])
    
    # Add buffer (telomeres might not be at absolute ends)
    for chr_name in chr_lengths:
        chr_lengths[chr_name] = int(chr_lengths[chr_name] * 1.1)
    
    return chr_lengths

def plot_chromosome_end_enrichment(annotations, output_prefix):
    """Plot telomere enrichment at chromosome ends"""
    try:
        import matplotlib.pyplot as plt
        
        print("📍 Generating chromosome end enrichment plot...")
        
        chr_lengths = calculate_chromosome_lengths(annotations)
        
        # Calculate distances from chromosome ends
        distances_from_ends = []
        
        for ann in annotations:
            chr_len = chr_lengths.get(ann['chr'], ann['end'] * 2)
            
            # Distance from start
            dist_from_start = ann['start']
            # Distance from end  
            dist_from_end = chr_len - ann['end']
            
            # Use minimum distance (closest end)
            min_distance = min(dist_from_start, dist_from_end)
            distances_from_ends.append(min_distance)
        
        # Create bins
        bins = [0, 1000, 5000, 10000, 50000, 100000, float('inf')]
        bin_labels = ['0-1kb', '1-5kb', '5-10kb', '10-50kb', '50-100kb', '>100kb']
        
        # Count annotations in each bin
        bin_counts = [0] * len(bin_labels)
        for dist in distances_from_ends:
            for i, bin_edge in enumerate(bins[1:]):
                if dist <= bin_edge:
                    bin_counts[i] += 1
                    break
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.bar(bin_labels, bin_counts, color='darkgreen', alpha=0.7, edgecolor='black')
        
        # Add value labels
        for bar, count in zip(bars, bin_counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(bin_counts)*0.01,
                   str(count), ha='center', va='bottom', fontweight='bold')
        
        ax.set_xlabel('Distance from Chromosome End')
        ax.set_ylabel('Number of Telomere Annotations')
        ax.set_title('Telomere Motif Enrichment at Chromosome Ends')
        ax.grid(True, alpha=0.3)
        
        # Add interpretation text
        end_enrichment = (bin_counts[0] + bin_counts[1]) / sum(bin_counts) * 100
        ax.text(0.02, 0.98, f'End enrichment (0-5kb): {end_enrichment:.1f}%', 
               transform=ax.transAxes, bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7),
               verticalalignment='top', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_end_enrichment.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   ✓ End enrichment plot saved: {output_prefix}_end_enrichment.png")
        
        return end_enrichment
        
    except ImportError:
        print("❌ Matplotlib not available for end enrichment plot")
        return 0

def plot_strand_bias_validation(annotations, motif_info, output_prefix):
    """Plot strand bias validation"""
    try:
        import matplotlib.pyplot as plt
        
        print("⚡ Generating strand bias validation plot...")
        
        # Group by motif and calculate bias
        motif_bias = {}
        for ann in annotations:
            motif = ann['motif']
            if motif not in motif_bias:
                motif_bias[motif] = {'forward': 0, 'reverse': 0}
            
            if ann['strand'] == '+':
                motif_bias[motif]['forward'] += 1
            else:
                motif_bias[motif]['reverse'] += 1
        
        # Create scatter plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        motifs = list(motif_bias.keys())
        forward_counts = [motif_bias[m]['forward'] for m in motifs]
        reverse_counts = [motif_bias[m]['reverse'] for m in motifs]
        
        # Color by bias ratio
        colors = []
        bias_ratios = []
        for f, r in zip(forward_counts, reverse_counts):
            if r > 0:
                ratio = f / r
            else:
                ratio = float('inf')
            bias_ratios.append(ratio)
            
            if ratio > 2:
                colors.append('red')  # Forward biased
            elif ratio < 0.5:
                colors.append('blue')  # Reverse biased
            else:
                colors.append('gray')  # Balanced
        
        scatter = ax.scatter(forward_counts, reverse_counts, c=colors, s=100, alpha=0.7, edgecolors='black')
        
        # Add diagonal lines for bias ratios
        max_count = max(max(forward_counts), max(reverse_counts))
        x_line = np.linspace(0, max_count, 100)
        
        ax.plot(x_line, x_line, 'k--', alpha=0.5, label='1:1 (balanced)')
        ax.plot(x_line, x_line/2, 'r--', alpha=0.5, label='2:1 (forward bias)')
        ax.plot(x_line, x_line/5, 'r:', alpha=0.5, label='5:1 (strong forward)')
        ax.plot(x_line*2, x_line, 'b--', alpha=0.5, label='1:2 (reverse bias)')
        ax.plot(x_line*5, x_line, 'b:', alpha=0.5, label='1:5 (strong reverse)')
        
        # Add motif labels
        for i, motif in enumerate(motifs):
            ax.annotate(motif, (forward_counts[i], reverse_counts[i]), 
                       xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        ax.set_xlabel('Forward Strand Count')
        ax.set_ylabel('Reverse Strand Count')
        ax.set_title('Strand Bias Validation for Telomere Motifs')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_strand_bias.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   ✓ Strand bias plot saved: {output_prefix}_strand_bias.png")
        
    except ImportError:
        print("❌ Matplotlib not available for strand bias plot")

def plot_motif_length_distribution(annotations, output_prefix):
    """Plot distribution of telomeric tract lengths"""
    try:
        import matplotlib.pyplot as plt
        
        print("📏 Generating motif length distribution plot...")
        
        # Group consecutive annotations to find tract lengths
        chr_data = defaultdict(list)
        for ann in annotations:
            chr_data[ann['chr']].append(ann)
        
        tract_lengths = []
        
        for chr_name, chr_annotations in chr_data.items():
            # Sort by position
            chr_annotations.sort(key=lambda x: x['start'])
            
            # Group consecutive annotations
            current_tract_start = None
            current_tract_end = None
            
            for ann in chr_annotations:
                if current_tract_start is None:
                    current_tract_start = ann['start']
                    current_tract_end = ann['end']
                elif ann['start'] - current_tract_end <= 1000:  # Within 1kb = same tract
                    current_tract_end = ann['end']
                else:
                    # End of current tract
                    tract_lengths.append(current_tract_end - current_tract_start)
                    current_tract_start = ann['start']
                    current_tract_end = ann['end']
            
            # Add final tract
            if current_tract_start is not None:
                tract_lengths.append(current_tract_end - current_tract_start)
        
        # Create histogram
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Linear scale
        ax1.hist(tract_lengths, bins=30, color='purple', alpha=0.7, edgecolor='black')
        ax1.set_xlabel('Telomeric Tract Length (bp)')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Telomeric Tract Length Distribution')
        ax1.grid(True, alpha=0.3)
        
        # Add statistics
        mean_length = np.mean(tract_lengths)
        median_length = np.median(tract_lengths)
        ax1.axvline(mean_length, color='red', linestyle='--', label=f'Mean: {mean_length:.0f} bp')
        ax1.axvline(median_length, color='orange', linestyle='--', label=f'Median: {median_length:.0f} bp')
        ax1.legend()
        
        # Log scale for wide range
        ax2.hist(tract_lengths, bins=30, color='purple', alpha=0.7, edgecolor='black')
        ax2.set_xlabel('Telomeric Tract Length (bp)')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Telomeric Tract Length Distribution (Log Scale)')
        ax2.set_yscale('log')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_tract_lengths.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   ✓ Tract length plot saved: {output_prefix}_tract_lengths.png")
        
        return tract_lengths
        
    except ImportError:
        print("❌ Matplotlib not available for length distribution plot")
        return []

def plot_rotational_variants(primary_motif, annotations, output_prefix):
    """Plot frequency of rotational variants"""
    try:
        import matplotlib.pyplot as plt
        
        print("🔄 Generating rotational variants plot...")
        
        # Generate all expected rotational variants
        expected_variants = generate_rotational_variants(primary_motif)
        
        # Count actual variants found
        motif_counts = Counter(ann['motif'] for ann in annotations)
        
        # Separate found vs expected
        found_variants = []
        found_counts = []
        missing_variants = []
        
        for variant in expected_variants:
            if variant in motif_counts:
                found_variants.append(variant)
                found_counts.append(motif_counts[variant])
            else:
                missing_variants.append(variant)
        
        # Create plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Found variants
        bars1 = ax1.bar(range(len(found_variants)), found_counts, 
                       color='green', alpha=0.7, edgecolor='black')
        ax1.set_xlabel('Rotational Variant')
        ax1.set_ylabel('Count')
        ax1.set_title(f'Found Rotational Variants of {primary_motif}')
        ax1.set_xticks(range(len(found_variants)))
        ax1.set_xticklabels(found_variants, rotation=45, ha='right')
        
        # Add value labels
        for bar, count in zip(bars1, found_counts):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(found_counts)*0.01,
                    str(count), ha='center', va='bottom', fontweight='bold')
        
        # Coverage pie chart
        found_ratio = len(found_variants) / len(expected_variants)
        missing_ratio = len(missing_variants) / len(expected_variants)
        
        ax2.pie([found_ratio, missing_ratio], 
               labels=[f'Found ({len(found_variants)})', f'Missing ({len(missing_variants)})'],
               colors=['green', 'lightcoral'],
               autopct='%1.1f%%',
               startangle=90)
        ax2.set_title(f'Rotational Variant Coverage\\n({len(found_variants)}/{len(expected_variants)} variants found)')
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_rotational_variants.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   ✓ Rotational variants plot saved: {output_prefix}_rotational_variants.png")
        
        return found_ratio
        
    except ImportError:
        print("❌ Matplotlib not available for rotational variants plot")
        return 0

def plot_motif_density_heatmap(annotations, output_prefix):
    """Create heatmap of motif density along chromosomes"""
    try:
        import matplotlib.pyplot as plt
        
        print("🎯 Generating motif density heatmap...")
        
        # Group by chromosome
        chr_data = defaultdict(list)
        for ann in annotations:
            chr_data[ann['chr']].append(ann)
        
        # Calculate density for each chromosome
        chromosomes = sorted(chr_data.keys())[:20]  # Limit to top 20 chromosomes
        
        if not chromosomes:
            print("   No chromosomes to plot")
            return
        
        # Create density matrix
        bin_size = 100000  # 100kb bins
        max_bins = 100     # Maximum bins per chromosome
        
        density_matrix = []
        chr_labels = []
        
        for chr_name in chromosomes:
            chr_annotations = chr_data[chr_name]
            chr_length = max(ann['end'] for ann in chr_annotations)
            
            # Create bins
            num_bins = min(max_bins, max(10, chr_length // bin_size))
            bins = np.linspace(0, chr_length, num_bins + 1)
            
            # Count annotations in each bin
            bin_counts = np.zeros(num_bins)
            for ann in chr_annotations:
                bin_idx = int(ann['start'] / chr_length * num_bins)
                bin_idx = min(bin_idx, num_bins - 1)
                bin_counts[bin_idx] += 1
            
            density_matrix.append(bin_counts)
            chr_labels.append(chr_name)
        
        # Pad shorter chromosomes
        max_bins_actual = max(len(row) for row in density_matrix)
        for i, row in enumerate(density_matrix):
            if len(row) < max_bins_actual:
                density_matrix[i] = np.pad(row, (0, max_bins_actual - len(row)), 'constant')
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(15, max(8, len(chromosomes) * 0.5)))
        
        im = ax.imshow(density_matrix, cmap='Reds', aspect='auto', interpolation='nearest')
        
        # Labels
        ax.set_xlabel('Relative Position Along Chromosome')
        ax.set_ylabel('Chromosome')
        ax.set_title('Telomere Motif Density Heatmap')
        ax.set_yticks(range(len(chr_labels)))
        ax.set_yticklabels(chr_labels)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Motif Count per Bin')
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_density_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   ✓ Density heatmap saved: {output_prefix}_density_heatmap.png")
        
    except ImportError:
        print("❌ Matplotlib not available for density heatmap")

def generate_verification_report(annotations, motif_info, output_prefix):
    """Generate comprehensive verification report"""
    print("\n📋 GENERATING VERIFICATION REPORT...")
    
    report_file = f"{output_prefix}_verification_report.txt"
    with open(report_file, 'w') as f:
        f.write("# TeloX Telomere Motif Verification Report\\n")
        f.write(f"# Generated for motif verification\\n\\n")
        
        # Basic statistics
        f.write("## Basic Statistics\\n")
        f.write(f"Total annotations: {len(annotations)}\\n")
        f.write(f"Unique chromosomes: {len(set(ann['chr'] for ann in annotations))}\\n")
        f.write(f"Unique motifs: {len(set(ann['motif'] for ann in annotations))}\\n")
        
        # Primary motif info
        if motif_info:
            f.write(f"\\n## Primary Motif Information\\n")
            f.write(f"Canonical sequence: {motif_info.get('canonical_sequence', 'N/A')}\\n")
            f.write(f"Rotational variants: {motif_info.get('rotational_variants', [])}\\n")
            
            stats = motif_info.get('statistics', {})
            f.write(f"Forward count: {stats.get('forward_count', 'N/A')}\\n")
            f.write(f"Reverse count: {stats.get('reverse_count', 'N/A')}\\n")
            f.write(f"Total frequency: {stats.get('total_frequency', 'N/A')}\\n")
            f.write(f"Longest stretch: {stats.get('longest_continuous_stretch', 'N/A')}\\n")
        
        # Strand bias analysis
        forward_count = sum(1 for ann in annotations if ann['strand'] == '+')
        reverse_count = sum(1 for ann in annotations if ann['strand'] == '-')
        bias_ratio = forward_count / reverse_count if reverse_count > 0 else float('inf')
        
        f.write(f"\\n## Strand Bias Analysis\\n")
        f.write(f"Forward annotations: {forward_count}\\n")
        f.write(f"Reverse annotations: {reverse_count}\\n")
        f.write(f"Bias ratio: {bias_ratio:.2f}\\n")
        
        bias_strength = "Strong" if bias_ratio > 2 or bias_ratio < 0.5 else "Weak"
        f.write(f"Bias strength: {bias_strength}\\n")
        
        # Length statistics
        lengths = [ann['length'] for ann in annotations]
        f.write(f"\\n## Length Statistics\\n")
        f.write(f"Mean motif length: {np.mean(lengths):.1f} bp\\n")
        f.write(f"Median motif length: {np.median(lengths):.1f} bp\\n")
        f.write(f"Length range: {min(lengths)} - {max(lengths)} bp\\n")
        
    print(f"   ✓ Verification report saved: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='Verify TeloX telomere motif identification')
    parser.add_argument('input_file', help='Input BED file from TeloX annotation')
    parser.add_argument('--output-prefix', default='verification',
                       help='Output file prefix (default: verification)')
    
    args = parser.parse_args()
    
    try:
        # Load data
        annotations = load_bed_file(args.input_file)
        motif_info = load_motif_info()
        
        print(f"📁 Loaded {len(annotations)} annotations for verification")
        
        if not annotations:
            print("❌ No annotations found. Run TeloX with --annotate first.")
            return
        
        primary_motif = motif_info.get('canonical_sequence', 
                                     Counter(ann['motif'] for ann in annotations).most_common(1)[0][0])
        
        print(f"🧬 Verifying primary motif: {primary_motif}")
        print(f"📊 Generating verification plots...")
        
        # Generate verification plots
        end_enrichment = plot_chromosome_end_enrichment(annotations, args.output_prefix)
        plot_strand_bias_validation(annotations, motif_info, args.output_prefix)
        tract_lengths = plot_motif_length_distribution(annotations, args.output_prefix)
        variant_coverage = plot_rotational_variants(primary_motif, annotations, args.output_prefix)
        plot_motif_density_heatmap(annotations, args.output_prefix)
        
        # Generate verification report
        generate_verification_report(annotations, motif_info, args.output_prefix)
        
        # Summary assessment
        print(f"\\n✅ VERIFICATION SUMMARY:")
        print(f"   Primary motif: {primary_motif}")
        print(f"   Total annotations: {len(annotations)}")
        print(f"   End enrichment: {end_enrichment:.1f}% (good if >50%)")
        print(f"   Rotational coverage: {variant_coverage:.1f} (good if >0.5)")
        
        # Overall assessment
        if end_enrichment > 50 and variant_coverage > 0.5:
            print(f"   🎉 ASSESSMENT: Strong telomere motif candidate!")
        elif end_enrichment > 30 or variant_coverage > 0.3:
            print(f"   ⚠️  ASSESSMENT: Moderate telomere motif candidate")
        else:
            print(f"   ❌ ASSESSMENT: Weak telomere motif candidate")
        
        print(f"\\n📁 Verification files generated:")
        print(f"   - {args.output_prefix}_end_enrichment.png")
        print(f"   - {args.output_prefix}_strand_bias.png") 
        print(f"   - {args.output_prefix}_tract_lengths.png")
        print(f"   - {args.output_prefix}_rotational_variants.png")
        print(f"   - {args.output_prefix}_density_heatmap.png")
        print(f"   - {args.output_prefix}_verification_report.txt")
        
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    main()
