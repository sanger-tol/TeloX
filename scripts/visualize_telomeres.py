#!/usr/bin/env python3
"""
TeloX Telomere Visualization Script
Standalone visualization tool for TeloX output files

Usage:
    python3 visualize_telomeres.py telox_genome_telomeres.bed
    python3 visualize_telomeres.py telox_genome_annotations.tsv --format tsv
    python3 visualize_telomeres.py telox_genome_telomeres.bed --interactive
    python3 visualize_telomeres.py telox_genome_telomeres.bed --text-only
"""

import sys
import argparse
from collections import defaultdict, Counter

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
                    'color': cols[5] if len(cols) > 5 else '255,0,0',
                    'score': 0.8  # Default score for BED files
                })
    
    return annotations

def load_tsv_file(tsv_file):
    """Load telomere annotations from TSV file"""
    annotations = []
    
    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            cols = line.split('\t')
            if len(cols) >= len(header):
                annotations.append({
                    'chr': cols[0],
                    'start': int(cols[1]),
                    'end': int(cols[2]),
                    'motif': cols[3],
                    'strand': cols[4],
                    'score': float(cols[5]) if len(cols) > 5 else 0.5,
                    'color': '255,0,0' if cols[4] == '+' else '0,0,255'
                })
    
    return annotations

def create_text_summary(annotations):
    """Create comprehensive text-based summary"""
    print("\n" + "="*60)
    print("           TELOMERE VISUALIZATION SUMMARY")
    print("="*60)
    
    # Basic stats
    print(f"\n📊 OVERVIEW:")
    print(f"   Total annotations: {len(annotations)}")
    
    # Group by chromosome
    chr_data = defaultdict(list)
    for ann in annotations:
        chr_data[ann['chr']].append(ann)
    
    print(f"   Chromosomes with telomeres: {len(chr_data)}")
    print(f"   Average per chromosome: {len(annotations)/len(chr_data):.1f}")
    
    # Chromosome summary
    print(f"\n🧬 CHROMOSOME DISTRIBUTION:")
    print(f"{'Chromosome':<15} {'Count':<8} {'Forward':<8} {'Reverse':<8} {'Position Range':<20}")
    print("-" * 70)
    
    for chr_name in sorted(chr_data.keys()):
        chr_annotations = chr_data[chr_name]
        forward_count = sum(1 for ann in chr_annotations if ann['strand'] == '+')
        reverse_count = sum(1 for ann in chr_annotations if ann['strand'] == '-')
        
        positions = [ann['start'] for ann in chr_annotations]
        if positions:
            pos_range = f"{min(positions):,} - {max(positions):,}"
        else:
            pos_range = "N/A"
        
        print(f"{chr_name:<15} {len(chr_annotations):<8} {forward_count:<8} {reverse_count:<8} {pos_range:<20}")
    
    # Motif frequency
    motif_counts = Counter(ann['motif'] for ann in annotations)
    print(f"\n🎯 MOTIF FREQUENCY:")
    print(f"{'Motif':<12} {'Count':<8} {'Percentage':<12} {'Bar':<20}")
    print("-" * 55)
    
    max_count = max(motif_counts.values()) if motif_counts else 1
    for motif, count in motif_counts.most_common(10):
        percentage = (count / len(annotations)) * 100
        bar_length = int((count / max_count) * 20)
        bar = "█" * bar_length + "░" * (20 - bar_length)
        print(f"{motif:<12} {count:<8} {percentage:<11.1f}% {bar}")
    
    # Strand distribution
    forward_count = sum(1 for ann in annotations if ann['strand'] == '+')
    reverse_count = sum(1 for ann in annotations if ann['strand'] == '-')
    
    print(f"\n⚡ STRAND DISTRIBUTION:")
    print(f"   Forward (+): {forward_count:>6} ({forward_count/len(annotations)*100:>5.1f}%) {'█' * int(forward_count/len(annotations)*30)}")
    print(f"   Reverse (-): {reverse_count:>6} ({reverse_count/len(annotations)*100:>5.1f}%) {'█' * int(reverse_count/len(annotations)*30)}")

def create_ascii_ideogram(annotations):
    """Create ASCII ideogram of telomere distribution"""
    print(f"\n🗺️  ASCII CHROMOSOME IDEOGRAM")
    print("="*60)
    
    # Group by chromosome
    chr_data = defaultdict(list)
    for ann in annotations:
        chr_data[ann['chr']].append(ann)
    
    for chr_name in sorted(chr_data.keys()):
        chr_annotations = chr_data[chr_name]
        if not chr_annotations:
            continue
            
        chr_length = max(ann['end'] for ann in chr_annotations)
        
        print(f"\n{chr_name}: ({len(chr_annotations)} telomeres, {chr_length:,} bp)")
        
        # Create ASCII representation
        bar_length = 50
        bar = [' '] * bar_length
        
        for ann in chr_annotations:
            start_pos = int(ann['start'] / chr_length * bar_length)
            end_pos = int(ann['end'] / chr_length * bar_length)
            
            symbol = '>' if ann['strand'] == '+' else '<'
            
            for i in range(start_pos, min(end_pos + 1, bar_length)):
                if bar[i] == ' ':
                    bar[i] = symbol
                else:
                    bar[i] = '|'  # Overlapping regions
        
        bar_str = ''.join(bar)
        print(f"   |{bar_str}|")
        print(f"   0{chr_length:>49,}")
    
    print(f"\n   Legend: > = forward (+), < = reverse (-), | = overlapping")

def create_matplotlib_plots(annotations, output_prefix="telomere"):
    """Create matplotlib plots with memory optimization"""
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        
        print(f"\n📈 GENERATING MATPLOTLIB PLOTS...")
        print(f"   Processing {len(annotations)} annotations...")
        
        # Memory optimization: limit number of annotations if too many
        if len(annotations) > 10000:
            print(f"   Large dataset detected ({len(annotations)} annotations)")
            print(f"   Sampling top 5000 annotations by score for visualization...")
            
            # Sort by score and take top annotations
            sorted_annotations = sorted(annotations, key=lambda x: x.get('score', 0.5), reverse=True)
            annotations = sorted_annotations[:5000]
        
        # Group by chromosome
        chr_data = defaultdict(list)
        for ann in annotations:
            chr_data[ann['chr']].append(ann)
        
        chromosomes = sorted(chr_data.keys())
        
        # Limit number of chromosomes to prevent memory issues
        if len(chromosomes) > 50:
            print(f"   Too many chromosomes ({len(chromosomes)}), showing top 20 by annotation count")
            chr_counts = {chr: len(anns) for chr, anns in chr_data.items()}
            top_chromosomes = sorted(chr_counts.keys(), key=lambda x: chr_counts[x], reverse=True)[:20]
            chromosomes = sorted(top_chromosomes)
            chr_data = {chr: chr_data[chr] for chr in chromosomes}
        
        # 1. Create ideogram with memory-efficient settings
        plt.rcParams['figure.max_open_warning'] = 0  # Disable warning
        fig, axes = plt.subplots(len(chromosomes), 1, 
                                figsize=(12, min(20, max(6, len(chromosomes) * 0.6))),
                                dpi=150)  # Lower DPI to save memory
        
        if len(chromosomes) == 1:
            axes = [axes]
        
        for i, chr_name in enumerate(chromosomes):
            ax = axes[i]
            chr_annotations = chr_data[chr_name]
            
            # Calculate chromosome length
            chr_length = max(ann['end'] for ann in chr_annotations) if chr_annotations else 1000000
            
            # Draw chromosome backbone
            ax.barh(0, chr_length, height=0.3, color='lightgray', alpha=0.3, edgecolor='gray')
            
            # Add telomere regions
            for ann in chr_annotations:
                color = 'red' if ann['strand'] == '+' else 'blue'
                width = ann['end'] - ann['start']
                ax.barh(0, width, left=ann['start'], height=0.3, 
                       color=color, alpha=0.8, edgecolor='black', linewidth=0.5)
            
            # Formatting
            ax.set_ylabel(chr_name, rotation=0, ha='right', va='center', fontsize=10)
            ax.set_xlim(0, chr_length * 1.05)
            ax.set_ylim(-0.2, 0.2)
            ax.set_yticks([])
            
            # Clean up spines
            for spine in ax.spines.values():
                spine.set_visible(False)
            
            # Add telomere count
            ax.text(chr_length * 1.02, 0, f'n={len(chr_annotations)}', 
                   va='center', fontsize=8, color='gray')
        
        # Final formatting
        axes[-1].set_xlabel("Genomic Position (bp)", fontsize=12)
        plt.suptitle("Telomere Distribution Along Chromosomes", fontsize=16, fontweight='bold')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', alpha=0.8, label='Forward strand (+)'),
            Patch(facecolor='blue', alpha=0.8, label='Reverse strand (-)')
        ]
        fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.95))
        
        plt.tight_layout()
        ideogram_file = f"{output_prefix}_ideogram.png"
        plt.savefig(ideogram_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"   ✓ Ideogram saved: {ideogram_file}")
        
        # 2. Create summary plots
        create_summary_plots(annotations, output_prefix)
        
    except ImportError as e:
        print(f"❌ Matplotlib not available: {e}")
        print("   Install with: pip install matplotlib numpy")
        print("   Or: conda install matplotlib numpy")

def create_summary_plots(annotations, output_prefix):
    """Create 4-panel summary plots with memory optimization"""
    import matplotlib.pyplot as plt
    import numpy as np
    
    print(f"   Creating summary plots...")
    
    # Memory optimization: sample data if too large
    plot_annotations = annotations
    if len(annotations) > 5000:
        print(f"   Sampling {min(5000, len(annotations))} annotations for summary plots...")
        plot_annotations = annotations[:5000]
    
    # Extract data
    chromosomes = [ann['chr'] for ann in plot_annotations]
    motifs = [ann['motif'] for ann in plot_annotations]
    strands = [ann['strand'] for ann in plot_annotations]
    scores = [ann.get('score', 0.5) for ann in plot_annotations]
    
    # Create subplots with smaller figure size
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8), dpi=150)
    
    # 1. Chromosome distribution (limit to top 20)
    chr_counts = Counter(chromosomes)
    top_chrs = chr_counts.most_common(20)  # Limit to top 20
    chrs = [item[0] for item in top_chrs]
    counts = [item[1] for item in top_chrs]
    
    bars1 = ax1.bar(range(len(chrs)), counts, color='skyblue', edgecolor='navy', alpha=0.7)
    ax1.set_xlabel('Chromosome')
    ax1.set_ylabel('Telomere Count')
    ax1.set_title(f'Top {len(chrs)} Chromosomes by Telomere Count')
    ax1.set_xticks(range(len(chrs)))
    ax1.set_xticklabels(chrs, rotation=45, ha='right', fontsize=8)
    
    # Add value labels on bars
    for bar, count in zip(bars1, counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01,
                str(count), ha='center', va='bottom', fontsize=8)
    
    # 2. Motif frequency
    motif_counts = Counter(motifs)
    motif_names = list(motif_counts.keys())[:10]  # Top 10
    motif_vals = [motif_counts[m] for m in motif_names]
    
    bars2 = ax2.barh(range(len(motif_names)), motif_vals, color='lightcoral', edgecolor='darkred', alpha=0.7)
    ax2.set_xlabel('Count')
    ax2.set_ylabel('Motif')
    ax2.set_title('Top Telomere Motifs')
    ax2.set_yticks(range(len(motif_names)))
    ax2.set_yticklabels(motif_names)
    
    # Add value labels
    for bar, val in zip(bars2, motif_vals):
        ax2.text(bar.get_width() + max(motif_vals)*0.01, bar.get_y() + bar.get_height()/2,
                str(val), ha='left', va='center', fontsize=8)
    
    # 3. Strand bias
    strand_counts = Counter(strands)
    strand_labels = list(strand_counts.keys())
    strand_values = list(strand_counts.values())
    colors = ['lightcoral' if s == '+' else 'lightblue' for s in strand_labels]
    
    wedges, texts, autotexts = ax3.pie(strand_values, 
                                      labels=[f"{s} ({v})" for s, v in zip(strand_labels, strand_values)], 
                                      colors=colors, autopct='%1.1f%%')
    ax3.set_title('Strand Distribution')
    
    # 4. Score distribution (if available)
    if any(s > 0 for s in scores):
        try:
            n, bins, patches = ax4.hist(scores, bins=20, color='gold', edgecolor='orange')
            ax4.set_xlabel('Quality Score')
            ax4.set_ylabel('Frequency')
            ax4.set_title('Annotation Quality Scores')
            
            mean_score = np.mean(scores)
            ax4.axvline(mean_score, color='red', linestyle='--', linewidth=2,
                       label=f'Mean: {mean_score:.3f}')
            ax4.legend()
            
            # Color bars by score (simplified)
            for patch, bin_center in zip(patches, (bins[:-1] + bins[1:]) / 2):
                if bin_center > mean_score:
                    patch.set_facecolor('darkgreen')
        except Exception as e:
            print(f"   Warning: Could not create score histogram: {e}")
            ax4.text(0.5, 0.5, f'Score plot error\\n{str(e)[:50]}...', 
                    ha='center', va='center', transform=ax4.transAxes, fontsize=10)
            ax4.set_title('Quality Scores (Error)')
    else:
        ax4.text(0.5, 0.5, 'No score data available\\n(BED file format)', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12)
        ax4.set_title('Quality Scores (N/A)')
    
    plt.tight_layout()
    summary_file = f"{output_prefix}_summary.png"
    plt.savefig(summary_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"   ✓ Summary plots saved: {summary_file}")

def create_interactive_html(annotations, output_file="telomere_interactive.html"):
    """Create interactive HTML visualization using Plotly"""
    try:
        import plotly.graph_objects as go
        from plotly.offline import plot
        import plotly.express as px
        
        print(f"\n🌐 GENERATING INTERACTIVE HTML...")
        
        # Group by chromosome
        chr_data = defaultdict(list)
        for ann in annotations:
            chr_data[ann['chr']].append(ann)
        
        fig = go.Figure()
        
        # Add traces for each chromosome
        for i, (chr_name, chr_annotations) in enumerate(sorted(chr_data.items())):
            # Forward strand
            forward_anns = [ann for ann in chr_annotations if ann['strand'] == '+']
            if forward_anns:
                fig.add_trace(go.Scatter(
                    x=[ann['start'] for ann in forward_anns],
                    y=[chr_name] * len(forward_anns),
                    mode='markers',
                    marker=dict(
                        color='red',
                        size=10,
                        symbol='triangle-right',
                        line=dict(width=1, color='darkred')
                    ),
                    name='Forward (+)' if i == 0 else '',
                    text=[f"{ann['motif']} ({ann['start']:,}-{ann['end']:,})" for ann in forward_anns],
                    hovertemplate='<b>%{y}</b><br>Position: %{x:,}<br>%{text}<br>Strand: +<extra></extra>',
                    showlegend=(i == 0),
                    legendgroup='forward'
                ))
            
            # Reverse strand
            reverse_anns = [ann for ann in chr_annotations if ann['strand'] == '-']
            if reverse_anns:
                fig.add_trace(go.Scatter(
                    x=[ann['start'] for ann in reverse_anns],
                    y=[chr_name] * len(reverse_anns),
                    mode='markers',
                    marker=dict(
                        color='blue',
                        size=10,
                        symbol='triangle-left',
                        line=dict(width=1, color='darkblue')
                    ),
                    name='Reverse (-)' if i == 0 else '',
                    text=[f"{ann['motif']} ({ann['start']:,}-{ann['end']:,})" for ann in reverse_anns],
                    hovertemplate='<b>%{y}</b><br>Position: %{x:,}<br>%{text}<br>Strand: -<extra></extra>',
                    showlegend=(i == 0),
                    legendgroup='reverse'
                ))
        
        fig.update_layout(
            title=dict(
                text="Interactive Telomere Distribution",
                font=dict(size=18)
            ),
            xaxis=dict(
                title="Genomic Position (bp)",
                tickformat=',d'
            ),
            yaxis=dict(
                title="Chromosome"
            ),
            height=max(500, len(chr_data) * 60),
            hovermode='closest',
            template='plotly_white',
            showlegend=True
        )
        
        plot(fig, filename=output_file, auto_open=False)
        print(f"   ✓ Interactive plot saved: {output_file}")
        
    except ImportError:
        print("❌ Plotly not available. Install with: pip install plotly")

def main():
    parser = argparse.ArgumentParser(
        description='Visualize TeloX telomere annotations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python3 visualize_telomeres.py telox_genome_telomeres.bed
    python3 visualize_telomeres.py telox_genome_annotations.tsv --format tsv
    python3 visualize_telomeres.py telox_genome_telomeres.bed --interactive
    python3 visualize_telomeres.py telox_genome_telomeres.bed --text-only
        """
    )
    
    parser.add_argument('input_file', help='Input file (BED or TSV format)')
    parser.add_argument('--format', choices=['bed', 'tsv'], default='bed', 
                       help='Input file format (default: bed)')
    parser.add_argument('--text-only', action='store_true', 
                       help='Generate text summary only (no plots)')
    parser.add_argument('--interactive', action='store_true',
                       help='Generate interactive HTML plot (requires plotly)')
    parser.add_argument('--output-prefix', default='telomere',
                       help='Output file prefix (default: telomere)')
    parser.add_argument('--max-annotations', type=int, default=10000,
                       help='Maximum annotations to plot (default: 10000)')
    parser.add_argument('--max-chromosomes', type=int, default=50,
                       help='Maximum chromosomes to show (default: 50)')
    
    args = parser.parse_args()
    
    # Check if file exists
    try:
        # Load data
        if args.format == 'bed':
            annotations = load_bed_file(args.input_file)
        else:
            annotations = load_tsv_file(args.input_file)
        
        print(f"📁 Loaded {len(annotations)} telomere annotations from {args.input_file}")
        
        if not annotations:
            print("❌ No annotations found. Check file format and content.")
            print("\nExpected BED format:")
            print("   chromosome\\tstart\\tend\\tmotif\\tstrand\\tcolor")
            print("   CM060242.1\\t169\\t174\\tGGTGT\\t+\\t255,0,0")
            return
        
        # Generate visualizations
        create_text_summary(annotations)
        create_ascii_ideogram(annotations)
        
        if not args.text_only:
            # Apply memory limits
            plot_annotations = annotations
            if len(annotations) > args.max_annotations:
                print(f"⚠️  Large dataset: limiting to {args.max_annotations} annotations for plotting")
                # Sort by score and take top annotations
                sorted_annotations = sorted(annotations, key=lambda x: x.get('score', 0.5), reverse=True)
                plot_annotations = sorted_annotations[:args.max_annotations]
            
            create_matplotlib_plots(plot_annotations, args.output_prefix)
            
            if args.interactive:
                # Use smaller subset for interactive plots
                interactive_annotations = plot_annotations[:2000] if len(plot_annotations) > 2000 else plot_annotations
                if len(plot_annotations) > 2000:
                    print(f"⚠️  Using top 2000 annotations for interactive plot to prevent memory issues")
                create_interactive_html(interactive_annotations, f"{args.output_prefix}_interactive.html")
        
        print(f"\n✅ Visualization complete!")
        
    except FileNotFoundError:
        print(f"❌ Error: File '{args.input_file}' not found.")
        print("   Make sure to run TeloX with --annotate first:")
        print("   ./target/release/telox genome.fasta --annotate")
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    main()