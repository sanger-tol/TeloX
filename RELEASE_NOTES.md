# TeloX Release Notes

## Version 1.0.0 - Initial Release (January 2025)

### 🎉 Core Features

#### **Two-Stage Telomere Discovery Pipeline**
- Intelligent two-stage approach for comprehensive telomere identification
- **Stage 1**: Searches 23 predefined telomere motifs using optimized `telo_finder` algorithm
- **Stage 2**: If no predefined motifs found, automatically launches k-mer discovery pipeline
- **Result**: Dramatically improved success rate for both known and novel telomere sequences

#### **Comprehensive K-mer Discovery Range**
- K-mer analysis covers **5-15 nucleotides** for comprehensive motif detection
- **Benefit**: Captures shorter telomeric motifs that were previously missed
- **Impact**: Better detection of species with compact telomere repeat units

#### **Tandem Repeat K-mer Filtering**
- Advanced filter prevents composite k-mers from inflating results
- **Problem Solved**: Eliminates cases where `AACCGAACCG` (2×`AACCG`) gets higher scores than the true `AACCG` motif
- **Algorithm**: Detects tandem repeats of 2 to k/2 length within k-mers ≥6bp
- **Result**: More accurate identification of fundamental telomere repeat units

#### **Rotational Motif Consolidation**
- Intelligent grouping of rotational variants
- **Example**: `TTAGGG`, `TAGGGT`, `AGGGTT`, `GGGTTA` → consolidated as single motif group
- **Statistics**: Combined frequency counts and maximum stretch values across variants
- **Output**: Clean, non-redundant motif identification

#### **Final Motif Identification with JSON Output**
- `telomere_motif_final.json` - structured results for the most likely telomere motif
- **Contains**: Canonical sequence, all rotational variants, comprehensive statistics
- **Rationale**: Detailed explanation of selection criteria and confidence metrics
- **Format**: Machine-readable JSON for downstream analysis and integration

### 🔧 Technical Improvements

#### **Indel-Tolerant Stretch Analysis**
- **ENHANCED**: Configurable indel tolerance for longest stretch calculations
- **Parameters**: `--max-indels` (default: 5) and `--max-gap-size` (default: 5)
- **Modes**: Default indel-tolerant mode vs `--strict` exact matching
- **Biological Relevance**: Accounts for sequencing errors and natural variation

#### **Comprehensive Debug Output**
- **NEW**: K-mer count tables for each size (5-15 mers)
- **NEW**: Top longest stretch results with debugging information
- **NEW**: Processing statistics and filtering summaries
- **Benefit**: Transparent analysis process for troubleshooting and validation

#### **Enhanced Biological Filtering**
- **EXPANDED**: 7 comprehensive biological filters (previously 6)
- **NEW Filter**: Tandem repeat k-mer detection and removal
- **Improved**: More stringent quality control for telomere motif candidates
- **Result**: Higher precision in motif identification

#### **Optimized Performance**
- **FOCUS**: Processes only scaffolds ≥1MB for computational efficiency
- **REGION**: Analyzes last 5000bp where telomeres are typically located
- **PARALLEL**: Multi-threaded processing using Rayon for faster analysis
- **MEMORY**: Optimized data structures for large genome processing

### 📊 Output Enhancements

#### **Consolidated Results Table**
```
Canonical_kmer  kmer_group                     forward_total   reverse_total   longeststretch
AGGGTT          TTAGGG,TAGGGT,AGGGTT,GGGTTA  1250            890             15
AACCCT          CCCTAA,CCTAAC                 456             234             12
```

#### **Enhanced File Structure**
- `initial_anno.txt` - Results from predefined motif search (Stage 1)
- `anno.txt` - Results from discovered motifs (Stage 2)
- `strand_bias_Nmer.tsv` - Strand bias analysis for k-mers (N=5-15)
- `rank.tsv` - All filtered k-mers ranked by stretch and frequency
- `telomere_motif_final.json` - **NEW**: Most likely telomere motif with full statistics

#### **Improved Ranking System**
- **Primary**: Longest continuous stretch (with indel tolerance)
- **Secondary**: Total frequency (forward + reverse counts)
- **Quality Filters**: Stretch ≥2, significance ≠ "weak"
- **Result**: More biologically relevant motif prioritization


### 📈 Performance Features

- **Optimized** k-mer processing through efficient filtering pipeline
- **Memory efficient** algorithms for large genome analysis
- **Scalable** design for genomes with many scaffolds
- **Parallel processing** across multiple k-mer sizes using Rayon

### 🧪 Algorithm Details

#### **Telomere Detection Scoring**
- Match: +1 point
- Mismatch: -1 point  
- Max drop: 2000 points
- Min score: 300 points
- **Result**: More sensitive detection of telomeric regions

#### **Strand Bias Classification**
- Forward biased: ratio >2.0
- Reverse biased: ratio <0.5
- Balanced: 0.5 ≤ ratio ≤ 2.0
- Strong significance: ratio >1.5 or <0.66
- **Result**: More accurate biological interpretation

### 📚 Documentation Updates

- **COMPREHENSIVE**: Completely rewritten README with detailed examples
- **NEW**: Algorithm overview and two-stage workflow explanation
- **ENHANCED**: Usage examples for different scenarios
- **ADDED**: Troubleshooting guide and best practices
- **UPDATED**: Logic diagram reflecting all new features

### 🔮 Future Compatibility

- **Rust 1.70+**: Minimum supported Rust version
- **JSON Output**: Structured format for integration with other tools
- **Extensible**: Modular design for future enhancements
- **Standards**: Follows bioinformatics best practices for reproducibility

### 🚀 Getting Started

```bash
# Standard de novo telomere discovery
telox genome.fasta

# With custom indel tolerance
telox genome.fasta --max-indels 3 --max-gap-size 3

# Strict mode for highly conserved sequences  
telox genome.fasta --strict

# Extract telomere-enriched regions
telox genome.fasta extract_lastN 10000 telomere_regions.fasta
```

### 🙏 Acknowledgments

Special thanks to the Wellcome Sanger Institute Tree of Life Programme for supporting this development, and to the bioinformatics community for valuable feedback and testing.

### 📞 Support

- **Issues**: Report bugs and feature requests on GitHub
- **Documentation**: See updated README.md and logic diagram
- **Examples**: Check the examples/ directory for usage scenarios

---

**Authors**: Yumi Sims, Chenxi Zhou, and Will Eagles  
**Institution**: Wellcome Sanger Institute  
**Date**: January 2025  
**License**: MIT License
