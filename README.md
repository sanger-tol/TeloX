# TeloX

**TeloX** (Telomere Motif Extraction Tool) is a high-performance bioinformatics tool for **de novo identification and analysis of telomere motifs** in DNA sequences. It provides comprehensive analysis of strand bias, motif distribution, and telomere identification with optimized performance using parallel processing and efficient algorithms.

## Features

- **Two-stage telomere discovery** - Predefined motif search followed by k-mer discovery pipeline
- **De novo telomere motif discovery** - Automatically identify novel telomere motifs in species genomes
- **Comprehensive k-mer analysis** (5-15 mers) with biological filtering
- **Strand bias analysis** with statistical significance testing
- **Longest continuous stretch analysis** with indel tolerance options
- **Tandem repeat filtering** - Prevents composite k-mers from inflating results
- **Rotational motif consolidation** - Groups rotational variants of the same motif
- **Parallel processing** using Rayon for multi-threaded performance
- **Memory-efficient algorithms** optimized for large genomes (≥1MB scaffolds)
- **Comprehensive output formats** including TSV tables, JSON results, and debug information
- **Intelligent motif ranking** based on stretch length and frequency

## Installation

### Prerequisites

- Rust (version 1.70 or higher)
- Optional: KMC (K-Mer Counter) for external k-mer counting - [Download from GitHub](https://github.com/refresh-bio/KMC)

### Building from Source

```bash
# Clone the repository
git clone https://github.com/your-username/telox.git
cd telox/telox

# Build the project
cargo build --release

# Install globally (optional)
cargo install --path .
```

The compiled binary will be available at `target/release/telox`.

## Usage

### Basic Usage

```bash
# Standard telomere analysis (two-stage approach)
telox genome.fasta

# With indel tolerance parameters
telox genome.fasta --max-indels 5 --max-gap-size 5

# Strict mode (exact matching only)
telox genome.fasta --strict
```

### Advanced Usage

```bash
# Extract last N bp of large scaffolds to separate FASTA
telox genome.fasta extract_lastN 10000 last10000.fasta

# KMC pipeline mode (external k-mer counting)
telox kmc input.fasta 7 kmc_db kmc_output.txt
```

### Command Line Options

```
USAGE:
    telox <FASTA_FILE> [OPTIONS]
    telox kmc <input_fasta> <k> <db_prefix> <output_txt>
    telox <fasta_file> extract_lastN <N> <output_fasta>

ARGS:
    <FASTA_FILE>    Input FASTA file to analyze

OPTIONS:
    --max-indels <N>      Maximum indels allowed in stretch analysis [default: 5]
    --max-gap-size <N>    Maximum gap size before resetting stretch [default: 5]
    --strict              Use exact matching only (no indel tolerance)
    -h, --help           Print help information
    -V, --version        Print version information

EXAMPLES:
    telox genome.fasta                           # Standard analysis with indel tolerance
    telox genome.fasta --max-indels 3 --max-gap-size 3  # More strict indel tolerance
    telox genome.fasta --strict                  # Exact matching only
    telox genome.fasta extract_lastN 10000 last10000.fasta
    telox kmc last5000.fasta 7 kmc_db kmc_dump.txt
```

## Algorithm Overview

### Two-Stage Approach

TeloX uses a sophisticated two-stage approach for telomere discovery:

**Stage 1: Predefined Motif Search**
- Searches for 23 known telomere motifs using the `telo_finder` algorithm
- Uses optimized scoring system (match: +1, mismatch: -1, max drop: 2000, min score: 300)
- If telomeres are found, analysis completes with results in `initial_anno.txt`

**Stage 2: K-mer Discovery Pipeline** (if Stage 1 finds no telomeres)
- Analyzes k-mers from 5-15 nucleotides in length
- Processes last 5000bp of scaffolds ≥1MB for telomere-enriched regions
- Applies comprehensive biological filters
- Performs strand bias analysis and longest stretch calculations
- Ranks and consolidates results by rotational equivalence
- Generates new motif database for final telomere search

### Biological Filtering

TeloX applies rigorous biological filters to ensure high-quality results:

1. **N-base removal** - Excludes k-mers containing ambiguous bases
2. **Homopolymer filtering** - Removes simple repeats (AAAA, TTTT, etc.)
3. **Dinucleotide repeat filtering** - Excludes alternating patterns (ATATAT, etc.)
4. **Simple repeat filtering** - Removes tandem repeats of 1-4bp units
5. **Tandem repeat k-mer filtering** - NEW: Prevents composite k-mers (e.g., AACCGAACCG = 2×AACCG)
6. **G/C content filtering** - Requires >28% G or C content on either strand
7. **Quality thresholds** - Longest stretch ≥2, significance ≠ "weak"

### Longest Stretch Analysis

Two modes available:

**Exact Matching** (--strict):
- Counts consecutive, non-overlapping k-mer occurrences
- Requires perfect sequence matches

**Indel-Tolerant** (default):
- Allows small gaps and indels in stretch calculation
- Configurable parameters: `--max-indels` and `--max-gap-size`
- More biologically relevant for real sequencing data

### Rotational Consolidation

Groups rotational variants of the same motif:
- Example: TTAGGG, TAGGGT, AGGGTT, GGGTTA → canonical AGGGTT
- Combines statistics: forward counts, reverse counts, max stretch
- Prevents redundant motif identification

## Output Files

TeloX generates comprehensive output files:

### Primary Output Files
- `initial_anno.txt` - Results from predefined motif search (Stage 1)
- `anno.txt` - Results from discovered motifs (Stage 2)
- `telomere_candidate.json` - Structured results for machine processing
- `telomere_motif_final.json` - Most likely telomere motif with statistics

### Analysis Files
- `strand_bias_Nmer.tsv` - Strand bias analysis for each k-mer size (N=5-15)
- `rank.tsv` - All filtered k-mers ranked by stretch and frequency

### Debug Files (when enabled)
- K-mer count tables for each size
- Top longest stretch results
- Processing statistics

### Output Format

The strand bias analysis includes:
- **Kmer**: The k-mer sequence
- **Forward**: Count in forward orientation
- **RC**: Count in reverse complement orientation  
- **Total**: Total count across both orientations
- **BiasRatio**: Ratio of forward to reverse complement counts
- **Direction**: Direction of bias (forward/reverse/balanced)
- **Significance**: Statistical significance (strong/weak)
- **LongestStretchIndels**: Longest continuous stretch (with indel tolerance)

### Consolidated Results Format

```
Canonical_kmer  kmer_group                     forward_total   reverse_total   longeststretch
AGGGTT          TTAGGG,TAGGGT,AGGGTT,GGGTTA  1250            890             15
AACCCT          CCCTAA,CCTAAC                 456             234             12
```

## Telomere Motif Database

TeloX includes a comprehensive database of 23 known telomere motifs:

```
[ 1] AAAATTGTCCGTCC    [13] AACCCAGACGC
[ 2] AAACCACCCT        [14] AACCCCAACCT
[ 3] AAACCC            [15] AACCCGAACCT
[ 4] AAACCCC           [16] AACCCT
[ 5] AAACCCT           [17] AACCCTG
[ 6] AAAGAACCT         [18] AACCCTGACGC
[ 7] AAATGTGGAGG       [19] AACCT
[ 8] AACAGACCCG        [20] AAGGAC
[ 9] AACCATCCCT        [21] ACCCAG
[10] AACCC             [22] ACCTG
[11] AACCCAGACCC       [23] ACGGCAGCG
[12] AACCCAGACCT
```

## De Novo Telomere Motif Discovery

TeloX excels at **de novo identification** of telomere motifs in species genomes where the telomere sequence is unknown. This is particularly valuable for:

- **Non-model organisms** with unknown telomere sequences
- **Novel species** where telomere motifs haven't been characterized
- **Comparative genomics** studies across different taxa
- **Evolutionary studies** of telomere sequence diversity

### Discovery Workflow

```bash
# Standard de novo discovery
telox unknown_species.fasta

# With custom indel tolerance
telox unknown_species.fasta --max-indels 3 --max-gap-size 3

# Strict mode for highly conserved sequences
telox unknown_species.fasta --strict
```

### Example Output

```bash
# For a species with unknown telomere sequence
telox new_species.fasta

# Console output shows:
Processing 5-mers...
Top 20 5-mers by count:
Rank      K-mer           Total Count     Forward         RC             
----------------------------------------------------------------------
1         AACCG           1250            800             450            

Top 10 5-mers by longest stretch:
  1: AACCG (stretch: 25, count: 1250)

# Final consolidated results:
Canonical_kmer  kmer_group     forward_total   reverse_total   longeststretch
AACCG           AACCG,CCGAA    1250            890             25

# JSON output in telomere_motif_final.json:
{
  "most_likely_telomere_motif": {
    "canonical_sequence": "AACCG",
    "rotational_variants": ["AACCG", "CCGAA"],
    "statistics": {
      "forward_count": 1250,
      "reverse_count": 890,
      "total_frequency": 2140,
      "longest_continuous_stretch": 25
    }
  }
}
```

## Examples

### Basic Telomere Analysis

```bash
# Analyze telomere motifs in a genome assembly
telox genome.fasta

# With custom parameters
telox genome.fasta --max-indels 3 --max-gap-size 2
```

### Strict Analysis

```bash
# Exact matching only (no indel tolerance)
telox genome.fasta --strict
```

### Utility Functions

```bash
# Extract last 10kb of large scaffolds
telox genome.fasta extract_lastN 10000 telomere_regions.fasta

# External KMC analysis
telox kmc genome.fasta 6 kmc_db kmc_output.txt
```

## Troubleshooting

### Common Issues

**No telomeres found**: 
- Try strict mode: `--strict`
- Check if scaffolds are ≥1MB
- Verify FASTA format

**Composite k-mers identified**:
- The tandem repeat filter should prevent this
- Check debug output for filtering statistics

**Memory issues**:
- Process smaller genome chunks
- Use external KMC mode for very large genomes

## Citation

If you use TeloX in your research, please cite:

```
Yumi Sims, Chenxi Zhou and Will Eagles TeloX: Telomere Motif Extraction Tool. 
Wellcome Sanger Institute, 2025.
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Issues

If you encounter any issues or have questions, please file an issue on the GitHub repository.

## About

TeloX is developed at the Wellcome Sanger Institute for high-performance telomere motif analysis in genomic sequences.