# TeloX

**TeloX** (Telomere Motif Extraction Tool) is a high-performance bioinformatics tool for **de novo identification and analysis of telomere motifs** in DNA sequences. It provides comprehensive analysis of strand bias, motif distribution, and telomere identification with optimized performance using parallel processing and efficient algorithms.

## Features

- **De novo telomere motif discovery** - Automatically identify novel telomere motifs in species genomes
- **High-performance k-mer counting** using KMC (K-Mer Counter) for fast processing of large genomes
- **Strand bias analysis** with statistical significance testing
- **Telomere motif identification** with customizable motif databases
- **Longest continuous stretch analysis** for motif clustering assessment
- **Parallel processing** using Rayon for multi-threaded performance
- **Memory-efficient algorithms** with optimized data structures
- **Comprehensive output formats** including TSV tables and summary statistics
- **Hybrid approach** combining KMC for counting and Rust for positional analysis

## Installation

### Prerequisites

- Rust (version 1.70 or higher)
- KMC (K-Mer Counter) - [Download from GitHub](https://github.com/refresh-bio/KMC)
- Make sure KMC binaries (`kmc` and `kmc_dump`) are in your PATH

### Building from Source

```bash
# Clone the repository
git clone https://github.com/your-username/telox.git
cd telox

# Build the project
cargo build --release

# Install globally (optional)
cargo install --path .
```

The compiled binary will be available at `target/release/telox`.

## Usage

### Basic Usage

```bash
# Analyze telomere motifs in a FASTA file
telox input.fasta

# Specify custom telomere motif
telox --telo-motif "TTAGGG" input.fasta

# List available telomere motifs
telox --print-telo-motifs
```

### Advanced Usage

```bash
# K-mer analysis with strand bias
telox --kmer-size 6 --strand-bias input.fasta

# Analyze last 5000bp of each scaffold for strand-specific counts
telox --last-5000bp --kmer-size 6 input.fasta

# Filter results by significance and stretch length
telox --filter-significant --min-stretch 2 input.fasta

# Output to specific files
telox --output-prefix results input.fasta
```

### Command Line Options

```
USAGE:
    telox [OPTIONS] <FASTA_FILE>

ARGS:
    <FASTA_FILE>    Input FASTA file to analyze

OPTIONS:
    -k, --kmer-size <SIZE>           K-mer size for analysis [default: 6]
    --telo-motif <MOTIF>             Custom telomere motif to search for
    --print-telo-motifs              Print available telomere motifs and exit
    --strand-bias                    Perform strand bias analysis
    --last-5000bp                    Analyze last 5000bp of each scaffold
    --filter-significant             Filter results by significance
    --min-stretch <LENGTH>           Minimum stretch length for filtering [default: 2]
    --output-prefix <PREFIX>         Output file prefix [default: telox]
    -h, --help                       Print help information
    -V, --version                    Print version information
```

## Output Files

TeloX generates several output files depending on the analysis performed:

- `{prefix}_strand_bias.tsv` - Strand bias analysis results
- `{prefix}_filtered.tsv` - Filtered results (significant motifs only)
- `{prefix}_summary.txt` - Summary statistics
- `{prefix}_last_5000bp.tsv` - Last 5000bp analysis results

### Output Format

The strand bias analysis includes:
- **kmer**: The k-mer sequence
- **forward_count**: Count in forward orientation
- **rc_count**: Count in reverse complement orientation
- **total_count**: Total count across both orientations
- **bias_ratio**: Ratio of forward to reverse complement counts
- **bias_direction**: Direction of bias (forward/reverse/balanced)
- **significance**: Statistical significance (strong/moderate/weak)
- **longest_stretch**: Longest continuous stretch of the motif

## Telomere Motif Database

TeloX includes a comprehensive database of known telomere motifs:

```
[ 1] AAAATTGTCCGTCC
[ 2] AAACCACCCT
[ 3] AAACCC
[ 4] AAACCCC
[ 5] AAACCCT
[ 6] AAAGAACCT
[ 7] AAATGTGGAGG
[ 8] AACAGACCCG
[ 9] AACCATCCCT
[10] AACCC
[11] AACCCAGACCC
[12] AACCCAGACCT
[13] AACCCAGACGC
[14] AACCCCAACCT
[15] AACCCGAACCT
[16] AACCCT
[17] AACCCTG
[18] AACCCTGACGC
[19] AACCT
[20] AAGGAC
[21] ACCCAG
[22] ACCTG
[23] ACGGCAGCG
```

## Performance Optimizations

TeloX implements several performance optimizations:

1. **KMC Integration**: Uses KMC for ultra-fast k-mer counting
2. **Parallel Processing**: Multi-threaded analysis using Rayon
3. **Memory Efficiency**: Optimized data structures and reduced allocations
4. **Fast Hashing**: Uses ahash for high-performance hash tables
5. **Byte-level Operations**: Efficient string processing with byte slices

## Algorithm Details

### Strand Bias Analysis

The strand bias analysis calculates:
- **Bias Ratio**: `forward_count / rc_count`
- **Significance**: Based on chi-square test with Bonferroni correction
- **Direction**: Determined by bias ratio thresholds

### Longest Continuous Stretch

For each k-mer, TeloX finds the longest continuous stretch of non-overlapping occurrences in the sequence, considering both forward and reverse complement orientations.

### Hybrid Approach

TeloX uses a hybrid approach combining:
- **KMC**: For fast k-mer counting across entire genomes
- **Rust**: For positional analysis and strand-specific counting in specific regions

## De Novo Telomere Motif Discovery

TeloX excels at **de novo identification** of telomere motifs in species genomes where the telomere sequence is unknown. This is particularly valuable for:

- **Non-model organisms** with unknown telomere sequences
- **Novel species** where telomere motifs haven't been characterized
- **Comparative genomics** studies across different taxa
- **Evolutionary studies** of telomere sequence diversity

### De Novo Discovery Workflow

```bash
# 1. Perform comprehensive k-mer analysis to identify candidate motifs
telox --kmer-size 6 --strand-bias genome.fasta

# 2. Filter for significant strand bias patterns
telox --kmer-size 6 --strand-bias --filter-significant genome.fasta

# 3. Analyze last 5000bp of scaffolds for telomere-enriched regions
telox --last-5000bp --kmer-size 6 genome.fasta

# 4. Generate filtered results for candidate telomere motifs
telox --filter-significant --min-stretch 2 genome.fasta
```

### Identifying Novel Telomere Motifs

TeloX identifies novel telomere motifs by analyzing:

1. **Strand Bias Patterns**: Telomere motifs typically show strong strand bias
2. **End Enrichment**: Motifs enriched in the last 5000bp of scaffolds
3. **Continuous Stretches**: Long continuous stretches of motif repeats
4. **Statistical Significance**: Statistically significant bias ratios
5. **Motif Conservation**: Consistent patterns across multiple scaffolds

### Example: Discovering Unknown Telomere Motifs

```bash
# For a species with unknown telomere sequence
telox --kmer-size 6 --strand-bias --last-5000bp --filter-significant genome.fasta

# This will output candidates like:
# kmer    forward_count  rc_count  bias_ratio  significance  longest_stretch
# TTAGGG  1250          45        27.78       strong        15
# AACCCT  42            1180      0.036       strong        12
# CCCTAA  1180          42        28.10       strong        12
```

The top candidates with high bias ratios and long stretches are likely the telomere motifs for that species.

## Examples

### De Novo Telomere Discovery

```bash
# Discover telomere motifs in a new species
telox --kmer-size 6 --strand-bias --last-5000bp --filter-significant new_species.fasta

# Analyze multiple k-mer sizes for comprehensive discovery
telox --kmer-size 4 --strand-bias --filter-significant new_species.fasta
telox --kmer-size 5 --strand-bias --filter-significant new_species.fasta
telox --kmer-size 6 --strand-bias --filter-significant new_species.fasta
telox --kmer-size 7 --strand-bias --filter-significant new_species.fasta
```

### Basic Telomere Analysis

```bash
# Analyze telomere motifs in a genome assembly
telox genome.fasta

# Use custom telomere motif
telox --telo-motif "TTAGGG" genome.fasta
```

### K-mer Analysis with Strand Bias

```bash
# Analyze 6-mers with strand bias
telox --kmer-size 6 --strand-bias genome.fasta

# Filter for significant results
telox --kmer-size 6 --strand-bias --filter-significant genome.fasta
```

### Last 5000bp Analysis

```bash
# Analyze last 5000bp of each scaffold
telox --last-5000bp --kmer-size 6 genome.fasta
```

## Limitations

- Requires KMC to be installed and available in PATH
- Memory usage scales with genome size and k-mer size
- Maximum k-mer size is limited by available memory
- KMC output is limited to 255 counts by default (can be increased with `-cs` parameter)

## Citation

If you use TeloX in your research, please cite:

```
Yumi Sims. TeloX: Telomere Motif Extraction Tool. 
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
