# TeloX

**TeloX** (Telomere Motif Extraction Tool) is a high-performance bioinformatics tool for analyzing telomere motifs and k-mers in DNA sequences. It provides comprehensive analysis of strand bias, motif distribution, and telomere identification with optimized performance using parallel processing and efficient algorithms.

## Features

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

## Examples

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
