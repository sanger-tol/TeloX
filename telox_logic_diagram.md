# TeloX - Telomere Motif Extraction Tool Logic Diagram

## Program Overview
TeloX is a Rust-based tool for identifying telomere motifs in genomic sequences through a two-stage approach: predefined motif searching and k-mer discovery.

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              TELOX MAIN WORKFLOW                                │
└─────────────────────────────────────────────────────────────────────────────────┘

INPUT: FASTA File
        │
        ▼
┌─────────────────────┐
│  Command Line Args  │
│  Processing         │
└─────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           EXECUTION MODES                                      │
├─────────────────────┬─────────────────────┬─────────────────────────────────────┤
│  KMC Pipeline Mode  │  Extract LastN Mode │    Default Analysis Mode          │
│  (kmc command)      │  (extract_lastN)    │    (main telomere analysis)       │
└─────────────────────┴─────────────────────┴─────────────────────────────────────┘
                                                        │
                                                        ▼
                                            ┌─────────────────────┐
                                            │  STEP 1: Initial   │
                                            │  Predefined Motif  │
                                            │  Search            │
                                            └─────────────────────┘
                                                        │
                                                        ▼
                                            ┌─────────────────────┐
                                            │  Load FASTA file    │
                                            │  using SequenceDict │
                                            └─────────────────────┘
                                                        │
                                                        ▼
                                            ┌─────────────────────┐
                                            │  Run telo_finder    │
                                            │  with TELO_MOTIF_DB │
                                            │  (23 predefined     │
                                            │   motifs)           │
                                            └─────────────────────┘
                                                        │
                                                        ▼
                                            ┌─────────────────────┐
                                            │  Check if any       │
                                            │  telomeres found?   │
                                            └─────────────────────┘
                                                    │
                                            ┌───────┴───────┐
                                            │ YES           │ NO
                                            ▼               ▼
                                ┌─────────────────────┐  ┌─────────────────────┐
                                │  Generate Report    │  │  STEP 2: K-mer     │
                                │  & Exit Successfully│  │  Discovery Pipeline │
                                └─────────────────────┘  └─────────────────────┘
                                                                    │
                                                                    ▼
                                                        ┌─────────────────────┐
                                                        │  For k = 5 to 15:  │
                                                        │  Process k-mers     │
                                                        └─────────────────────┘
                                                                    │
                                                                    ▼
                                                        ┌─────────────────────┐
                                                        │  Count k-mers in    │
                                                        │  last 5000bp of     │
                                                        │  scaffolds ≥1MB     │
                                                        └─────────────────────┘
                                                                    │
                                                                    ▼
                                                        ┌─────────────────────┐
                                                        │  Apply Biological   │
                                                        │  Filters            │
                                                        └─────────────────────┘
                                                                    │
                                                                    ▼
                                                        ┌─────────────────────┐
                                                        │  Analyze Strand     │
                                                        │  Bias               │
                                                        └─────────────────────┘
                                                                    │
                                                                    ▼
                                                        ┌─────────────────────┐
                                                        │  Calculate Longest  │
                                                        │  Continuous Stretch │
                                                        └─────────────────────┘
                                                                    │
                                                                    ▼
                                                        ┌─────────────────────┐
                                                        │  Filter & Rank      │
                                                        │  by stretch & count │
                                                        └─────────────────────┘
                                                                    │
                                                                    ▼
                                                        ┌─────────────────────┐
                                                        │  STEP 3: Generate   │
                                                        │  New Motif Database │
                                                        └─────────────────────┘
                                                                    │
                                                                    ▼
                                                        ┌─────────────────────┐
                                                        │  Consolidate by     │
                                                        │  Rotation &         │
                                                        │  Generate Final     │
                                                        │  Motif JSON         │
                                                        └─────────────────────┘
                                                                    │
                                                                    ▼
                                                        ┌─────────────────────┐
                                                        │  Run telo_finder    │
                                                        │  with discovered    │
                                                        │  motifs             │
                                                        └─────────────────────┘
                                                                    │
                                                                    ▼
                                                        ┌─────────────────────┐
                                                        │  Generate Final     │
                                                        │  Report             │
                                                        └─────────────────────┘

```

## Detailed Component Breakdown

### 1. TELO_FINDER Core Algorithm

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           TELO_FINDER_CORE ALGORITHM                           │
└─────────────────────────────────────────────────────────────────────────────────┘

Input: DNA sequence, motif table, scoring parameters
        │
        ▼
┌─────────────────────┐
│  Generate Motif     │
│  Hash Table         │
│  (all rotations)    │
└─────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                          SCORING ALGORITHM                                     │
├─────────────────────────────────────────────────────────────────────────────────┤
│  For each position in sequence:                                                │
│    1. Convert nucleotide to 2-bit encoding                                     │
│    2. Slide k-mer window, check against motif table                           │
│    3. Update score: +1 for match, -penalty for mismatch                       │
│    4. Track maximum score and position                                         │
│    5. Stop if score drops > max_drop from maximum                             │
└─────────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────┐     ┌─────────────────────┐
│  5' End Analysis    │     │  3' End Analysis    │
│  (forward strand)   │     │  (reverse strand)   │
└─────────────────────┘     └─────────────────────┘
        │                           │
        ▼                           ▼
┌─────────────────────┐     ┌─────────────────────┐
│  Score ≥ min_score? │     │  Score ≥ min_score? │
│  Mark as telomere   │     │  Mark as telomere   │
└─────────────────────┘     └─────────────────────┘
        │                           │
        └─────────────┬─────────────┘
                      ▼
            ┌─────────────────────┐
            │  Return telomere    │
            │  positions (5', 3') │
            └─────────────────────┘

Constants:
- TELO_PENALTY = 1
- TELO_MAX_DROP = 2000  
- TELO_MIN_SCORE = 300
```

### 2. K-mer Discovery Pipeline

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           K-MER DISCOVERY PIPELINE                             │
└─────────────────────────────────────────────────────────────────────────────────┘

Input: FASTA file
        │
        ▼
┌─────────────────────┐
│  Filter Scaffolds   │
│  (≥1MB only)        │
└─────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Extract Last 5000bp│
│  from each scaffold │
└─────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Generate k-mers    │
│  (k = 5 to 15)      │
└─────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                          BIOLOGICAL FILTERS                                    │
├─────────────────────────────────────────────────────────────────────────────────┤
│  1. Remove k-mers containing 'N'                                               │
│  2. Remove homopolymers (AAAA, TTTT, etc.)                                     │
│  3. Remove dinucleotide repeats (ATATAT, CGCGCG, etc.)                        │
│  4. Remove simple repeats (tandem repeats of 1-4 bp)                          │
│  5. NEW: Remove tandem repeat k-mers (e.g., AACCGAACCG = 2×AACCG)             │
│  6. Remove low G/C content (<28% G or C on either strand)                     │
│  7. Remove k-mers with whitespace                                              │
└─────────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Count Forward &    │
│  Reverse Complement │
│  occurrences        │
└─────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                          STRAND BIAS ANALYSIS                                  │
├─────────────────────────────────────────────────────────────────────────────────┤
│  For each k-mer:                                                               │
│    - Calculate bias_ratio = forward_count / rc_count                          │
│    - Direction: forward (>2.0), reverse (<0.5), balanced (0.5-2.0)           │
│    - Significance: strong (>1.5 or <0.66), weak (0.66-1.5)                   │
└─────────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Calculate Longest  │
│  Continuous Stretch │
│  for each k-mer     │
└─────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Filter by:         │
│  - stretch ≥ 2      │
│  - significance ≠   │
│    "weak"           │
└─────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Rank by:           │
│  1. Longest stretch │
│     (descending)    │
│  2. Total count     │
│     (descending)    │
└─────────────────────┘
```

### 3. Data Structures and Key Functions

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           KEY DATA STRUCTURES                                  │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  SequenceDict                    KmerPair                                       │
│  ├── sequences: Vec<Sequence>    ├── forward: u32                              │
│                                  └── rc: u32                                    │
│  Sequence                                                                       │
│  ├── name: String                StrandBiasAnalysis                            │
│  └── seq: String                 ├── kmer: String                              │
│                                  ├── forward_count: u32                        │
│  TelomereCandidate              ├── rc_count: u32                              │
│  ├── status: String             ├── total_count: u32                           │
│  ├── primary_motif: Option       ├── bias_ratio: f64                           │
│  ├── frequency: Option           ├── bias_direction: String                    │
│  ├── percentage: Option          ├── significance: String                      │
│  ├── source: Option              └── longest_stretch: usize                    │
│  ├── total_occurrences: Option                                                 │
│  ├── unique_motifs: Option                                                     │
│  ├── top_motifs: Vec<MotifFreq>                                                │
│  └── message: Option                                                           │
└─────────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────────┐
│                           KEY ALGORITHMS                                       │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  longest_continuous_stretch_for_kmers()                                        │
│  ├── Input: sequence, k-mer list                                               │
│  ├── Method: Single-pass sliding window                                        │
│  ├── Tracks: consecutive non-overlapping k-mer occurrences                     │
│  └── Output: HashMap<kmer, max_stretch>                                        │
│                                                                                 │
│  count_kmers_last_n_bp_parallel()                                              │
│  ├── Input: FASTA, k-size, bp_count                                           │
│  ├── Method: Parallel processing with rayon                                    │
│  ├── Filters: scaffolds ≥1MB only                                             │
│  └── Output: HashMap<kmer, KmerPair>                                           │
│                                                                                 │
│  consolidate_ranked_motifs_by_rotation()                                       │
│  ├── Input: ranked TSV file                                                    │
│  ├── Method: Minimal rotation canonicalization                                 │
│  ├── Purpose: Remove rotational duplicates                                     │
│  └── Output: Vec<unique_motifs>                                                │
└─────────────────────────────────────────────────────────────────────────────────┘
```

### 4. Output Files and Results

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              OUTPUT FILES                                      │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  initial_anno.txt                   strand_bias_Nmer.tsv (N=5-15)              │
│  ├── Results from predefined motifs ├── K-mer strand bias analysis             │
│  └── Format: seqname\tstart\tend\tmotif └── Columns: Kmer, Forward, RC, etc.   │
│                                                                                 │
│  anno.txt                           rank.tsv                                   │
│  ├── Results from discovered motifs ├── All filtered k-mers ranked             │
│  └── Same format as initial_anno.txt └── Sorted by stretch, then count         │
│                                                                                 │
│  telomere_candidate.json            telomere_motif_final.json                  │
│  ├── Primary motif identification   ├── Most likely telomere motif             │
│  ├── Structured JSON output         ├── Rotational variants & statistics       │
│  └── Machine-readable results       └── Consolidated analysis results          │
└─────────────────────────────────────────────────────────────────────────────────┘
```

### 5. Scoring and Parameter Details

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           SCORING PARAMETERS                                   │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  Telomere Detection:                Biological Filters:                        │
│  ├── Match: +1 point               ├── G/C content: >28%                       │
│  ├── Mismatch: -1 point            ├── No homopolymers                         │
│  ├── Max drop: 2000 points         ├── No dinuc repeats                        │
│  └── Min score: 300 points         ├── No simple repeats                       │
│                                     └── No ambiguous bases                      │
│  Strand Bias Classification:                                                   │
│  ├── Forward biased: ratio >2.0    Significance Levels:                        │
│  ├── Reverse biased: ratio <0.5    ├── Strong: ratio >1.5 or <0.66           │
│  ├── Balanced: 0.5 ≤ ratio ≤ 2.0   └── Weak: 0.66 ≤ ratio ≤ 1.5             │
│                                                                                 │
│  Quality Thresholds:                                                           │
│  ├── Scaffold size: ≥1MB                                                      │
│  ├── Analysis region: last 5000bp                                             │
│  ├── Min stretch: ≥2 consecutive                                              │
│  └── Exclude weak significance                                                 │
└─────────────────────────────────────────────────────────────────────────────────┘
```

This diagram shows the complete logic flow of TeloX, from input processing through the two-stage analysis (predefined motifs → k-mer discovery) to final output generation. The tool is designed to be robust, using parallel processing for efficiency and multiple biological filters to ensure high-quality telomere motif identification.
