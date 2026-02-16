WGS - Advanced Gene Knockout Confirmation with SAMtools Flag Filtering
================
Neha Sahu
16 February, 2026

- [Advanced Gene Knockout Confirmation Workflow Using Whole Genome
  Sequencing](#advanced-gene-knockout-confirmation-workflow-using-whole-genome-sequencing)
  - [Overview](#overview)
  - [Prerequisites](#prerequisites)
- [Workflow Steps](#workflow-steps)
  - [Step 1: Data Upload and Initial
    Processing](#step-1-data-upload-and-initial-processing)
  - [Step 2: Advanced FASTQ to Multiple BAM
    Conversion](#step-2-advanced-fastq-to-multiple-bam-conversion)
    - [Script Configuration](#script-configuration)
  - [Step 3: Understanding SAMtools Flags and Filtering
    Categories](#step-3-understanding-samtools-flags-and-filtering-categories)
    - [SAM Flag Reference](#sam-flag-reference)
    - [The Nine Samtools Flags for BAM File
      Categories](#the-nine-samtools-flags-for-bam-file-categories)
  - [Step 4: Knockout Analysis
    Strategy](#step-4-knockout-analysis-strategy)
    - [Automated Region Analysis](#automated-region-analysis)
  - [Step 5: IGV Visualization
    Strategy](#step-5-igv-visualization-strategy)
    - [Recommended BAM File Loading
      Order](#recommended-bam-file-loading-order)
  - [Step 6: Quality Control Metrics](#step-6-quality-control-metrics)
    - [Key Statistics](#key-statistics)
  - [Other Useful Info](#other-useful-info)
    - [Contamination Detection](#contamination-detection)
  - [Expected Output Files](#expected-output-files)
    - [Primary BAM Files](#primary-bam-files)
    - [Statistics Files](#statistics-files)
    - [Storage Recommendations](#storage-recommendations)
  - [Contact and Support](#contact-and-support)

# Advanced Gene Knockout Confirmation Workflow Using Whole Genome Sequencing

## Overview

This advanced workflow describes the process for confirming gene
knockouts in *Magnaporthe oryzae* using whole genome sequencing (WGS)
data with sophisticated SAMtools flag filtering. Unlike the basic
approach, this pipeline creates multiple filtered BAM files to provide
comprehensive analysis of different read categories, enabling more
precise knockout confirmation and quality assessment.

The pipeline processes raw FASTQ files through to multiple filtered BAM
files, each targeting specific read characteristics that are crucial for
gene deletion validation.

## Prerequisites

- **Access to HPC** - Norwich Biosciences Institute/The Sainsbury
  Laboratory HPC system or any other
- **Raw WGS sequencing data** (paired-end FASTQ files) - from sequencing
  provider (e.g., Novogene)
- ***Magnaporthe oryzae*** **reference genome** (MGGv8_genome.fasta)
  from [Ensembl
  fungi](https://fungi.ensembl.org/Magnaporthe_oryzae/Info/Index)
- ***Magnaporthe oryzae*** **GTF/GFF3 annotation file** - from [Ensembl
  fungi](https://fungi.ensembl.org/info/data/ftp/index.html)
- **IGV software** installed locally
- **Understanding of SAM/BAM flags** for proper interpretation

# Workflow Steps

## Step 1: Data Upload and Initial Processing

1.  **Upload sequencing data**: Raw sequencing results must be uploaded
    to <https://sequences.tsl.ac.uk/>
2.  **Wait for processing**: Only proceed with downstream analysis after
    data has been fully processed (mandatory for projects run at TSL)
3.  **Access HPC**: Log into the Norwich Biosciences Institute/The
    Sainsbury Laboratory HPC system

## Step 2: Advanced FASTQ to Multiple BAM Conversion

The advanced script creates multiple filtered BAM files for
comprehensive analysis, using various SAMtools flags to separate
different read categories.

### Script Configuration

Run the script using `sbatch 02_advanced_fastq_to_bam_with_flags.sh`

⚠️ **Before running the script, ensure the following paths are correctly
configured:**

``` bash
#!/bin/bash
#SBATCH -p tsl-medium
#SBATCH --wckey=talbot_core
#SBATCH -N 1
#SBATCH --cpus-per-task=28
#SBATCH --mem 100000
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=neha.sahu@tsl.ac.uk

source package /tgac/software/production/bin/bowtie-2.2.5
source package c92263ec-95e5-43eb-a527-8f1496d56f1a  # Samtools 1.18

# Setup paths - EDIT ACCORDINGLY
reference_genome="MGGv8_genome.fasta"
index_prefix="MGGv8_index"
input_dir1="/path/to/your/sample1/"
input_dir2="/path/to/your/sample2/" #for multiple samples
output_dir="./alignments_wgs"

# Define knockout gene regions for analysis - this is Optional
declare -A GENES
GENES[mst12]="3:4220800-4225271"
GENES[bip1]="2:5750363-5753268"
# Add more genes as needed
```

## Step 3: Understanding SAMtools Flags and Filtering Categories

This advanced script creates nine different BAM files for each sample,
each filtered based on specific SAM flags. Understanding these flags is
crucial for proper knockout analysis. It is not mandatory to run all the
eight, but having them available can be very helpful during the
validation and troubleshooting process.

### SAM Flag Reference

SAM flags are bit-wise flags that encode information about each read.
More details can be found at
<https://samtools.github.io/hts-specs/SAMv1.pdf> (Page 7).

Here are the key flags used in this pipeline:

| Flag  | Hex  | Description   | Meaning                                      |
|-------|------|---------------|----------------------------------------------|
| 0x1   | 1    | PAIRED        | Read is paired                               |
| 0x2   | 2    | PROPER_PAIR   | Read mapped in proper pair (both read files) |
| 0x4   | 4    | UNMAP         | Query sequence is unmapped                   |
| 0x8   | 8    | MUNMAP        | Mate is unmapped                             |
| 0x10  | 16   | REVERSE       | Query on reverse strand                      |
| 0x20  | 32   | MREVERSE      | Mate on reverse strand                       |
| 0x40  | 64   | READ1         | First read in pair                           |
| 0x80  | 128  | READ2         | Second read in pair                          |
| 0x100 | 256  | SECONDARY     | Secondary alignment                          |
| 0x200 | 512  | QCFAIL        | QC failure                                   |
| 0x400 | 1024 | DUP           | PCR or optical duplicate                     |
| 0x800 | 2048 | SUPPLEMENTARY | Supplementary alignment                      |

### The Nine Samtools Flags for BAM File Categories

#### 1. All Reads (`*_all.bam`)

``` bash
samtools view -@ 28 -bS "$sam_file" | samtools sort -@ 28 -o "$bam_all" -
```

- **Filter**: None applied
- **Purpose**: Contains every read from the alignment
- **Use case**: Overall coverage assessment and general visualization
- **Includes**: Everything, all mapped, unmapped, primary, secondary,
  and supplementary alignments

#### 2. Unique Reads (`*_unique.bam`)

``` bash
samtools view -@ 28 -bS -q 10 -F 4 -F 256 -F 2048 "$sam_file"
```

- **Filters Applied**:
  - `-q 10`: MAPQ ≥ 10 (mapping quality score ≥ 10)
  - `-F 4`: Exclude unmapped reads (NOT flag 0x4)
  - `-F 256`: Exclude secondary alignments (NOT flag 0x100)
  - `-F 2048`: Exclude supplementary alignments (NOT flag 0x800)
- **Purpose**: High-confidence, uniquely mapped reads only
- **Use case**: High-quality deletion detection, eliminates ambiguous
  mappings
- **Quality**: High specificity for deletion detection

#### 3. Properly Paired Reads (`*_properly_paired.bam`)

``` bash
samtools view -@ 28 -bS -f 2 "$sam_file"
```

- **Filter Applied**: `-f 2` (REQUIRE flag 0x2)
- **Purpose**: Reads where both mates mapped at expected distance and
  orientation
- **Use case**: Structural variant detection and quality control
- **Includes**: Only reads flagged as properly paired by the aligner

#### 4. Unique + Properly Paired Reads (`*_unique_paired.bam`)

``` bash
samtools view -@ 28 -bS -q 10 -f 2 -F 4 -F 256 -F 2048 "$sam_file"
```

- **Filters Applied**:
  - `-q 10`: MAPQ ≥ 10 (mapping quality score ≥ 10)
  - `-f 2`: REQUIRE properly paired (flag 0x2)
  - `-F 4`: Exclude unmapped reads (NOT flag 0x4)
  - `-F 256`: Exclude secondary alignments (NOT flag 0x100)
  - `-F 2048`: Exclude supplementary alignments (NOT flag 0x800)
- **Purpose**: Combines unique mapping with proper pairing in both read
  files
- **Use case**: Highest confidence deletion detection
- **Quality**: Maximum specificity and reliability for knockout
  confirmation

#### 5. First in Pair (`*_first_in_pair.bam`)

``` bash
samtools view -@ 28 -bS -f 64 "$sam_file"
```

- **Filter Applied**: `-f 64` (REQUIRE flag 0x40)
- **Purpose**: Only the first read (R1) from each paired-end fragment
- **Use case**: Strand-specific analysis and coverage uniformity
  assessment

#### 6. Second in Pair (`*_second_in_pair.bam`)

``` bash
samtools view -@ 28 -bS -f 128 "$sam_file"
```

- **Filter Applied**: `-f 128` (REQUIRE flag 0x80)
- **Purpose**: Only the second read (R2) from each paired-end fragment
- **Use case**: Complement to first-in-pair for comprehensive paired-end
  analysis

#### 7. Unmapped Reads (`*_unmapped.bam`)

``` bash
samtools view -@ 28 -bS -f 4 "$sam_file"
```

- **Filter Applied**: `-f 4` (REQUIRE flag 0x4)
- **Purpose**: Reads that could not be mapped to the reference genome
- **Use case**: Quality control and identification of potential
  contamination or novel sequences

#### 8. Secondary Alignments (`*_secondary.bam`)

``` bash
samtools view -@ 28 -bS -f 256 "$sam_file"
```

- **Filter Applied**: `-f 256` (REQUIRE flag 0x100)
- **Purpose**: Alternative alignments for reads that map to multiple
  locations
- **Use case**: Identifying repetitive regions and multi-mapping reads
- **Note**: Bowtie2 reports one primary + multiple secondary alignments
  for multi-mappers

#### 9. Supplementary Alignments (`*_supplementary.bam`)

``` bash
samtools view -@ 28 -bS -f 2048 "$sam_file"
```

- **Filter Applied**: `-f 2048` (REQUIRE flag 0x800)
- **Purpose**: Chimeric alignments where parts of a read map to
  different locations
- **Use case**: Structural variant detection, including large
  insertions/deletions

## Step 4: Knockout Analysis Strategy

### Automated Region Analysis

The script automatically analyzes predefined gene regions:

``` bash
# Example knockout analysis output
Knockout region analysis:
  mst12:
    ALL: 1247 | UNIQUE: 856 | PAIRED: 1180 | UNIQUE+PAIRED: 834 | SECONDARY: 234
  bip1:
    ALL: 0 | UNIQUE: 0 | PAIRED: 0 | UNIQUE+PAIRED: 0 | SECONDARY: 0
```

## Step 5: IGV Visualization Strategy

### Recommended BAM File Loading Order

For checking gene deletions in IGV, load BAM and BAI files in this
priority order:

1.  **`*_unique_paired.bam`**

    - **Most stringent** for knockout confirmation
    - Shows only high-confidence, uniquely mapped AND properly paired
      reads

2.  **`*_unique.bam`**

    - High-confidence, uniquely mapped reads (less stringent than \#1)
    - Good for secondary confirmation

3.  **`*_all.bam`**

    - Complete coverage picture
    - Useful for overall context

4.  **`*_properly_paired.bam`**

    - Validates structural integrity
    - Confirms proper paired-end mapping

5.  **`*_secondary.bam` and `*_supplementary.bam`**

    - Load only if investigating further complex rearrangements
    - May show junction reads or partial deletions

## Step 6: Quality Control Metrics

### Key Statistics

The script generates detailed statistics for each BAM category:

``` bash
samtools flagstat "$bam_unique" > "${sample_name}_unique_stats.txt"
```

**Critical QC Metrics:**

\- **Total reads processed**: Should be consistent across all categories

\- **Mapping rate**: Check `*_all_stats.txt` for overall mapping
efficiency

\- **Properly paired rate**: Indicates library quality

\- **Secondary alignment rate**: High rates may indicate repetitive
regions

## Other Useful Info

### Contamination Detection

Use `*_unmapped.bam` to identify:

\- Non-target organism contamination

\- Adapter sequences (if any)

\- Novel sequences not in reference

## Expected Output Files

For each sample, the pipeline generates:

### Primary BAM Files

- `{sample}_all.bam` + `.bai`
- `{sample}_unique.bam` + `.bai`
- `{sample}_properly_paired.bam` + `.bai`
- `{sample}_unique_paired.bam` + `.bai`
- `{sample}_first_in_pair.bam` + `.bai`
- `{sample}_second_in_pair.bam` + `.bai`
- `{sample}_unmapped.bam` + `.bai`
- `{sample}_secondary.bam` + `.bai`
- `{sample}_supplementary.bam` + `.bai`

### Statistics Files

- `{sample}_all_stats.txt`
- `{sample}_unique_stats.txt`
- `{sample}_paired_stats.txt`
- `{sample}_unique_paired_stats.txt`

### Storage Recommendations

**Essential Files (Keep):**

\- `*_unique_paired.bam` and `*_unique_paired.bam.bai`

\- `*_unique.bam` and `*_unique.bam.bai`

\- `*_all.bam` and `*_all.bam.bai`

\- `*_unique_paired_stats.txt`, `*_unique_stats.txt` and
`*_all_stats.txt`

**Optional Files (Archive after analysis):** - All other filtered BAM
files as well as intermediate SAM files as these are uncompressed and
can be huge.

## Contact and Support

For questions about this advanced workflow or issues with flag
interpretation, contact: <neha.sahu.tsl@gmail.com>

------------------------------------------------------------------------

*This enhanced workflow was developed for the most precise gene knockout
confirmation using comprehensive SAMtools flag filtering in Magnaporthe
oryzae research projects. The nine-BAM approach with unique+properly
paired filtering provides maximum sensitivity and specificity for
deletion detection.*
