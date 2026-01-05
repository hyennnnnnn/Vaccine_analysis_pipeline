# GC Pipeline - Influenza Sequencing Data Analysis

A robust, resumable bioinformatics pipeline for analyzing influenza sequencing data with comprehensive progress tracking and error handling.

## Features

- **Resumability**: Automatically skip completed samples when re-running the pipeline
- **Progress Tracking**: Real-time progress updates with emoji indicators and timestamps
- **Comprehensive Logging**: Individual log files for each sample and analysis step
- **Error Handling**: Robust error detection and reporting with detailed failure messages
- **Summary Reports**: Automatic generation of pipeline execution summaries
- **Protein Analysis**: Consensus generation, ORF detection, translation, and multiple sequence alignment
- **Variant Impact**: Analyze amino acid changes from SNPs and INDELs in protein sequences

## Directory Structure

```
GC_analysis/
├── setup.sh                              # 1. Project setup script
├── GC_pipeline.sh                        # 2. Analysis pipeline script (sequential)
├── GC_pipeline_parallel.sh               # 2. Analysis pipeline script (parallel)
├── README.md                             # Documentation
├── .pipeline/                            # Pipeline state markers (auto-generated)
└── scripts/                              # Internal helper scripts
    ├── config.sh                         # Configuration file
    ├── utils/                            # Utility scripts
    │   └── setup_directories.sh          # Directory structure creation
    ├── pipeline/                         # Analysis pipeline scripts (Bash)
    │   ├── 01_trimmomatic.sh             # Step 1: Quality trimming
    │   ├── 02_bwa_alignment.sh           # Step 2: Read alignment
    │   ├── 03_bcftools_snp.sh            # Step 3: SNP calling
    │   ├── 04_bcftools_indel.sh          # Step 4: INDEL calling
    │   ├── 05_consensus_translation_msa.sh  # Step 7: Consensus + Translation + MSA
    │   └── 06_run_summary.sh             # Step 8: Summary generation
    └── python/                           # Python analysis scripts
        ├── generate_consensus_translate_msa.py  # MSA and translation
        └── generate_overall_summary.py   # Overall summary for all samples
```

**You only need to run two scripts:**
1. `setup.sh` - One-time project setup
2. `GC_pipeline.sh` - Analysis pipeline (can be run multiple times)

The numbered scripts in `scripts/` are internal helpers called automatically.

## Prerequisites

- BWA (Burrows-Wheeler Aligner)
- SAMtools
- BCFtools
- Trimmomatic
- FastQC
- MAFFT (for multiple sequence alignment)
- Java (for Trimmomatic)
- Conda (for Python environment management)
- Python 3.9+ with BioPython, pandas, matplotlib (for summary generation and MSA)

### Setting up Conda Environment

Create a conda environment with required Python packages:

```bash
conda create -n vac_env python=3.9
conda activate vac_env
conda install -c conda-forge biopython pandas matplotlib pillow
```

Or use a different environment name by setting the `GC_CONDA_ENV` variable:

```bash
export GC_CONDA_ENV=my_custom_env
```

## Configuration

Edit `scripts/config.sh` to set the correct paths for your system:

```bash
# Tool Paths
# If tools are in your PATH, you can use just the command name
# Otherwise, specify the full path to each tool
export BWA_PATH="bwa" # Or specify full path: /usr/local/bin/bwa
export SAMTOOLS_PATH="samtools" # Or specify full path: /usr/local/bin/samtools
export BCFTOOLS_PATH="bcftools" # Or specify full path: /usr/local/bin/bcftools
export TRIMMOMATIC_JAR="/path/to/trimmomatic.jar"
export TRIMMOMATIC_ADAPTERS="/path/to/trimmomatic/adapters/TruSeq3-PE-2.fa"
export FASTQC_PATH="fastqc" # Or specify full path: /usr/local/bin/fastqc
export MAFFT_PATH="mafft"  # MAFFT for MSA (assumes in PATH)

# Global Settings
export THREADS=8  # Adjust based on your system

# Segment Configuration for Translation and MSA
# Comma-separated list of segments to analyze (empty = analyze all)
# Common segments: PB2,PB1,PA,HA,NP,NA,M,NS
export SEGMENTS_TO_ANALYZE="HA,NA"  # Customize as needed
```

## Usage

### Option A: Fully Automated Setup (Recommended)

One command does everything - copies files, indexes genome, generates sample list:

```bash
bash setup.sh \
  --region avianflu \
  --year 2025 \
  --strain A_H5N1 \
  --raw-data /path/to/raw_data \
  --reference /path/to/2025_A_H5N1_8segments_spikein_ATGC.fasta
```

**This will automatically:**
1. ✓ Validate all required tools
2. ✓ Create project directory structure with strain-specific subdirectories
3. ✓ Copy all FASTQ files (expects `*_1.fastq.gz` and `*_2.fastq.gz` format)
4. ✓ Generate sample list automatically from FASTQ prefixes
5. ✓ Copy reference genome
6. ✓ Index reference genome with BWA

**Then run the pipeline:**
```bash
bash GC_pipeline.sh \
  --region avianflu \
  --year 2025 \
  --strain A_H5N1 \
  --samples /path/to/project/GreenCross_flu_avianflu_2025/A_H5N1/samples.txt
```

---

### Option B: Manual Setup

**Step 1: Create directory structure only:**
```bash
bash setup.sh --region southern --year 2024 --strain A_H3N2
```

**Step 2: Manually prepare data:**
```bash
# Copy FASTQ files
cp /path/to/sample*.fastq.gz /path/to/project/GreenCross_flu_southern_2024/A_H3N2/0_fastq/

# Copy reference genome
cp /path/to/2024_A_H3N2_8segments_spikein_ATGC.fasta /path/to/project/GreenCross_flu_southern_2024/A_H3N2/reference/

# Index reference genome
bwa index /path/to/project/GreenCross_flu_southern_2024/A_H3N2/reference/2024_A_H3N2_8segments_spikein_ATGC.fasta

# Create sample list
ls /path/to/project/GreenCross_flu_southern_2024/A_H3N2/0_fastq/*_1.fastq.gz | \
  xargs -n1 basename | sed 's/_1\.fastq\.gz$//' > /path/to/project/GreenCross_flu_southern_2024/A_H3N2/samples.txt
```

**Step 3: Run the pipeline:**
```bash
bash GC_pipeline.sh \
  --region southern \
  --year 2024 \
  --strain A_H3N2 \
  --samples /path/to/project/GreenCross_flu_southern_2024/A_H3N2/samples.txt
```

### Resume Interrupted Pipeline

If the pipeline is interrupted, simply re-run the same command:

```bash
bash GC_pipeline.sh \
  --region southern \
  --year 2024 \
  --strain A_H3N2 \
  --samples /path/to/project/GreenCross_flu_southern_2024/A_H3N2/samples.txt
```

Already completed samples will be automatically skipped:
```
[14:23:45] ⏭️  Sample sample001 (1/3) - Sample already processed
[14:23:46] ⏭️  Sample sample002 (2/3) - Sample already processed
========================================================================
[PIPELINE] 2024-01-02 14:23:47 - Processing Sample: sample003 (3/3)
========================================================================
[14:23:47] ▶️  [1/8] Initial FastQC - Sample: sample003
...
```

## Pipeline Steps

For each sample, the pipeline performs:

1. **Initial FastQC** - Quality assessment of raw reads
2. **Trimmomatic** - Adapter removal and quality trimming
3. **Post-Trimming FastQC** - Quality assessment after trimming
4. **BWA Alignment** - Map reads to reference genome
5. **BCFtools SNP Calling** - Identify SNP variants
6. **BCFtools INDEL Calling** - Identify insertion/deletion variants
7. **Consensus Generation, Translation & MSA** - Apply variants to reference, translate to protein, perform multiple sequence alignment
8. **Summary Report** - Generate analysis summary

## Detailed Pipeline Description

### Step 1: Initial FastQC
**Purpose**: Assess raw sequencing data quality before processing

**Input**:
- Raw paired-end FASTQ files (`{sample}_1.fastq.gz`, `{sample}_2.fastq.gz`)

**Output**:
- HTML quality reports in `0_fastq/fastqc_pre_trim/`
- Metrics: per-base quality scores, GC content, adapter contamination

### Step 2: Trimmomatic
**Purpose**: Remove adapters and low-quality bases

**Parameters**:
- `ILLUMINACLIP`: Adapter trimming (seed mismatches:2, palindrome threshold:30, simple clip:10)
- `SLIDINGWINDOW:4:15`: Trim when average quality in 4-base window drops below 15
- `LEADING:10`: Cut bases from start if quality < 10
- `TRAILING:10`: Cut bases from end if quality < 10
- `MINLEN:50`: Discard reads shorter than 50bp

**Input**: Raw FASTQ files

**Output**:
- Trimmed paired FASTQ files in `1_trimming/{sample}/`
- Trimming log with survival statistics

### Step 3: Post-Trimming FastQC
**Purpose**: Verify quality improvement after trimming

**Input**: Trimmed FASTQ files

**Output**:
- HTML quality reports in `1_trimming/fastqc_post_trim/`
- Confirms adapter removal and quality improvement

### Step 4: BWA Alignment
**Purpose**: Map reads to influenza reference genome

**Algorithm**: BWA-MEM (optimized for 70bp-1Mbp reads)

**Process**:
1. Align trimmed reads to reference
2. Convert SAM to sorted BAM
3. Index BAM file
4. Generate alignment statistics (flagstat, idxstats)

**Input**: Trimmed FASTQ files + indexed reference genome

**Output**:
- `{sample}.aln.sorted.bam` - Sorted alignment file
- `{sample}.aln.sorted.bam.bai` - BAM index
- `{sample}.aln.flagstat` - Alignment statistics (mapped reads, properly paired)
- `{sample}.aln.idxstat` - Per-chromosome alignment counts

### Step 5: BCFtools SNP Calling
**Purpose**: Identify single nucleotide polymorphisms

**Process**:
1. **Mpileup**: Generate pileup of bases at each position
   - Minimum mapping quality: 30 (`-q30`)
   - Max depth: 1,000,000 (`--max-depth 1000000`)
2. **Call**: Identify SNP variants
   - Consensus variant calling (`-cv`)
   - Exclude INDELs (`-V indels`)
3. **Filter**: Require minimum depth of 10,000 (`DP>=10000`)

**Input**: Sorted BAM file + reference genome

**Output**:
- `{sample}.aln.snp.vcf` - Filtered SNP variants

### Step 6: BCFtools INDEL Calling
**Purpose**: Identify insertions and deletions

**Process**:
1. **Mpileup**: Generate pileup (same as SNP calling)
2. **Call**: Identify INDEL variants
   - Exclude SNPs (`-V snps`)
3. **No filtering applied** (keeps all INDELs)

**Input**: Sorted BAM file + reference genome

**Output**:
- `{sample}.aln.indel.vcf` - All INDEL variants

### Step 7: Consensus Generation, Translation & MSA
**Purpose**: Generate consensus sequences, translate to protein, and align with reference

**Process**:
1. **Consensus Generation**:
   - Merge SNP and INDEL VCF files
   - Apply variants to reference genome using `bcftools consensus`
   - Output: `{sample}_consensus.fasta`

2. **Segment Extraction**:
   - Extract configured segments (default: HA, NA)
   - Parse multi-FASTA reference and consensus

3. **ORF Detection & Translation**:
   - Search for first ATG (start codon) in each segment
   - Translate DNA to protein using BioPython
   - Handle reading frames and stop codons
   - Output: `{sample}_{segment}_proteins.fasta`

4. **Multiple Sequence Alignment**:
   - Align reference and sample proteins using MAFFT
   - Generate color-coded visualization (PNG)
   - Identify amino acid differences
   - Output:
     - `{sample}_{segment}_aligned.fasta` - Aligned sequences
     - `{sample}_{segment}_alignment.png` - Visual alignment
     - `{sample}_{segment}_differences.csv` - Amino acid changes (position, ref, sample)

**Input**:
- Reference genome
- SNP VCF file
- INDEL VCF file

**Output** (per segment):
- Consensus FASTA
- Protein FASTA (reference + sample)
- MAFFT alignment FASTA
- Alignment visualization PNG
- Differences table CSV

### Step 8: Summary Report
**Purpose**: Aggregate statistics for all processed samples

**Process**:
1. Parse Trimmomatic logs for read counts
2. Parse flagstat files for alignment metrics
3. Count variants from VCF files
4. Generate unified summary CSV with statistics
5. Create read count visualization

**Input**:
- All sample logs and output files
- Sample list file

**Output**:
- `overall_summary.csv` - All samples + average statistics
  - Columns: SampleID, Total_Reads, After_Trimmomatic (%), Properly_Paired (%), Mapped_Reads (%), SNPs, INDELs
  - Final row: AVERAGE with mean values and totals
- `overall_read_counts.png` - Stacked bar chart showing read counts across all samples

## Output Structure

After running the pipeline, your project directory will contain:

```
/path/to/project/GreenCross_flu_southern_2024/
├── A_H3N2/                    # Strain-specific directory
│   ├── 0_fastq/               # Raw FASTQ files
│   │   └── fastqc_pre_trim/  # Pre-trimming QC reports
│   ├── 1_trimming/            # Trimmed reads
│   │   ├── fastqc_post_trim/ # Post-trimming QC reports
│   │   └── {sampleID}/       # Sample-specific trimmed files
│   ├── 2_align/               # Alignment results
│   │   ├── {sampleID}/       # Sample BAM files and indices
│   │   ├── flagstat/         # Alignment statistics
│   │   └── idxstat/          # Index statistics
│   ├── 3_variantCall/         # Variant calling results
│   │   └── {sampleID}/       # SNP and INDEL VCF files
│   ├── 4_msa/                 # MSA results (Step 7)
│   │   └── {sampleID}/       # Per-sample MSA results
│   │       ├── {sampleID}_consensus.fasta
│   │       ├── {sampleID}_{segment}_proteins.fasta
│   │       ├── {sampleID}_{segment}_aligned.fasta
│   │       ├── {sampleID}_{segment}_alignment.png
│   │       └── {sampleID}_{segment}_differences.csv
│   ├── 5_summary/             # Summary reports (Step 8)
│   │   ├── overall_summary.csv
│   │   └── overall_read_counts.png
│   ├── logs/                  # Pipeline execution logs
│   │   └── {sampleID}/       # Per-sample step logs
│   │       ├── 01_fastqc_pre.log
│   │       ├── 02_trimmomatic.log
│   │       ├── 03_fastqc_post.log
│   │       ├── 04_bwa_alignment.log
│   │       ├── 05_bcftools_snp.log
│   │       ├── 06_bcftools_indel.log
│   │       ├── 07_consensus_translation_msa.log
│   │       └── 08_summary.log
│   ├── reference/             # Reference genome files
│   └── samples.txt            # Sample list
├── A_H1N1/                    # Another strain (if processing multiple strains)
│   └── (same structure)
└── B_Victoria/                # Another strain (if processing multiple strains)
    └── (same structure)
```

## Summary Report

At the end of execution, you'll see a comprehensive summary:

```
================ Pipeline Execution Summary ================
Total samples in list: 10
Successfully processed: 8
Skipped (already done): 2
Failed: 0
⏱️  Total execution time: 45m 32s
============================================================
```

## Troubleshooting

### Pipeline Fails for a Sample

1. Check the logs in `logs/{sampleID}/` directory
2. Review the specific step that failed
3. Fix the issue (e.g., missing file, insufficient disk space)
4. Remove the marker file: `rm .pipeline/sample_{sampleID}_done`
5. Re-run the pipeline

### Reset Entire Pipeline

To start fresh and re-process all samples:

```bash
rm -rf .pipeline/
```

### Disable Verbose Output

```bash
VERBOSE=0 bash GC_pipeline.sh \
  --region southern \
  --year 2024 \
  --strain A_H3N2 \
  --samples /path/to/project/GreenCross_flu_southern_2024/A_H3N2/samples.txt
```

## Advanced Usage

### Process Specific Sample Subset

Create a subset sample list and run:

```bash
# Create subset in the strain directory
echo "sample005" > /path/to/project/GreenCross_flu_southern_2024/A_H3N2/subset.txt
echo "sample006" >> /path/to/project/GreenCross_flu_southern_2024/A_H3N2/subset.txt

# Run pipeline
bash GC_pipeline.sh \
  --region southern \
  --year 2024 \
  --strain A_H3N2 \
  --samples /path/to/project/GreenCross_flu_southern_2024/A_H3N2/subset.txt
```

### Custom Thread Count

Edit `scripts/config.sh` and adjust:
```bash
export THREADS=16  # Use 16 threads
```

## Reference Genome Naming Convention

Reference genomes must follow this naming pattern:
```
YYYY_STRAIN_8segments_spikein_ATGC.fasta
```

Examples:
- `2024_A_H3N2_8segments_spikein_ATGC.fasta`
- `2024_B_Victoria_8segments_spikein_ATGC.fasta`

## Contributing

When modifying scripts, ensure:
- All bash scripts use `#!/usr/bin/env bash` and `set -euo pipefail`
- Error messages are sent to stderr (`>&2`)
- Proper exit codes are returned
- New features are documented in this README

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this pipeline in your research, please cite:

```
Influenza Vaccine NGS Analysis Pipeline
https://github.com/YOUR_USERNAME/vaccine-ngs-pipeline
```
