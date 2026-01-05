#!/usr/bin/env python3
"""
Overall Summary Generator
Collects summary data from all samples and generates:
1. Comprehensive summary table (CSV)
2. Visualization plots (PNG)
"""

import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

def parse_trimmomatic_log(log_file):
    """Parses a Trimmomatic log file to extract read counts."""
    if not os.path.exists(log_file):
        return 0, 0

    with open(log_file, 'r') as f:
        content = f.read()

    input_reads_match = re.search(r'Input Read Pairs: (\d+)', content)
    both_surviving_match = re.search(r'Both Surviving: (\d+)', content)

    input_pairs = int(input_reads_match.group(1)) if input_reads_match else 0
    trimmed_pairs = int(both_surviving_match.group(1)) if both_surviving_match else 0

    # Convert pairs to individual reads
    return input_pairs * 2, trimmed_pairs * 2

def parse_flagstat_file(flagstat_file):
    """Parses a Samtools flagstat file to extract mapped and properly paired read counts."""
    if not os.path.exists(flagstat_file):
        return 0, 0, 0

    with open(flagstat_file, 'r') as f:
        content = f.read()

    # Total reads
    total_reads_match = re.search(r'(\d+) \+ \d+ in total', content)
    # Mapped reads
    mapped_reads_match = re.search(r'(\d+) \+ \d+ mapped', content)
    # Properly paired reads
    properly_paired_match = re.search(r'(\d+) \+ \d+ properly paired', content)

    total_reads = int(total_reads_match.group(1)) if total_reads_match else 0
    mapped_reads = int(mapped_reads_match.group(1)) if mapped_reads_match else 0
    properly_paired = int(properly_paired_match.group(1)) if properly_paired_match else 0

    return total_reads, mapped_reads, properly_paired

def collect_sample_data(sample_id, strain_dir):
    """Collect all relevant statistics for a single sample."""

    # Define file paths
    trimmomatic_log = os.path.join(strain_dir, f"1_trimming/{sample_id}/{sample_id}_trim.log")
    flagstat_file = os.path.join(strain_dir, f"2_align/{sample_id}/{sample_id}.aln.flagstat")
    snp_vcf = os.path.join(strain_dir, f"3_variantCall/{sample_id}/{sample_id}.aln.snp.vcf")
    indel_vcf = os.path.join(strain_dir, f"3_variantCall/{sample_id}/{sample_id}.aln.indel.vcf")

    # Get read counts
    input_reads, trimmed_reads = parse_trimmomatic_log(trimmomatic_log)
    total_reads, mapped_reads, properly_paired = parse_flagstat_file(flagstat_file)

    # Count variants
    snp_count = count_vcf_variants(snp_vcf)
    indel_count = count_vcf_variants(indel_vcf)

    # Calculate percentages
    trim_pct = (trimmed_reads / input_reads * 100) if input_reads > 0 else 0
    mapped_pct = (mapped_reads / trimmed_reads * 100) if trimmed_reads > 0 else 0
    paired_pct = (properly_paired / trimmed_reads * 100) if trimmed_reads > 0 else 0

    return {
        'SampleID': sample_id,
        'Total_Reads': input_reads,
        'After_Trimmomatic': trimmed_reads,
        'Trimmomatic_Pct': trim_pct,
        'Mapped_Reads': mapped_reads,
        'Mapped_Pct': mapped_pct,
        'Properly_Paired': properly_paired,
        'Properly_Paired_Pct': paired_pct,
        'SNP_Count': snp_count,
        'INDEL_Count': indel_count
    }

def count_vcf_variants(vcf_file):
    """Counts non-header variants in a VCF file."""
    if not os.path.exists(vcf_file):
        return 0

    count = 0
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                count += 1
    return count

def format_number_with_pct(value, pct):
    """Format number with percentage like: 32,149,428 (98.04%)"""
    return f"{value:,} ({pct:.2f}%)"

def generate_overall_summary(sample_list_file, strain_dir, output_dir):
    """Generate overall summary for all samples."""

    # Read sample list
    with open(sample_list_file, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]

    print(f"Collecting data for {len(samples)} samples...")

    # Collect data for all samples
    all_data = []
    for sample_id in samples:
        try:
            data = collect_sample_data(sample_id, strain_dir)
            all_data.append(data)
            print(f"  ✓ {sample_id}")
        except Exception as e:
            print(f"  ✗ {sample_id}: {e}")

    if not all_data:
        print("No data collected. Exiting.")
        return

    # Create DataFrame
    df = pd.DataFrame(all_data)

    # --- 1. Generate Unified Summary CSV (formatted data + statistics) ---
    # Sample data with formatting
    formatted_df = pd.DataFrame({
        'SampleID': df['SampleID'],
        'Total_Reads': df.apply(lambda x: f"{int(x['Total_Reads']):,}", axis=1),
        'After_Trimmomatic': df.apply(lambda x: format_number_with_pct(int(x['After_Trimmomatic']), x['Trimmomatic_Pct']), axis=1),
        'Properly_Paired': df.apply(lambda x: format_number_with_pct(int(x['Properly_Paired']), x['Properly_Paired_Pct']), axis=1),
        'Mapped_Reads': df.apply(lambda x: format_number_with_pct(int(x['Mapped_Reads']), x['Mapped_Pct']), axis=1),
        'SNPs': df['SNP_Count'],
        'INDELs': df['INDEL_Count']
    })

    # Add summary statistics row
    summary_row = pd.DataFrame([{
        'SampleID': f'AVERAGE (n={len(df)})',
        'Total_Reads': f"{df['Total_Reads'].mean():,.0f}",
        'After_Trimmomatic': format_number_with_pct(df['After_Trimmomatic'].mean(), df['Trimmomatic_Pct'].mean()),
        'Properly_Paired': format_number_with_pct(df['Properly_Paired'].mean(), df['Properly_Paired_Pct'].mean()),
        'Mapped_Reads': format_number_with_pct(df['Mapped_Reads'].mean(), df['Mapped_Pct'].mean()),
        'SNPs': f"{df['SNP_Count'].mean():.1f} (Total: {df['SNP_Count'].sum():,})",
        'INDELs': f"{df['INDEL_Count'].mean():.1f} (Total: {df['INDEL_Count'].sum():,})"
    }])

    # Combine sample data and summary
    combined_df = pd.concat([formatted_df, summary_row], ignore_index=True)

    summary_csv_path = os.path.join(output_dir, "overall_summary.csv")
    combined_df.to_csv(summary_csv_path, index=False)
    print(f"\n✓ Overall summary saved to: {summary_csv_path}")

    # --- 2. Generate Read Count Visualization ---
    fig, ax = plt.subplots(figsize=(14, 8))

    x_pos = range(len(df))
    width = 0.6

    # Create stacked bars
    p1 = ax.bar(x_pos, df['Total_Reads'], width, label='Total Reads', color='#3498db')
    p2 = ax.bar(x_pos, df['After_Trimmomatic'], width, label='After Trimmomatic', color='#2ecc71')
    p3 = ax.bar(x_pos, df['Properly_Paired'], width, label='Properly Paired', color='#f39c12')

    ax.set_xlabel('Sample ID', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Reads', fontsize=12, fontweight='bold')
    ax.set_title('Read Count Summary Across All Samples', fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df['SampleID'], rotation=45, ha='right')
    ax.legend(loc='upper right')
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    read_plot_path = os.path.join(output_dir, "overall_read_counts.png")
    plt.savefig(read_plot_path, dpi=300)
    plt.close()
    print(f"✓ Read count plot saved to: {read_plot_path}")

    # --- 3. Print Summary Statistics to Console ---
    summary_stats = {
        'Total_Samples': len(df),
        'Avg_Total_Reads': df['Total_Reads'].mean(),
        'Avg_After_Trimmomatic': df['After_Trimmomatic'].mean(),
        'Avg_Trimmomatic_Retention_Pct': df['Trimmomatic_Pct'].mean(),
        'Avg_Properly_Paired': df['Properly_Paired'].mean(),
        'Avg_Properly_Paired_Pct': df['Properly_Paired_Pct'].mean(),
        'Avg_Mapped_Reads': df['Mapped_Reads'].mean(),
        'Avg_Mapped_Pct': df['Mapped_Pct'].mean(),
        'Total_SNPs': df['SNP_Count'].sum(),
        'Total_INDELs': df['INDEL_Count'].sum(),
        'Avg_SNPs_Per_Sample': df['SNP_Count'].mean(),
        'Avg_INDELs_Per_Sample': df['INDEL_Count'].mean()
    }

    print("\n" + "="*70)
    print("OVERALL SUMMARY STATISTICS")
    print("="*70)
    print(f"Total Samples: {summary_stats['Total_Samples']}")
    print(f"Average Total Reads: {summary_stats['Avg_Total_Reads']:,.0f}")
    print(f"Average After Trimmomatic: {summary_stats['Avg_After_Trimmomatic']:,.0f} ({summary_stats['Avg_Trimmomatic_Retention_Pct']:.2f}%)")
    print(f"Average Properly Paired: {summary_stats['Avg_Properly_Paired']:,.0f} ({summary_stats['Avg_Properly_Paired_Pct']:.2f}%)")
    print(f"Average Mapped Reads: {summary_stats['Avg_Mapped_Reads']:,.0f} ({summary_stats['Avg_Mapped_Pct']:.2f}%)")
    print(f"Total SNPs Detected: {summary_stats['Total_SNPs']:,} (Avg: {summary_stats['Avg_SNPs_Per_Sample']:.1f} per sample)")
    print(f"Total INDELs Detected: {summary_stats['Total_INDELs']:,} (Avg: {summary_stats['Avg_INDELs_Per_Sample']:.1f} per sample)")
    print("="*70)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 generate_overall_summary.py <sample_list_file> <strain_dir> <output_dir>")
        print("\nExample:")
        print("  python3 generate_overall_summary.py \\")
        print("    /path/to/project/samples.txt \\")
        print("    /path/to/project/A_H5N1 \\")
        print("    /path/to/project/A_H5N1/5_summary")
        sys.exit(1)

    sample_list_file = sys.argv[1]
    strain_dir = sys.argv[2]
    output_dir = sys.argv[3]

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    generate_overall_summary(sample_list_file, strain_dir, output_dir)
