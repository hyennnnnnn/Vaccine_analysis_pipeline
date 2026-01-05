#!/usr/bin/env python3
"""
Consensus Generation, Translation, and Multiple Sequence Alignment
This script:
1. Generates consensus sequence from VCF files (applying SNPs and INDELs)
2. Extracts specified segments from reference and consensus
3. Translates DNA to protein (auto-detecting ORF from first ATG)
4. Performs MSA using MAFFT
5. Generates visualizations and summary statistics
"""

import sys
import os
import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from collections import defaultdict

def run_command(cmd, description=""):
    """Run a shell command and handle errors."""
    print(f"Running: {description if description else cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {cmd}")
        print(f"STDERR: {result.stderr}")
        raise RuntimeError(f"Command failed: {cmd}")
    return result.stdout

def generate_consensus(reference_file, snp_vcf, indel_vcf, output_file, bcftools_path, samtools_path):
    """Generate consensus sequence by applying SNPs and INDELs to reference."""
    print(f"\n=== Generating Consensus Sequence ===")

    # Merge SNP and INDEL VCF files
    merged_vcf = output_file.replace(".fasta", "_merged.vcf.gz")

    # Check if VCF files exist and have variants
    snp_has_variants = os.path.exists(snp_vcf) and os.path.getsize(snp_vcf) > 0
    indel_has_variants = os.path.exists(indel_vcf) and os.path.getsize(indel_vcf) > 0

    if not snp_has_variants and not indel_has_variants:
        print("Warning: No variants found in either SNP or INDEL VCF. Using reference as consensus.")
        # Copy reference as consensus
        run_command(f"cp {reference_file} {output_file}", "Copying reference as consensus")
        return output_file

    # Compress and index VCF files if they have variants
    vcf_files = []
    if snp_has_variants:
        snp_vcf_gz = snp_vcf + ".gz"
        run_command(f"{bcftools_path} view {snp_vcf} -Oz -o {snp_vcf_gz}", "Compressing SNP VCF")
        run_command(f"{bcftools_path} index {snp_vcf_gz}", "Indexing SNP VCF")
        vcf_files.append(snp_vcf_gz)

    if indel_has_variants:
        indel_vcf_gz = indel_vcf + ".gz"
        run_command(f"{bcftools_path} view {indel_vcf} -Oz -o {indel_vcf_gz}", "Compressing INDEL VCF")
        run_command(f"{bcftools_path} index {indel_vcf_gz}", "Indexing INDEL VCF")
        vcf_files.append(indel_vcf_gz)

    # Merge VCF files if we have both
    if len(vcf_files) == 2:
        run_command(f"{bcftools_path} concat -a {' '.join(vcf_files)} -Oz -o {merged_vcf}",
                   "Merging SNP and INDEL VCFs")
        run_command(f"{bcftools_path} index {merged_vcf}", "Indexing merged VCF")
        vcf_to_use = merged_vcf
    else:
        vcf_to_use = vcf_files[0]

    # Index reference if not already indexed
    if not os.path.exists(reference_file + ".fai"):
        run_command(f"{samtools_path} faidx {reference_file}", "Indexing reference")

    # Generate consensus
    run_command(f"{bcftools_path} consensus -f {reference_file} {vcf_to_use} > {output_file}",
               "Generating consensus sequence")

    print(f"✓ Consensus sequence generated: {output_file}")
    return output_file

def extract_segments(fasta_file, segment_names):
    """Extract specified segments from FASTA file.
    Assumes segment names are in the FASTA headers.
    Returns dict: {segment_name: SeqRecord}
    """
    segments = {}
    records = list(SeqIO.parse(fasta_file, "fasta"))

    if not segment_names:  # If empty, analyze all segments
        for record in records:
            segments[record.id] = record
        return segments

    # Search for specified segments
    for seg_name in segment_names:
        found = False
        for record in records:
            # Case-insensitive search in the header
            if seg_name.upper() in record.id.upper() or seg_name.upper() in record.description.upper():
                segments[seg_name] = record
                found = True
                break

        if not found:
            print(f"Warning: Segment '{seg_name}' not found in {fasta_file}")

    return segments

def find_first_orf(seq_record):
    """Find the first ORF starting with ATG.
    Returns (start_position, translated_protein) or (None, None) if not found.
    """
    sequence = str(seq_record.seq).upper()

    # Search for ATG in all three frames
    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon == "ATG":
                # Found start codon, translate from here
                orf_seq = sequence[i:]
                # Trim to multiple of 3
                orf_seq = orf_seq[:len(orf_seq) - (len(orf_seq) % 3)]

                if len(orf_seq) >= 3:
                    protein = Seq(orf_seq).translate(to_stop=False)
                    return i, protein

    # If no ATG found, translate entire sequence from position 0
    print(f"  Warning: No ATG found in {seq_record.id}, translating from position 0")
    orf_seq = sequence[:len(sequence) - (len(sequence) % 3)]
    if len(orf_seq) >= 3:
        protein = Seq(orf_seq).translate(to_stop=False)
        return 0, protein

    return None, None

def translate_segments(segments_dict, prefix):
    """Translate DNA segments to protein.
    Returns dict: {segment_name: SeqRecord of protein}
    """
    proteins = {}

    for seg_name, seq_record in segments_dict.items():
        print(f"  Translating {seg_name}...")
        start_pos, protein_seq = find_first_orf(seq_record)

        if protein_seq is None:
            print(f"    Warning: Could not translate {seg_name}")
            continue

        protein_record = SeqRecord(
            protein_seq,
            id=f"{prefix}_{seg_name}",
            description=f"{prefix} {seg_name} translation (ORF start: {start_pos})"
        )
        proteins[seg_name] = protein_record
        print(f"    ✓ Translated {len(protein_seq)} amino acids (ORF start: {start_pos})")

    return proteins

def run_mafft(input_fasta, output_fasta, mafft_path="mafft"):
    """Run MAFFT for multiple sequence alignment."""
    print(f"\n=== Running MAFFT on {input_fasta} ===")
    cmd = f"{mafft_path} --auto --thread -1 {input_fasta} > {output_fasta}"
    run_command(cmd, "MAFFT alignment")
    print(f"✓ Alignment saved: {output_fasta}")
    return output_fasta

def visualize_alignment(alignment_file, output_image, max_seq_display=50):
    """Visualize protein alignment with color-coded amino acids."""
    print(f"\n=== Generating alignment visualization ===")

    alignment = AlignIO.read(alignment_file, "fasta")

    # Amino acid color scheme (similar to ClustalX)
    aa_colors = {
        'G': '#FF9900', 'P': '#FF9900',  # Orange (small/proline)
        'S': '#33CC00', 'T': '#33CC00', 'N': '#33CC00', 'Q': '#33CC00',  # Green (hydroxyl/amine)
        'D': '#FF0000', 'E': '#FF0000',  # Red (acidic)
        'K': '#0000FF', 'R': '#0000FF', 'H': '#6666FF',  # Blue (basic)
        'F': '#00CCCC', 'Y': '#00CCCC', 'W': '#00CCCC',  # Cyan (aromatic)
        'I': '#33CC00', 'L': '#33CC00', 'M': '#33CC00', 'V': '#33CC00',  # Green (aliphatic)
        'A': '#33CC00', 'C': '#FFFF00',  # Yellow (cysteine), green (alanine)
        '-': '#FFFFFF', 'X': '#999999'  # Gap (white), unknown (gray)
    }

    num_seqs = len(alignment)
    align_length = alignment.get_alignment_length()

    # Limit display length if too long
    display_length = min(align_length, max_seq_display)

    fig, ax = plt.subplots(figsize=(max(20, display_length * 0.3), num_seqs * 0.5))

    # Draw alignment
    for i, record in enumerate(alignment):
        for j, aa in enumerate(str(record.seq)[:display_length]):
            color = aa_colors.get(aa.upper(), '#CCCCCC')
            rect = mpatches.Rectangle((j, num_seqs - i - 1), 1, 1,
                                     facecolor=color, edgecolor='black', linewidth=0.5)
            ax.add_patch(rect)

            # Add amino acid letter
            ax.text(j + 0.5, num_seqs - i - 0.5, aa,
                   ha='center', va='center', fontsize=8, fontweight='bold')

    # Set axis properties
    ax.set_xlim(0, display_length)
    ax.set_ylim(0, num_seqs)
    ax.set_yticks([i + 0.5 for i in range(num_seqs)])
    ax.set_yticklabels([record.id for record in reversed(alignment)])
    ax.set_xlabel('Position', fontsize=12, fontweight='bold')
    ax.set_title(f'Protein Sequence Alignment (showing first {display_length} positions)',
                fontsize=14, fontweight='bold')
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Visualization saved: {output_image}")

def analyze_alignment_differences(alignment_file, output_csv):
    """Analyze differences between reference and consensus in alignment."""
    print(f"\n=== Analyzing alignment differences ===")

    alignment = AlignIO.read(alignment_file, "fasta")

    if len(alignment) < 2:
        print("Warning: Less than 2 sequences in alignment, cannot compare")
        return

    differences = []
    ref_seq = None
    cons_seq = None

    # Identify reference and consensus sequences
    for record in alignment:
        if 'reference' in record.id.lower() or 'ref' in record.id.lower():
            ref_seq = str(record.seq)
            ref_id = record.id
        elif 'consensus' in record.id.lower() or 'cons' in record.id.lower():
            cons_seq = str(record.seq)
            cons_id = record.id

    if ref_seq is None or cons_seq is None:
        # If not found by name, assume first is reference, second is consensus
        ref_seq = str(alignment[0].seq)
        ref_id = alignment[0].id
        cons_seq = str(alignment[1].seq)
        cons_id = alignment[1].id

    # Find differences
    for pos in range(len(ref_seq)):
        ref_aa = ref_seq[pos]
        cons_aa = cons_seq[pos]

        if ref_aa != cons_aa:
            differences.append({
                'Position': pos + 1,
                'Reference_AA': ref_aa,
                'Consensus_AA': cons_aa,
                'Change': f"{ref_aa}{pos + 1}{cons_aa}",
                'Type': 'Substitution' if ref_aa != '-' and cons_aa != '-' else 'Indel'
            })

    # Save to CSV
    if differences:
        df = pd.DataFrame(differences)
        df.to_csv(output_csv, index=False)
        print(f"✓ Found {len(differences)} differences")
        print(f"✓ Differences saved: {output_csv}")

        # Print summary
        print(f"\n  Summary:")
        print(f"    Total differences: {len(differences)}")
        print(f"    Substitutions: {sum(1 for d in differences if d['Type'] == 'Substitution')}")
        print(f"    Indels: {sum(1 for d in differences if d['Type'] == 'Indel')}")
    else:
        print("  ✓ No differences found (sequences are identical)")
        # Create empty CSV
        pd.DataFrame(columns=['Position', 'Reference_AA', 'Consensus_AA', 'Change', 'Type']).to_csv(output_csv, index=False)

def main():
    if len(sys.argv) != 10:
        print("Usage: python3 generate_consensus_translate_msa.py \\")
        print("  <sample_id> <reference_fasta> <snp_vcf> <indel_vcf> \\")
        print("  <bcftools_path> <samtools_path> <mafft_path> \\")
        print("  <segments_to_analyze> <output_dir>")
        print("\nArguments:")
        print("  segments_to_analyze: comma-separated (e.g., 'HA,NA') or empty for all")
        sys.exit(1)

    sample_id = sys.argv[1]
    reference_fasta = sys.argv[2]
    snp_vcf = sys.argv[3]
    indel_vcf = sys.argv[4]
    bcftools_path = sys.argv[5]
    samtools_path = sys.argv[6]
    mafft_path = sys.argv[7]
    segments_str = sys.argv[8]
    output_dir = sys.argv[9]

    # Parse segments
    if segments_str and segments_str.strip():
        segments_to_analyze = [s.strip() for s in segments_str.split(',')]
    else:
        segments_to_analyze = []  # Empty = analyze all

    print("="*70)
    print(f"Consensus Generation, Translation, and MSA for {sample_id}")
    print("="*70)
    print(f"Segments to analyze: {segments_to_analyze if segments_to_analyze else 'ALL'}")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Generate consensus sequence
    consensus_fasta = os.path.join(output_dir, f"{sample_id}_consensus.fasta")
    generate_consensus(reference_fasta, snp_vcf, indel_vcf, consensus_fasta,
                      bcftools_path, samtools_path)

    # Step 2: Extract segments from reference and consensus
    print("\n=== Extracting Segments ===")
    print("From reference...")
    ref_segments = extract_segments(reference_fasta, segments_to_analyze)
    print(f"  ✓ Extracted {len(ref_segments)} segments from reference")

    print("From consensus...")
    cons_segments = extract_segments(consensus_fasta, segments_to_analyze)
    print(f"  ✓ Extracted {len(cons_segments)} segments from consensus")

    # Step 3: Translate segments
    print("\n=== Translating Segments to Protein ===")
    print("Reference segments:")
    ref_proteins = translate_segments(ref_segments, "Reference")

    print("\nConsensus segments:")
    cons_proteins = translate_segments(cons_segments, "Consensus")

    # Step 4: Run MAFFT for each segment
    print("\n=== Running Multiple Sequence Alignment (MAFFT) ===")

    for seg_name in ref_proteins.keys():
        if seg_name not in cons_proteins:
            print(f"Warning: Segment {seg_name} not in consensus, skipping")
            continue

        print(f"\nProcessing segment: {seg_name}")

        # Create input file with both reference and consensus
        input_fasta = os.path.join(output_dir, f"{sample_id}_{seg_name}_proteins.fasta")
        with open(input_fasta, 'w') as f:
            SeqIO.write([ref_proteins[seg_name], cons_proteins[seg_name]], f, "fasta")

        # Run MAFFT
        aligned_fasta = os.path.join(output_dir, f"{sample_id}_{seg_name}_aligned.fasta")
        run_mafft(input_fasta, aligned_fasta, mafft_path)

        # Visualize alignment
        alignment_image = os.path.join(output_dir, f"{sample_id}_{seg_name}_alignment.png")
        visualize_alignment(aligned_fasta, alignment_image)

        # Analyze differences
        diff_csv = os.path.join(output_dir, f"{sample_id}_{seg_name}_differences.csv")
        analyze_alignment_differences(aligned_fasta, diff_csv)

    print("\n" + "="*70)
    print(f"✓ Analysis complete for {sample_id}")
    print(f"✓ Output directory: {output_dir}")
    print("="*70)

if __name__ == "__main__":
    main()
