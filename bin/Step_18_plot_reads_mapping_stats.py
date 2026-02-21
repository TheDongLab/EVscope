#!/usr/bin/env python3
"""
Compute total read counts from featureCounts TSV files and plot a pie chart.
Also generates a summary TSV table of read counts per feature.
Legend order: 5'UTR, Exon, 3'UTR, Intron, Promoter, Downstream 2Kb, Intergenic, ENCODE_blacklist.
"""
import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set publication-quality parameters
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42

def sum_read_counts(tsv_file):
    """Sum counts in the last column of a TSV file, skipping headers and comments."""
    total = 0.0
    try:
        with open(tsv_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('Geneid'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    try:
                        total += float(parts[-1])
                    except (ValueError, IndexError):
                        continue # Skip lines that cannot be parsed
    except FileNotFoundError:
        print(f"Warning: File not found {tsv_file}. Counting as 0.", file=sys.stderr)
        return 0.0
    return total

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot read mapping statistics as a publication-quality pie chart.'
    )
    parser.add_argument('--input_5UTR_readcounts', required=True, help="TSV for 5'UTR counts")
    parser.add_argument('--input_exon_readcounts', required=True, help='TSV for exon counts')
    parser.add_argument('--input_3UTR_readcounts', required=True, help="TSV for 3'UTR counts")
    parser.add_argument('--input_intron_readcounts', required=True, help='TSV for intron counts')
    parser.add_argument('--input_promoters_readcounts', required=True, help='TSV for promoter counts')
    parser.add_argument('--input_downstream_2Kb_readcounts', required=True, help='TSV for downstream counts')
    parser.add_argument('--input_intergenic_readcounts', required=True, help='TSV for intergenic counts')
    parser.add_argument('--input_ENCODE_blacklist_readcounts', required=True, help='TSV for ENCODE blacklist counts')
    parser.add_argument('--sampleName', required=True, help='Sample name for the plot title')
    parser.add_argument('--output_dir', required=True, help='Directory to save output files')
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Sum reads for each genomic region
    counts = [
        sum_read_counts(args.input_5UTR_readcounts),
        sum_read_counts(args.input_exon_readcounts),
        sum_read_counts(args.input_3UTR_readcounts),
        sum_read_counts(args.input_intron_readcounts),
        sum_read_counts(args.input_promoters_readcounts),
        sum_read_counts(args.input_downstream_2Kb_readcounts),
        sum_read_counts(args.input_intergenic_readcounts),
        sum_read_counts(args.input_ENCODE_blacklist_readcounts)
    ]

    labels = ["5'UTR", 'Exon', "3'UTR", 'Intron', 'Promoter',
              'Downstream 2kb', 'Intergenic', 'ENCODE Blacklist']

    total = sum(counts)
    if total == 0:
        print('Error: Total read count is zero. Cannot generate plot.', file=sys.stderr)
        sys.exit(1)

    # ---------------------------------------------------------
    # Generate Summary Table (TSV)
    # ---------------------------------------------------------
    tsv_filename = f"{args.sampleName}_reads_mapping_readcounts_summary.tsv"
    tsv_path = os.path.join(args.output_dir, tsv_filename)
    
    try:
        with open(tsv_path, 'w') as f_out:
            # Header line: using 'Genomic Feature' instead of 'Meta-gene region type'
            f_out.write("Genomic_feature\tRead_counts\n")
            for label, count in zip(labels, counts):
                f_out.write(f"{label}\t{int(count)}\n")
        print(f"Summary table saved successfully to: {tsv_path}")
    except IOError as e:
        print(f"Error writing summary table: {e}", file=sys.stderr)

    # ---------------------------------------------------------
    # Generate Pie Chart
    # ---------------------------------------------------------
    
    # Prepare data for the pie chart
    sizes = [c / total for c in counts]
    colors = ['red', 'orange', 'cyan', 'green', 'purple', 'blue', 'pink', 'gray']

    # Create the plot
    fig, ax = plt.subplots(figsize=(7, 3.5), dpi=300)
    wedges, _ = ax.pie(
        sizes,
        startangle=90,
        colors=colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1},
        radius=1.0,
        center=(0, 0)  # Keep pie centered
    )
    ax.axis('equal')

    # Set title
    title = f"Read Distribution Across Genomic Features\n(Sample: {args.sampleName})"
    ax.set_title(title, fontsize=8, pad=8)

    # Position legend to the right of pie
    legend_labels = [f"{lbl} ({s*100:.1f}%)" for lbl, s in zip(labels, sizes)]
    ax.legend(wedges, legend_labels, loc='center left',
              bbox_to_anchor=(0.75, 0.5),  # Legend to the right of pie
              frameon=False, prop={'size': 6})

    # Save outputs (remove extra whitespace)
    pdf_path = os.path.join(args.output_dir, f"{args.sampleName}_reads_mapping_stats_pie.pdf")
    png_path = os.path.join(args.output_dir, f"{args.sampleName}_reads_mapping_stats_pie.png")

    fig.savefig(pdf_path, format='pdf', bbox_inches='tight', pad_inches=0.01)
    fig.savefig(png_path, format='png', bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)

    print(f"Plot saved successfully to: {png_path} and {pdf_path}")
