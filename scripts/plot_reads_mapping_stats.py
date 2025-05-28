#!/usr/bin/env python
"""
Compute total read counts from featureCounts TSV files and plot a pie chart.
Legend is placed on the right side with entries (top to bottom):
5'UTR, exon, 3'UTR, intron, promoter, downstream, intergenic, ENCODE_blacklist.
"""

import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl

# Function to sum read counts from a TSV file (ignores header lines)
def sum_read_counts(tsv_file):
    total = 0.0
    with open(tsv_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comment lines and header line
            if line.startswith('#') or line.startswith("Geneid"):
                continue
            parts = line.split('\t')
            try:
                count = float(parts[-1])
                total += count
            except ValueError:
                continue
    return total

def main():
    parser = argparse.ArgumentParser(
        description="Plot read mapping stats pie chart from featureCounts TSV files."
    )
    parser.add_argument("--input_5UTR_readcounts", required=True,
                        help="TSV file for 5'UTR read counts")
    parser.add_argument("--input_exon_readcounts", required=True,
                        help="TSV file for exon read counts")
    parser.add_argument("--input_3UTR_readcounts", required=True,
                        help="TSV file for 3'UTR read counts")
    parser.add_argument("--input_intron_readcounts", required=True,
                        help="TSV file for intron read counts")
    parser.add_argument("--input_promoters_readcounts", required=True,
                        help="TSV file for promoter read counts")
    parser.add_argument("--input_downstream_2Kb_readcounts", required=True,
                        help="TSV file for downstream read counts")
    parser.add_argument("--input_intergenic_readcounts", required=True,
                        help="TSV file for intergenic read counts")
    parser.add_argument("--input_ENCODE_blacklist_readcounts", required=True,
                        help="TSV file for ENCODE_blacklist read counts")
    # Added sampleName parameter for custom output filenames and figure title
    parser.add_argument("--sampleName", required=True,
                        help="Sample name used in output filenames and figure title")
    args = parser.parse_args()

    # Sum reads from each TSV file
    count_5UTR       = sum_read_counts(args.input_5UTR_readcounts)
    count_exon       = sum_read_counts(args.input_exon_readcounts)
    count_3UTR       = sum_read_counts(args.input_3UTR_readcounts)
    count_intron     = sum_read_counts(args.input_intron_readcounts)
    count_promoter   = sum_read_counts(args.input_promoters_readcounts)
    count_downstream = sum_read_counts(args.input_downstream_2Kb_readcounts)
    count_intergenic = sum_read_counts(args.input_intergenic_readcounts)
    count_blacklist  = sum_read_counts(args.input_ENCODE_blacklist_readcounts)

    # Compute total reads
    total_reads = (count_5UTR + count_exon + count_3UTR + count_intron +
                   count_promoter + count_downstream + count_intergenic + count_blacklist)
    if total_reads == 0:
        print("Total reads count is zero, cannot plot pie chart.")
        sys.exit(1)

    # Compute fraction for each region
    sizes = [
        count_5UTR / total_reads,
        count_exon / total_reads,
        count_3UTR / total_reads,
        count_intron / total_reads,
        count_promoter / total_reads,
        count_downstream / total_reads,
        count_intergenic / total_reads,
        count_blacklist / total_reads,
    ]
    # Order of labels as specified
    labels = ["5'UTR", "Exon", "3'UTR", "Intron", "Promoter", "Downstream 2 Kb", "Intergenic", "ENCODE_blacklist"]

    # Set style (SVG vector fonts and Arial)
    mpl.rcParams['svg.fonttype'] = 'none'
    plt.rcParams["font.family"] = "Arial"

    # Define eight colors for the eight regions
    colors = ['red', 'orange', 'cyan', 'green', 'purple', 'blue', 'pink', 'gray']

    # Create pie chart with same style as original code
    fig, ax = plt.subplots(figsize=(3, 3), dpi=300)
    wedges, texts = ax.pie(
        sizes,
        startangle=90,
        colors=colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1},
        textprops={'fontsize': 6},
        radius=0.2
    )
    ax.axis('equal')
    
    # Set the figure title with sample name at the top center
    plt.title(f"Percentage of reads mapping to meta-gene regions in {args.sampleName}", fontsize=8)

    # Create legend labels with percentages
    legend_labels = [f"{lbl} ({s*100:.1f}%)" for lbl, s in zip(labels, sizes)]
    # Place legend on the right of the pie chart
    ax.legend(wedges, legend_labels, loc='center left', bbox_to_anchor=(1, 0.5),
              ncol=1, frameon=False, prop={'size': 6})

    plt.tight_layout()

    # Save pie chart in PDF, PNG, and SVG formats using sampleName in filenames
    current_dir = os.getcwd()
    pie_pdf = os.path.join(current_dir, f"{args.sampleName}_reads_mapping_stats_pie.pdf")
    pie_png = os.path.join(current_dir, f"{args.sampleName}_reads_mapping_stats_pie.png")
    pie_svg = os.path.join(current_dir, f"{args.sampleName}_reads_mapping_stats_pie.svg")
    plt.savefig(pie_pdf, format='pdf', bbox_inches='tight')
    plt.savefig(pie_png, format='png', bbox_inches='tight')
    plt.savefig(pie_svg, format='svg', bbox_inches='tight')
    plt.close()

    print(f"Pie chart saved as: {pie_pdf}, {pie_png}, {pie_svg}")

if __name__ == "__main__":
    main()

