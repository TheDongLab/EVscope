#!/usr/bin/env python
"""
Compute total read counts from featureCounts TSV files and plot a pie chart.
Legend order: 5'UTR, Exon, 3'UTR, Intron, Promoter, Downstream 2Kb, Intergenic, ENCODE_blacklist.
"""
import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42

def sum_read_counts(tsv_file):
    """Sum counts in last column of TSV, skip headers and comments."""
    total = 0.0
    with open(tsv_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('Geneid'):
                continue
            parts = line.strip().split('\t')
            try:
                total += float(parts[-1])
            except ValueError:
                pass
    return total

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot read mapping stats pie chart.'
    )
    parser.add_argument('--input_5UTR_readcounts', required=True,
                        help="TSV for 5'UTR counts")
    parser.add_argument('--input_exon_readcounts', required=True,
                        help='TSV for exon counts')
    parser.add_argument('--input_3UTR_readcounts', required=True,
                        help="TSV for 3'UTR counts")
    parser.add_argument('--input_intron_readcounts', required=True,
                        help='TSV for intron counts')
    parser.add_argument('--input_promoters_readcounts', required=True,
                        help='TSV for promoter counts')
    parser.add_argument('--input_downstream_2Kb_readcounts', required=True,
                        help='TSV for downstream counts')
    parser.add_argument('--input_intergenic_readcounts', required=True,
                        help='TSV for intergenic counts')
    parser.add_argument('--input_ENCODE_blacklist_readcounts', required=True,
                        help='TSV for ENCODE blacklist counts')
    parser.add_argument('--sampleName', required=True,
                        help='Sample name for output')
    parser.add_argument('--output_dir', required=True,
                        help='Directory to save output files')
    args = parser.parse_args()

    # create output directory if needed
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)

    # sum reads for each region
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

    total = sum(counts)
    if total == 0:
        sys.exit('Error: total read count is zero.')

    # prepare pie chart data
    sizes = [c / total for c in counts]
    labels = ["5'UTR", 'Exon', "3'UTR", 'Intron', 'Promoter',
              'Downstream 2 Kb', 'Intergenic', 'ENCODE_blacklist']
    colors = ['red', 'orange', 'cyan', 'green', 'purple', 'blue', 'pink', 'gray']

    # plot pie
    fig, ax = plt.subplots(figsize=(3, 3), dpi=300)
    wedges, _ = ax.pie(
        sizes,
        startangle=90,
        colors=colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1},
        textprops={'fontsize': 6},
        radius=0.2
    )
    ax.axis('equal')
    plt.title(f"Read mapping in {args.sampleName}", fontsize=8)

    # legend with percentages
    legend_labels = [f"{lbl} ({s*100:.1f}%)" for lbl, s in zip(labels, sizes)]
    ax.legend(wedges, legend_labels, loc='center left',
              bbox_to_anchor=(1, 0.5), frameon=False, prop={'size': 6})

    plt.tight_layout()

    # save outputs
    pdf_path = os.path.join(args.output_dir, f"{args.sampleName}_reads_mapping_stats_pie.pdf")
    png_path = os.path.join(args.output_dir, f"{args.sampleName}_reads_mapping_stats_pie.png")
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight')
    fig.savefig(png_path, format='png', bbox_inches='tight')
    plt.close(fig)

    print(f"Outputs saved: {pdf_path}, {png_path}")

