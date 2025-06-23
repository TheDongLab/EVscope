#!/usr/bin/env python3
"""Generate RNA composition barplots split by Norm_Expr and ReadCounts filters"""

import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
from collections import defaultdict
import os

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42


# Configuration
CUSTOM_ORDER = [
    "SINEs", "LINEs", "ERVs", "piRNAs", "protein_coding",
    "lncRNAs", "pseudogenes", "misc-sncRNAs", "miRNAs",
    "TEC_protein_coding", "snRNAs", "IG_genes", "snoRNAs",
    "tRNAs", "Y_RNAs", "rRNAs", "TR_genes", "circRNAs",
    "scaRNAs", "vault_RNAs"
]
UNIFIED_COLOR = "#98DF8A"
LABEL_FONT = {'family': 'Arial', 'weight': 'bold', 'size': 9}
VALUE_COLOR = '#0000FF'

# Filter labels (5 thresholds each)
FILTER_LABELS_READCOUNTS = [
    "ReadCounts > 50", "ReadCounts > 20", "ReadCounts > 10", "ReadCounts > 5", "ReadCounts > 1"
]
FILTER_LABELS_NORM_EXPR = [
    "Norm_Expr > 0.001", "Norm_Expr > 0.1", "Norm_Expr > 0.2", "Norm_Expr > 0.5", "Norm_Expr > 1"
]

def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_name", required=True)
    parser.add_argument("--Expr_matrix", required=True)
    parser.add_argument("--out_plot", required=True,
                        help="Output PDF file name; PNG will be derived")
    return parser.parse_args()

def k_formatter(x, pos):
    """Format x-axis values with K notation"""
    if x == 0:
        return "0"
    if x < 1000:
        return f"{int(x)}"
    value = x / 1000
    return f"{value:.1f}K" if value != int(value) else f"{int(value)}K"

def generate_counts(data):
    """Calculate counts for each metric and threshold"""
    results = {}
    # ReadCounts thresholds: 50, 20, 10, 5, 1
    readcounts_thresholds = [50, 20, 10, 5, 1]
    readcounts_counters = []
    for thresh in readcounts_thresholds:
        counter = defaultdict(int)
        for entry in data:
            if float(entry["ReadCounts"]) > thresh:
                counter[entry["GeneType"]] += 1
        readcounts_counters.append(counter)
    results["ReadCounts"] = readcounts_counters

    # Norm_Expr thresholds: 0.001, 0.1, 0.2, 0.5, 1
    norm_expr_thresholds = [0.001, 0.1, 0.2, 0.5, 1]
    norm_expr_counters = []
    for thresh in norm_expr_thresholds:
        counter = defaultdict(int)
        for entry in data:
            if float(entry["Norm_Expr"]) > thresh:
                counter[entry["GeneType"]] += 1
        norm_expr_counters.append(counter)
    results["Norm_Expr"] = norm_expr_counters

    return results

def plot_bars(plot_data, sample, out_file):
    """Generate two-part barplots: top row for Norm_Expr and bottom row for ReadCounts"""
    plt.style.use('default')
    # Set figure size for 2 rows x 20 columns of subplots
    fig = plt.figure(figsize=(34, 12))
    fig.suptitle(
        f"{sample}: RNA Composition Under Filters\n"
        "Top: Norm_Expr filtering | Bottom: ReadCounts filtering\n"
        "(Norm_Expr: linear_RNAs / CircRNAs: TPM/CPM)",
        y=0.98, fontsize=14, fontweight='bold'
    )

    # Create GridSpec with 2 rows and 20 columns
    gs = plt.GridSpec(2, 20, figure=fig, hspace=0.5, wspace=0.05)

    # Loop through rows (0: Norm_Expr, 1: ReadCounts) and columns (RNA types)
    for row in range(2):
        metric = "Norm_Expr" if row == 0 else "ReadCounts"
        filter_labels = FILTER_LABELS_NORM_EXPR if row == 0 else FILTER_LABELS_READCOUNTS

        for idx, rna_type in enumerate(CUSTOM_ORDER):
            ax = fig.add_subplot(gs[row, idx])
            counts = plot_data[rna_type][metric]
            max_val = max(counts) if counts else 1

            # Set y positions for 5 bars
            y_pos = np.arange(5) * 0.6
            ax.barh(y_pos, counts, height=0.47, color=UNIFIED_COLOR)

            # Add value labels on bars
            for i, val in enumerate(counts):
                if val > 0:
                    ax.text(0.1, y_pos[i], f'({int(val)})',
                            ha='left', va='center',
                            color=VALUE_COLOR, fontdict=LABEL_FONT)

            ax.set_xlim(0, max_val * 1.25)
            ax.xaxis.set_major_formatter(FuncFormatter(k_formatter))
            ax.set_xlabel(rna_type, rotation=-45, ha='left', **LABEL_FONT)
            ax.grid(axis='x', alpha=0.2, linestyle='--')
            ax.spines[['right', 'top']].set_visible(False)

            # Show y-axis labels only in the first column
            if idx == 0:
                ax.set_yticks(y_pos)
                ax.set_yticklabels(filter_labels, fontdict=LABEL_FONT)
            else:
                ax.set_yticks([])

    # Save figure in PDF, PNG formats
    pdf_file = out_file
    base, _ = os.path.splitext(pdf_file)
    png_file = base + ".png"
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    args = parse_args()

    # Read merged expression data
    with open(args.Expr_matrix) as f:
        reader = csv.DictReader(f, delimiter='\t')
        merged_data = list(reader)

    # Generate counts for each threshold
    threshold_counts = generate_counts(merged_data)

    # Build plot data: for each RNA type, get counts for each metric
    plot_data = {
        rt: {
            "ReadCounts": [cnt.get(rt, 0) for cnt in threshold_counts["ReadCounts"]],
            "Norm_Expr": [cnt.get(rt, 0) for cnt in threshold_counts["Norm_Expr"]]
        } for rt in CUSTOM_ORDER
    }

    # Create plots and save files
    plot_bars(plot_data, args.sample_name, args.out_plot)

if __name__ == "__main__":
    main()

