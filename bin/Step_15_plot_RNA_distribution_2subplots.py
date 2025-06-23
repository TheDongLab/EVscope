#!/usr/bin/env python3
"""
Optimized version with improved layout:
 - Splits filters into two groups (Norm_Expr and ReadCounts) each with 5 thresholds.
 - Creates one figure with 8 subplots (4 rows x 2 columns).
 - Top two rows: Norm_Expr (counts row, percentage row).
 - Bottom two rows: ReadCounts (counts row, percentage row).
 - Legend placed on the right side to avoid overlap with plots.
 - Reduced spacing between subplots for better grouping of same-category plots.
 - Font set to Arial, and pdf/png output.
"""

import argparse
import csv
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from collections import defaultdict

# Set global font to Arial
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams['pdf.fonttype'] = 42

# ---------------- Configuration ----------------
SMALL_RNA_TYPES = ["vault_RNAs", "scaRNAs", "snoRNAs", "miRNAs", "Y_RNAs",
                   "misc-sncRNAs", "snRNAs", "tRNAs", "piRNAs"]

LONG_RNA_TYPES = ["IG_genes", "TR_genes", "TEC_protein_coding", "rRNAs",
                  "circRNAs", "pseudogenes", "lncRNAs", "LINEs", "SINEs",
                  "ERVs", "protein_coding"]

COLOR_MAP = {
    # Small RNAs (9)
    "vault_RNAs": "#1f77b4",
    "scaRNAs": "#aec7e8",
    "snoRNAs": "#ff7f0e",
    "miRNAs": "#ffbb78",
    "Y_RNAs": "#2ca02c",
    "misc-sncRNAs": "#98df8a",
    "snRNAs": "#d62728",
    "tRNAs": "#ff9896",
    "piRNAs": "#9467bd",

    # Long RNAs (11)
    "IG_genes": "#c5b0d5",
    "TR_genes": "#8c564b",
    "TEC_protein_coding": "#c49c94",
    "rRNAs": "#e377c2",
    "circRNAs": "#f7b6d2",
    "pseudogenes": "#7f7f7f",
    "lncRNAs": "#c7c7c7",
    "LINEs": "#bcbd22",
    "SINEs": "#dbdb8d",
    "ERVs": "#17becf",
    "protein_coding": "#9edae5"
}

NE_THRESHOLDS = [0.001, 0.1, 0.2, 0.5, 1]
NE_FILTER_LABELS = [
    "Norm_Expr > 0.001",
    "Norm_Expr > 0.1",
    "Norm_Expr > 0.2",
    "Norm_Expr > 0.5",
    "Norm_Expr > 1"
]

RC_THRESHOLDS = [50, 20, 10, 5, 1]
RC_FILTER_LABELS = [
    "ReadCounts > 50",
    "ReadCounts > 20",
    "ReadCounts > 10",
    "ReadCounts > 5",
    "ReadCounts > 1"
]

# ---------------- Core Functions ----------------
def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate RNA composition plots by filter type.")
    parser.add_argument("--Expr_matrix", required=True, help="Merged RNA expression matrix TSV")
    parser.add_argument("--out_plot", required=True, help="Output file prefix (no extension)")
    parser.add_argument("--sample_name", required=True, help="Sample name for titles")
    return parser.parse_args()

def read_merged_matrix(file_path):
    data = []
    with open(file_path, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            try:
                row["ReadCounts"] = float(row["ReadCounts"])
                row["Norm_Expr"] = float(row["Norm_Expr"])
            except ValueError:
                continue
            data.append(row)
    return data

def generate_filter_data_by_metric(data, metric, thresholds):
    filter_data = []
    for thresh in thresholds:
        type_counts = defaultdict(int)
        for row in data:
            if row[metric] > thresh:
                type_counts[row["GeneType"]] += 1
        filter_data.append(type_counts)
    return filter_data

def sort_rna_types(filter_data_list, rna_types):
    totals = {rna: sum(d.get(rna, 0) for d in filter_data_list) for rna in rna_types}
    sorted_types = sorted(rna_types, key=lambda r: totals[r], reverse=True)
    return sorted_types

def plot_side_counts(ax, filter_data, sorted_types, labels, is_left=True):
    n_filters = len(labels)
    y_positions = np.arange(n_filters)
    bar_height = 0.8
    bottom = np.zeros(n_filters)
    for rna in sorted_types:
        counts = [d.get(rna, 0) for d in filter_data]
        ax.barh(y_positions, counts, left=bottom, height=bar_height,
                color=COLOR_MAP[rna], edgecolor='white')
        bottom += counts
    ax.set_xlim(0, bottom.max()*1.15 if bottom.max() > 0 else 1)
    if is_left:
        ax.invert_xaxis()
    ax.grid(axis='x', alpha=0.3)
    ax.set_facecolor('#F5F5F5')
    ax.set_yticks(y_positions + bar_height/2)
    if is_left:
        ax.set_yticklabels(labels, fontsize=12, fontweight='bold')
        ax.tick_params(axis='y', pad=5)
    else:
        ax.set_yticklabels([])

def plot_side_percentage(ax, filter_data, sorted_types, labels, is_left=True):
    n_filters = len(labels)
    y_positions = np.arange(n_filters)
    bar_height = 0.8
    totals = [sum(d.get(r, 0) for r in sorted_types) for d in filter_data]
    totals = [t if t > 0 else 1 for t in totals]
    bottom = np.zeros(n_filters)
    for rna in sorted_types:
        counts = [d.get(rna, 0) for d in filter_data]
        percents = [(count/total)*100 for count, total in zip(counts, totals)]
        ax.barh(y_positions, percents, left=bottom, height=bar_height,
                color=COLOR_MAP[rna], edgecolor='white')
        bottom = [b + p for b, p in zip(bottom, percents)]
    ax.set_xlim(0, 100)
    if is_left:
        ax.invert_xaxis()
    ax.grid(axis='x', alpha=0.3)
    ax.set_facecolor('#F5F5F5')
    ax.set_yticks(y_positions + bar_height/2)
    if is_left:
        ax.set_yticklabels(labels, fontsize=12, fontweight='bold')
        ax.tick_params(axis='y', pad=5)
    else:
        ax.set_yticklabels([])
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])

def get_legend_handles(sorted_long, sorted_small):
    dummy_long = Line2D([], [], color='none', label='Long RNA Types')
    dummy_small = Line2D([], [], color='none', label='Small RNA Types')
    long_handles = [Patch(facecolor=COLOR_MAP[rt], label=rt) for rt in sorted_long]
    small_handles = [Patch(facecolor=COLOR_MAP[rt], label=rt) for rt in sorted_small]
    empty_handle = Patch(facecolor='none', label='')
    return [dummy_long] + long_handles + [empty_handle] + [dummy_small] + small_handles

def create_combined_plots(filter_data_ne, filter_data_rc, out_prefix, sample_name):
    # Combine filter data to determine sorted order for legend
    combined_data = filter_data_ne + filter_data_rc
    sorted_long = sort_rna_types(combined_data, LONG_RNA_TYPES)
    sorted_small = sort_rna_types(combined_data, SMALL_RNA_TYPES)

    # Figure with 4 rows x 2 columns
    fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(20, 28),
                            gridspec_kw={'hspace': 0.3, 'wspace': 0.15})
    # Adjust margins so legend can sit on the right without overlapping
    fig.subplots_adjust(left=0.08, right=0.78, top=0.88, bottom=0.06)

    # -- Upper block: Norm_Expr (rows 0,1) --
    # Row 0: Counts
    plot_side_counts(axs[0,0], filter_data_ne, sorted_long, NE_FILTER_LABELS, is_left=True)
    axs[0,0].set_xlabel("Long RNAs", fontsize=14, fontweight='bold', labelpad=10)
    plot_side_counts(axs[0,1], filter_data_ne, sorted_small, NE_FILTER_LABELS, is_left=False)
    axs[0,1].set_xlabel("Small RNAs", fontsize=14, fontweight='bold', labelpad=10)

    # Row 1: Percentage
    plot_side_percentage(axs[1,0], filter_data_ne, sorted_long, NE_FILTER_LABELS, is_left=True)
    axs[1,0].set_xlabel("Long RNAs", fontsize=14, fontweight='bold', labelpad=10)
    plot_side_percentage(axs[1,1], filter_data_ne, sorted_small, NE_FILTER_LABELS, is_left=False)
    axs[1,1].set_xlabel("Small RNAs", fontsize=14, fontweight='bold', labelpad=10)

    # -- Lower block: ReadCounts (rows 2,3) --
    # Row 2: Counts
    plot_side_counts(axs[2,0], filter_data_rc, sorted_long, RC_FILTER_LABELS, is_left=True)
    axs[2,0].set_xlabel("Long RNAs", fontsize=14, fontweight='bold', labelpad=10)
    plot_side_counts(axs[2,1], filter_data_rc, sorted_small, RC_FILTER_LABELS, is_left=False)
    axs[2,1].set_xlabel("Small RNAs", fontsize=14, fontweight='bold', labelpad=10)

    # Row 3: Percentage
    plot_side_percentage(axs[3,0], filter_data_rc, sorted_long, RC_FILTER_LABELS, is_left=True)
    axs[3,0].set_xlabel("Long RNAs", fontsize=14, fontweight='bold', labelpad=10)
    plot_side_percentage(axs[3,1], filter_data_rc, sorted_small, RC_FILTER_LABELS, is_left=False)
    axs[3,1].set_xlabel("Small RNAs", fontsize=14, fontweight='bold', labelpad=10)

    # -- Titles --
    fig.text(0.5, 0.92,
             f"{sample_name}: RNA type composition under Norm_Expr filters (TPM/CPM)",
             ha='center', fontsize=18, fontweight='bold')
    fig.text(0.5, 0.48,
             f"{sample_name}: RNA type composition under ReadCounts filters",
             ha='center', fontsize=18, fontweight='bold')

    # -- Legend on the right --
    handles = get_legend_handles(sorted_long, sorted_small)
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.82, 0.5),
               fontsize=12, frameon=False)

    # Save in multiple formats
    pdf_out = out_prefix + ".pdf"
    png_out = out_prefix + ".png"
    plt.savefig(pdf_out, dpi=600, bbox_inches='tight')
    plt.savefig(png_out, dpi=600, bbox_inches='tight')
    plt.close()

    print(f"Plots saved:\n {pdf_out}\n {png_out}\n")

def main():
    args = parse_arguments()
    merged_data = read_merged_matrix(args.Expr_matrix)

    # Generate filter data for Norm_Expr
    filter_data_ne = generate_filter_data_by_metric(merged_data, "Norm_Expr", NE_THRESHOLDS)
    # Generate filter data for ReadCounts
    filter_data_rc = generate_filter_data_by_metric(merged_data, "ReadCounts", RC_THRESHOLDS)

    create_combined_plots(filter_data_ne, filter_data_rc, args.out_plot, args.sample_name)
    print("Final plot generated.")

if __name__ == "__main__":
    main()

