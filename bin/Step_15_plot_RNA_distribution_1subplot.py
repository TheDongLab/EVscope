#!/usr/bin/env python3
"""Generate RNA composition plot from merged matrix, split into two subplots (Norm_Expr vs ReadCounts)."""

import argparse
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate RNA composition plot from merged matrix, split into two subplots."
    )
    parser.add_argument("--sample_name", required=True, help="Sample name for figure title")
    parser.add_argument("--Expr_matrix", required=True,
                        help="Merged matrix TSV (GeneID, GeneSymbol, GeneType, ReadCounts, Norm_Expr)")
    parser.add_argument("--out_plot", required=True,
                        help="Output plot path (will save .pdf, .png)")
    return parser.parse_args()

def read_merged_matrix(file_path):
    """Read merged TSV matrix into a list of dicts."""
    data = []
    with open(file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                row["ReadCounts"] = float(row["ReadCounts"])
                row["Norm_Expr"] = float(row["Norm_Expr"])
            except ValueError:
                continue
            data.append(row)
    return data

def count_genes_by_filter(data, threshold, metric):
    """Count how many entries (by GeneType) pass a given threshold in 'metric'."""
    counts = defaultdict(int)
    for row in data:
        try:
            if float(row[metric]) > threshold:
                counts[row["GeneType"]] += 1
        except (ValueError, KeyError):
            continue
    return counts

def generate_counts(data):
    """
    Return:
        filter_labels: list of filter labels
        counts_per_filter: list of dict, each dict maps GeneType->count
    """
    norm_thresholds = [0.001, 0.1, 0.2, 0.5, 1]
    # Updated read_thresholds to include 1,5,10,20,50
    read_thresholds = [1, 5, 10, 20, 50]

    filter_labels = []
    counts_per_filter = []

    # Norm_Expr filters
    for t in norm_thresholds:
        filter_labels.append(f"Norm_Expr > {t}")
        counts_per_filter.append(count_genes_by_filter(data, t, "Norm_Expr"))

    # ReadCounts filters
    for t in read_thresholds:
        filter_labels.append(f"ReadCounts > {t}")
        counts_per_filter.append(count_genes_by_filter(data, t, "ReadCounts"))

    return filter_labels, counts_per_filter

def plot_stacked_bar(ax, filter_labels, counts_per_filter, sorted_types, color_map, subplot_title):
    """
    Plot a stacked bar composition on a given Axes.
    - filter_labels: x-axis categories
    - counts_per_filter: list of dicts (GeneType->count)
    - sorted_types: RNA types sorted by total
    - color_map: dict for coloring
    - subplot_title: title for this subplot
    """
    # Prepare data
    plot_data = {rt: [c.get(rt, 0) for c in counts_per_filter] for rt in sorted_types}

    x = np.arange(len(filter_labels)) * 0.7
    bottom = np.zeros(len(filter_labels))

    for rt in sorted_types:
        counts = plot_data[rt]
        if sum(counts) == 0:
            continue
        ax.bar(
            x, counts, 0.5, bottom=bottom,
            color=color_map.get(rt, "#000000"),
            label=f"{rt} ({','.join(map(str, counts))})"
        )
        bottom += counts

    ax.set_xticks(x)
    ax.set_xticklabels(filter_labels, rotation=-45, ha="left", fontsize=10, fontweight='bold')
    ax.set_ylabel("Number of expressed genes/RNAs", fontsize=12, fontweight='bold')
    ax.set_title(subplot_title, fontsize=12, fontweight='bold')

    # Legend
    leg = ax.legend(
        title="RNA types and counts", fontsize=10,
        title_fontsize=12, loc="upper left",
        bbox_to_anchor=(1.01, 1.0)
    )
    leg.get_title().set_fontweight('bold')

def main():
    args = parse_arguments()

    # Read data
    merged_data = read_merged_matrix(args.Expr_matrix)
    filter_labels, counts_per_filter = generate_counts(merged_data)

    # Since each set (Norm_Expr vs ReadCounts) now has 5 filters each, we can split them 0..5 and 5..10
    norm_expr_labels = filter_labels[:5]
    norm_expr_counts = counts_per_filter[:5]
    read_counts_labels = filter_labels[5:]
    read_counts_counts = counts_per_filter[5:]

    # Determine all RNA types and their total usage
    totals = defaultdict(int)
    for cdict in counts_per_filter:
        for rna_type, val in cdict.items():
            totals[rna_type] += val
    sorted_types = sorted(totals, key=lambda x: totals[x], reverse=True)

    # Define color map
    color_map = {
        "vault_RNAs": "#66a61e", "IG_genes": "#7570b3",
        "scaRNAs": "#d95f02", "TR_genes": "#1b9e77",
        "snoRNAs": "#b3b3b3", "miRNAs": "#e5c494",
        "TEC_protein_coding": "#ffd92f", "rRNAs": "#a6d854",
        "Y_RNAs": "#e78ac3", "misc-sncRNAs": "#8da0cb",
        "snRNAs": "#fc8d62", "tRNAs": "#66c2a5",
        "circRNAs": "#f781bf", "pseudogenes": "#a65628",
        "lncRNAs": "#ffff33", "protein_coding": "#ff7f00",
        "ERVs": "#984ea3", "SINEs": "#4daf4a",
        "LINEs": "#377eb8", "piRNAs": "#e41a1c"
    }

    # Create figure with two subplots
    fig, axes = plt.subplots(2, 1, figsize=(8, 12))
    fig.suptitle(
        f"{args.sample_name}: RNA type composition under expression filters\n"
        "(Norm_Expr: linear RNAs or circRNAs in TPM/CPM; ReadCounts: raw read counts)",
        fontsize=14, fontweight='bold'
    )

    # Top subplot: Norm_Expr
    plot_stacked_bar(
        axes[0], norm_expr_labels, norm_expr_counts,
        sorted_types, color_map, "Norm_Expr Filters"
    )
    # Bottom subplot: ReadCounts
    plot_stacked_bar(
        axes[1], read_counts_labels, read_counts_counts,
        sorted_types, color_map, "ReadCounts Filters"
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Derive base name for saving multiple formats
    base_name = args.out_plot.rsplit('.', 1)[0]
    for ext in ["pdf", "png"]:
        out_path = f"{base_name}.{ext}"
        plt.savefig(out_path, dpi=500, bbox_inches='tight')
        print(f"Saved figure to {out_path}")

    plt.close()

if __name__ == "__main__":
    main()

