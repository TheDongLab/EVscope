#!/usr/bin/env python3
"""
This script merges gene and circRNA expression matrices and plots RNA type composition.
It reads:
  - Gene expression TSV (columns: GeneID, GeneSymbol, GeneType, ReadCounts, TPM)
  - circRNA expression TSV (columns: circRNA_ID, junction_reads, junction_reads_CPM, circRNA_type, strand, junction_reads_ID)
For circRNAs, it sets:
  GeneID = circRNA_ID, GeneSymbol = circRNA_ID, GeneType = "circRNAs",
  ReadCounts = junction_reads, TPM = junction_reads_CPM.
The merged matrix is output with the header "TPM/CPM".

The plot is a stacked bar chart of RNA type counts passing various filters:
  Norm_Expr > 0.001, > 0.1, > 0.2, > 0.5, > 1,
  and ReadCounts > 1, > 5, > 10.
The fixed color mapping for 21 RNA types is used (reversed order in the legend).
"""

import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Merge gene and circRNA matrices and plot RNA type composition."
    )
    parser.add_argument("--gene_expr", required=True,
                        help="Path to gene expression TSV (GeneID, GeneSymbol, GeneType, ReadCounts, TPM)")
    parser.add_argument("--circRNA_expr", required=True,
                        help="Path to circRNA expression TSV (circRNA_ID, junction_reads, junction_reads_CPM, circRNA_type, strand, junction_reads_ID)")
    parser.add_argument("--out_matrix", required=True,
                        help="Output path for merged matrix (TSV)")
    parser.add_argument("--out_plot", required=True,
                        help="Output path for the RNA distribution plot (PDF or SVG)")
    return parser.parse_args()

def read_gene_expr(file_path):
    """Read gene expression TSV and return list of dicts."""
    data = []
    with open(file_path, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            try:
                row["ReadCounts"] = float(row["ReadCounts"])
                row["TPM"] = float(row["TPM"])
            except ValueError:
                continue
            data.append(row)
    return data

def read_circRNA_expr(file_path):
    """Read circRNA TSV and return list of dicts with fixed keys."""
    data = []
    with open(file_path, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            try:
                readcounts = float(row["junction_reads"])
                norm_expr = float(row["junction_reads_CPM"])
            except ValueError:
                continue
            circ_dict = {
                "GeneID": row["circRNA_ID"],
                "GeneSymbol": row["circRNA_ID"],
                "GeneType": "circRNAs",
                "ReadCounts": readcounts,
                "TPM": norm_expr
            }
            data.append(circ_dict)
    return data

def write_merged_matrix(merged_data, out_file):
    """Write merged expression matrix (gene + circRNA) to TSV."""
    fieldnames = ["GeneID", "GeneSymbol", "GeneType", "ReadCounts", "TPM/CPM"]
    with open(out_file, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in merged_data:
            output_row = {
                "GeneID": row["GeneID"],
                "GeneSymbol": row["GeneSymbol"],
                "GeneType": row["GeneType"],
                "ReadCounts": row["ReadCounts"],
                "TPM/CPM": row["TPM"]
            }
            writer.writerow(output_row)

def count_genes_by_filter(data, threshold, metric):
    """Count entries with metric value > threshold; returns dict {RNA type: count}."""
    counts = defaultdict(int)
    for row in data:
        try:
            value = float(row[metric])
        except (ValueError, KeyError):
            continue
        if value > threshold:
            counts[row.get("GeneType", "NA")] += 1
    return counts

def generate_counts_by_filters(merged_data):
    """Apply filters and return filter labels and counts per filter."""
    norm_expr_thresholds = [0.001, 0.1, 0.2, 0.5, 1]
    read_count_thresholds = [1, 5, 10]

    filter_labels = []
    counts_per_filter = []

    for t in norm_expr_thresholds:
        label = f"Norm_Expr > {t}"
        filter_labels.append(label)
        counts = count_genes_by_filter(merged_data, t, "TPM")
        counts_per_filter.append(counts)

    for t in read_count_thresholds:
        label = f"ReadCounts > {t}"
        filter_labels.append(label)
        counts = count_genes_by_filter(merged_data, t, "ReadCounts")
        counts_per_filter.append(counts)

    return filter_labels, counts_per_filter

def plot_stacked_bar_composition(filter_labels, counts_per_filter, out_plot_file):
    """Plot stacked bar chart with fixed colors (reversed order in legend)."""
    # Fixed color mapping (reversed order)
    fixed_color_map = {
        "vault_RNAs": "#66a61e",            
        "artifact": "#e7298a",              
        "IG_genes": "#7570b3",             
        "scaRNAs": "#d95f02",              
        "TR_genes": "#1b9e77",             
        "snoRNAs": "#b3b3b3",              
        "miRNAs": "#e5c494",              
        "TEC_protein_coding": "#ffd92f",   
        "rRNAs": "#a6d854",               
        "Y_RNAs": "#e78ac3",              
        "misc-sncRNAs": "#8da0cb",         
        "snRNAs": "#fc8d62",              
        "tRNAs": "#66c2a5",               
        "circRNAs": "#f781bf",            
        "pseudogenes": "#a65628",         
        "lncRNAs": "#ffff33",             
        "protein_coding": "#ff7f00",      
        "ERVs": "#984ea3",                
        "SINEs": "#4daf4a",               
        "LINEs": "#377eb8",               
        "piRNAs": "#e41a1c"               
    }

    num_filters = len(filter_labels)
    # Compute overall counts per RNA type (for sorting)
    overall_totals = defaultdict(int)
    for counts in counts_per_filter:
        for rna_type, cnt in counts.items():
            overall_totals[rna_type] += cnt

    # Sort RNA types by overall counts (desc)
    sorted_rna_types = sorted(overall_totals, key=lambda x: overall_totals[x], reverse=True)

    # Prepare counts per RNA type per filter
    rna_type_counts = {}
    for rna in sorted_rna_types:
        rna_type_counts[rna] = [counts.get(rna, 0) for counts in counts_per_filter]

    bar_width = 0.5
    group_gap = 0.2
    x = np.arange(num_filters) * (bar_width + group_gap)

    fig, ax = plt.subplots(figsize=(11.69, 8.27))
    bottom = np.zeros(num_filters)

    for rna_type in sorted_rna_types:
        segment = np.array(rna_type_counts[rna_type])
        if np.all(segment == 0):
            continue
        label_str = f"{rna_type} ({', '.join(str(val) for val in segment)})"
        color = fixed_color_map.get(rna_type, "#000000")
        ax.bar(x, segment, bar_width, bottom=bottom, label=label_str, color=color)
        bottom += segment

    ax.set_xticks(x)
    ax.set_xticklabels(filter_labels, rotation=-45, ha="left", fontsize=10, fontweight='bold')
    ax.set_ylabel("Number of expressed genes/RNAs", fontsize=12, fontweight='bold')
    ax.set_title("RNA type composition under expression filters\n(Norm_Expr: linear_RNAs / CircRNAs: TPM/CPM)", fontsize=14, fontweight='bold')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    leg = ax.legend(
        title="RNA types and counts",
        fontsize=10,
        title_fontsize=12,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0)
    )
    leg.get_title().set_fontweight('bold')
    fig.tight_layout()
    plt.savefig(out_plot_file, dpi=500, bbox_inches='tight')
    plt.close()

def main():
    args = parse_arguments()
    gene_data = read_gene_expr(args.gene_expr)
    circRNA_data = read_circRNA_expr(args.circRNA_expr)
    merged_data = gene_data + circRNA_data
    write_merged_matrix(merged_data, args.out_matrix)
    print(f"Merged matrix written to {args.out_matrix}")

    filter_labels, counts_per_filter = generate_counts_by_filters(merged_data)
    plot_stacked_bar_composition(filter_labels, counts_per_filter, args.out_plot)
    print(f"RNA distribution plot saved to {args.out_plot}")

if __name__ == "__main__":
    main()

