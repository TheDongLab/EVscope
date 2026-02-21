#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step_15_plot_top_expressed_genes.py

Generates a horizontal bar chart of the top expressed genes from an RNA
expression matrix, colored by RNA biotype.

Usage:
    python Step_15_plot_top_expressed_genes.py --input <expr_matrix> \\
        --meta <gene_meta> --output <output_prefix> [--top_n 50]
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import argparse
import matplotlib as mpl
import re

# Set PDF font type for vector graphics
mpl.rcParams['pdf.fonttype'] = 42

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Plot top-expressing genes by type, sorted by mean expression.')
    parser.add_argument('--input_gene_expr_matrix', required=True, help='Input TSV with gene expression.')
    parser.add_argument('--output_pdf', required=True, help='Output PDF file.')
    parser.add_argument('--output_png', required=True, help='Output PNG file.')
    parser.add_argument('--gene_num_per_type', type=int, required=True, help='Top genes per gene type.')
    parser.add_argument('--total_gene_num', type=int, required=True, help='Total genes to display.')
    args = parser.parse_args()

    # Read TSV file without forcing string dtype
    df = pd.read_csv(args.input_gene_expr_matrix, sep='\t', header=0)

    # Filter for expressed genes (Norm_Expr > 0)
    df = df[df['Norm_Expr'].astype(float) > 0]

    # Check if enough genes are available
    total_available = len(df)
    target_num = min(args.total_gene_num, total_available)

    # Select top N genes per GeneType
    N = args.gene_num_per_type
    df_sorted = df.sort_values('Norm_Expr', ascending=False)
    top_per_type = df_sorted.groupby('GeneType').head(N).reset_index(drop=True)

    # Get IDs of selected genes
    selected_gene_ids = top_per_type['GeneID'].tolist()

    # Get remaining genes
    remaining_df = df[~df['GeneID'].isin(selected_gene_ids)]

    # Supplement with top overall genes if needed
    if len(top_per_type) < target_num:
        additional_genes = remaining_df.nlargest(target_num - len(top_per_type), 'Norm_Expr')
        final_selection = pd.concat([top_per_type, additional_genes])
    else:
        final_selection = top_per_type.head(target_num)

    # Calculate mean Norm_Expr per GeneType
    mean_expr_per_type = final_selection.groupby('GeneType')['Norm_Expr'].mean().reset_index()
    mean_expr_per_type = mean_expr_per_type.sort_values('Norm_Expr', ascending=False)

    # Sort final_selection by GeneType (ordered by mean Norm_Expr) and within type by Norm_Expr
    type_order = mean_expr_per_type['GeneType'].tolist()
    final_selection['GeneType'] = pd.Categorical(final_selection['GeneType'], categories=type_order, ordered=True)
    final_selection = final_selection.sort_values(['GeneType', 'Norm_Expr'], ascending=[True, False])

    # Create gene labels: GeneSymbol(GeneID) or GeneSymbol without isoform suffix if base IDs match
    def create_label(row):
        # Clean strings by removing trailing whitespace and converting to string
        gene_symbol_clean = str(row['GeneSymbol']).strip().rstrip()
        gene_id_clean = str(row['GeneID']).strip().rstrip()
        # Extract base ID by removing isoform suffix for comparison
        base_symbol = re.sub(r'\.\d+$', '', gene_symbol_clean)
        base_id = re.sub(r'\.\d+$', '', gene_id_clean)
        # If base IDs match, use cleaned GeneSymbol without suffix; otherwise, use original format
        if base_symbol == base_id:
            return base_symbol
        else:
            return f"{row['GeneSymbol']}({row['GeneID']})"

    final_selection['label'] = final_selection.apply(create_label, axis=1)

    # Assign distinct colors for up to 25 GeneTypes
    unique_types = mean_expr_per_type['GeneType'].tolist()  # Sorted by mean expression
    num_types = len(unique_types)
    # Combine tab20 and tab20b for up to 40 distinct colors
    colors = plt.cm.tab20(np.linspace(0, 1, 20))[:20]
    colors = np.vstack([colors, plt.cm.tab20b(np.linspace(0, 1, 20))[:20]])
    type_to_color = {t: colors[i % len(colors)] for i, t in enumerate(unique_types)}
    color_list = [type_to_color[t] for t in final_selection['GeneType']]

    # Create horizontal bar plot
    fig, ax = plt.subplots(figsize=(8, 0.1 * len(final_selection)))  # Dynamic height for spacing
    y_pos = np.arange(len(final_selection))
    ax.barh(y_pos, final_selection['Norm_Expr'], color=color_list, height=0.5)  # Thinner bars for gaps
    ax.set_yticks(y_pos)
    ax.set_yticklabels(final_selection['label'], fontsize=6)  # Small font for many labels
    ax.invert_yaxis()  # Highest expression at top
    ax.set_xlabel('log10(Normalized Expression)')
    ax.set_xscale('log')  # Log10 scale for x-axis

    # Set two-line plot title
    ax.set_title(f"Barplot for top {len(final_selection)} highly expressed genes by gene type\n"
                 f"(Top {args.gene_num_per_type} per Type, sorted by mean expression level)",
                 fontsize=10, pad=10, loc='center')

    # Create legend with "Gene Type" title, single column
    legend_patches = [mpatches.Patch(color=type_to_color[t], label=t) for t in unique_types]
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.05, 1), loc='upper left', 
              title='Gene type', fontsize=8, title_fontsize=9, 
              framealpha=0.8, edgecolor='gray')

    # Minimize whitespace
    plt.margins(y=0)  # Reduce vertical margins
    plt.tight_layout(pad=0.5)  # Tight layout with minimal padding

    # Save outputs
    plt.savefig(args.output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(args.output_png, format='png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()
