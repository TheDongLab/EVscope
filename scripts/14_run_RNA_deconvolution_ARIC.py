#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RNA Deconvolution using ARIC with Plotting

This script performs RNA deconvolution using the ARIC algorithm and generates plots for
visualizing the deconvolution results. It accepts:
  1. A bulk RNA expression CSV file containing exactly two columns:
     - Column 1: Gene identifier (either Gene Symbol or Ensembl ID). Ensure there are no extra spaces or tabs.
     - Column 2: Numeric expression values.
  2. A reference expression CSV file containing at least two columns:
     - One or two candidate gene identifier columns (non-numeric) followed by numeric expression columns.
     - The candidate gene identifier columns can be either Gene Symbols or Ensembl IDs.

Optional:
  --bulk_gene_type: Specify the gene identifier type for bulk RNA ("Gene_symbol" or "Ensembl_ID").
                   If not provided, the type is determined automatically by comparing gene overlaps.
  --bar_height:      (Optional) Set the bar height in the plot (default: 1).
"""

import argparse
import os
import re
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter  # for setting tick intervals and formatting
from ARIC import ARIC  # ARIC deconvolution function

mpl.rcParams['svg.fonttype'] = 'none'
def strip_ensembl_version(gene_id: str) -> str:
    """Remove version numbers from an Ensembl gene identifier."""
    return re.sub(r'\.\d+$', '', gene_id)

def is_gene_identifier_column(series: pd.Series) -> bool:
    """Determine if a Series represents a gene identifier column (mostly non-numeric)."""
    numeric_count = pd.to_numeric(series, errors='coerce').notnull().sum()
    return numeric_count < len(series)

def determine_bulk_gene_identifier_type(bulk_expression_file: str, ref_expression_file: str) -> (str, float):
    """
    Automatically determine the gene identifier type for bulk RNA data by comparing gene overlaps
    with candidate gene identifier columns from the reference expression file.
    Returns a tuple (gene_identifier_type, overlap_ratio).
    """
    bulk_df = pd.read_csv(bulk_expression_file, usecols=[0])
    bulk_genes = bulk_df.iloc[:, 0].astype(str).apply(
        lambda x: strip_ensembl_version(x.strip()) if re.match(r'^ENSG\d+', x.strip()) else x.strip()
    )
    bulk_gene_set = set(bulk_genes)

    ref_df = pd.read_csv(ref_expression_file, usecols=[0, 1])
    candidate_results = {}

    for col in ref_df.columns:
        col_data = ref_df[col].astype(str).apply(lambda x: x.strip())
        if is_gene_identifier_column(col_data):
            if re.match(r'^ENSG', col_data.iloc[0]):
                ref_genes = col_data.apply(lambda x: strip_ensembl_version(x))
            else:
                ref_genes = col_data
            ref_gene_set = set(ref_genes)
            overlap = len(bulk_gene_set & ref_gene_set)
            overlap_ratio = overlap / len(bulk_gene_set) if bulk_gene_set else 0
            candidate_results[col] = (overlap, overlap_ratio)
        else:
            candidate_results[col] = (0, 0)

    for col, (overlap, ratio) in candidate_results.items():
        print(f"Reference column '{col}': {overlap} overlaps, ratio = {ratio*100:.2f}%")

    valid_candidates = [col for col, (overlap, _) in candidate_results.items() if overlap > 0]
    chosen_column = valid_candidates[0] if len(valid_candidates) == 1 else max(candidate_results, key=lambda k: candidate_results[k][1])

    ref_first_value = str(ref_df[chosen_column].iloc[0]).strip()
    gene_identifier_type = "Ensembl_ID" if re.match(r'^ENSG', ref_first_value) else "Gene_symbol"
    chosen_overlap_ratio = candidate_results[chosen_column][1]
    print("Chosen bulk gene identifier type:", gene_identifier_type, f"with ratio = {chosen_overlap_ratio*100:.2f}%")
    return gene_identifier_type, chosen_overlap_ratio

def compute_overlap_ratio(bulk_expression_file: str, ref_expression_file: str, gene_type: str) -> float:
    """
    Compute the overlap ratio between bulk RNA gene identifiers and a candidate reference gene identifier column,
    based on the provided gene type.
    """
    bulk_df = pd.read_csv(bulk_expression_file, usecols=[0])
    bulk_genes = bulk_df.iloc[:, 0].astype(str).apply(
        lambda x: strip_ensembl_version(x.strip()) if gene_type == "Ensembl_ID" else x.strip()
    )
    bulk_gene_set = set(bulk_genes)

    ref_df = pd.read_csv(ref_expression_file, usecols=[0, 1])
    chosen_column = None
    for col in ref_df.columns:
        first_value = str(ref_df[col].iloc[0]).strip()
        if gene_type == "Ensembl_ID" and re.match(r'^ENSG', first_value):
            chosen_column = col
            break
        elif gene_type == "Gene_symbol" and not re.match(r'^ENSG', first_value):
            chosen_column = col
            break
    if chosen_column is None:
        return 0.0
    ref_genes = ref_df[chosen_column].astype(str).apply(lambda x: strip_ensembl_version(x.strip()) if gene_type=="Ensembl_ID" else x.strip())
    ref_gene_set = set(ref_genes)
    overlap = len(bulk_gene_set & ref_gene_set)
    overlap_ratio = overlap / len(bulk_gene_set) if bulk_gene_set else 0
    return overlap_ratio

def process_reference_expression_matrix(ref_expression_file: str, bulk_gene_type: str) -> str:
    """
    Process the reference expression matrix by retaining one gene identifier column and numeric expression columns.
    Returns the file path to the processed reference expression matrix.
    """
    ref_df = pd.read_csv(ref_expression_file)

    if ref_df.shape[1] < 2:
        print("Error: Reference file must have at least two columns.")
        sys.exit(1)

    col0, col1 = ref_df.columns[0], ref_df.columns[1]
    candidate_col0 = is_gene_identifier_column(ref_df[col0])
    candidate_col1 = is_gene_identifier_column(ref_df[col1])

    if candidate_col0 and candidate_col1:
        if bulk_gene_type == "Ensembl_ID":
            chosen_ref_column = col0 if re.match(r'^ENSG', str(ref_df[col0].iloc[0]).strip()) else col1
        else:
            chosen_ref_column = col0 if not re.match(r'^ENSG', str(ref_df[col0].iloc[0]).strip()) else col1
    elif candidate_col0:
        chosen_ref_column = col0
    elif candidate_col1:
        chosen_ref_column = col1
    else:
        print("Error: No valid gene identifier column found in the reference file.")
        sys.exit(1)

    # Retain the chosen gene identifier column and numeric columns.
    selected_columns = [chosen_ref_column]
    for col in ref_df.columns:
        if col == chosen_ref_column:
            continue
        numeric_series = pd.to_numeric(ref_df[col], errors='coerce')
        if numeric_series.notnull().mean() >= 0.5:
            selected_columns.append(col)
    processed_df = ref_df[selected_columns].copy()

    # Convert non-gene columns to numeric.
    for col in processed_df.columns[1:]:
        processed_df[col] = pd.to_numeric(processed_df[col], errors='coerce')

    if bulk_gene_type == "Ensembl_ID":
        processed_df[chosen_ref_column] = processed_df[chosen_ref_column].astype(str).apply(lambda x: strip_ensembl_version(x.strip()))

    processed_ref_expression_file = os.path.join(
        os.path.dirname(ref_expression_file), "processed_" + os.path.basename(ref_expression_file)
    )
    processed_df.to_csv(processed_ref_expression_file, index=False)
    return processed_ref_expression_file

def plot_deconvolution_results(deconv_csv_file: str, output_prefix: str, bar_height: float) -> None:
    """
    Generate and save a horizontal bar plot of the deconvolution results.
    The x-axis displays numbers starting from 0 with an interval of 5.
    """
    df = pd.read_csv(deconv_csv_file)

    if len(df.columns) < 2:
        print("Error: Deconvolution CSV file must have at least two columns.")
        return

    cell_type_column, fraction_column = df.columns[0], df.columns[1]
    # Convert fraction to percentage values for plotting
    df[fraction_column] = df[fraction_column] * 100
    df.sort_values(by=fraction_column, ascending=False, inplace=True)

    sns.set_style("whitegrid")
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 6

    num_cell_types = df.shape[0]
    figure_height = num_cell_types * 0.12
    plt.figure(figsize=(3, figure_height))

    ax = sns.barplot(x=fraction_column, y=cell_type_column, data=df, palette=sns.color_palette("bright"))
    # Set each bar's height
    for patch in ax.patches:
        patch.set_height(bar_height)

    ax.set_xlabel("Percentage contribution%", fontsize=6, fontname="Arial")
    ax.set_ylabel("Cell/tissue types", fontsize=6, fontname="Arial")
    ax.set_title("RNA tissue/cellular source", fontsize=6, fontname="Arial")
    
    # Set x-axis ticks starting at 0 with interval 5 and format labels to 2 decimals without '%'
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    plt.tight_layout()

    plt.savefig(f"{output_prefix}.png", dpi=600)
    plt.savefig(f"{output_prefix}.pdf", dpi=600)
    plt.savefig(f"{output_prefix}.svg", dpi=600)
    plt.show()

def main():
    parser = argparse.ArgumentParser(
        description="Perform RNA deconvolution using ARIC and visualize the results.\n\n"
                    "Input files:\n"
                    "  --bulk_expr_file: CSV with exactly 2 columns: gene identifier and expression for bulk RNA.\n"
                    "  --ref_expr_file: CSV with one or two candidate gene identifier columns (non-numeric) followed by numeric expression values."
    )
    parser.add_argument("--bulk_expr_file", required=True,
                        help="CSV file with 2 columns: gene identifier and expression for bulk RNA.")
    parser.add_argument("--ref_expr_file", required=True,
                        help="CSV file with one or two candidate gene identifier columns (non-numeric) followed by numeric expression values.")
    parser.add_argument("--bulk_gene_type", choices=["Gene_symbol", "Ensembl_ID"],
                        default=None,
                        help="Optional: specify gene identifier type for bulk RNA. If not provided, it will be determined automatically.")
    parser.add_argument("--bar_height", type=float, default=1.0,
                        help="Optional: set the bar height for the bar plot (default: 1.0).")
    args = parser.parse_args()

    # Determine or compute the gene identifier type and overlap ratio
    if args.bulk_gene_type is None:
        gene_identifier_type, chosen_overlap_ratio = determine_bulk_gene_identifier_type(args.bulk_expr_file, args.ref_expr_file)
    else:
        gene_identifier_type = args.bulk_gene_type
        chosen_overlap_ratio = compute_overlap_ratio(args.bulk_expr_file, args.ref_expr_file, gene_identifier_type)
        print("User provided bulk gene identifier type:", gene_identifier_type,
              f"with computed overlap ratio = {chosen_overlap_ratio*100:.2f}%")

    # Prepare the reference expression file
    processed_ref_expr_file = process_reference_expression_matrix(args.ref_expr_file, gene_identifier_type)

    # Build the output prefix
    bulk_base_name = os.path.splitext(os.path.basename(args.bulk_expr_file))[0]
    ref_base_name = os.path.splitext(os.path.basename(args.ref_expr_file))[0]
    ratio_str = f"{chosen_overlap_ratio*100:.2f}%" if chosen_overlap_ratio is not None else "NA"

    if gene_identifier_type == "Ensembl_ID":
        gene_suffix = "ensembl_IDs"
    elif gene_identifier_type == "Gene_symbol":
        gene_suffix = "gene_symbols"
    else:
        gene_suffix = "unknown"

    deconv_output_prefix = f"{bulk_base_name}_{ref_base_name}_deconvolution_fraction_{ratio_str}_overlaped_{gene_suffix}"
    print("Deconvolution output prefix:", deconv_output_prefix)

    # Check if the final deconvolution CSV already exists.
    deconv_csv_file = f"{deconv_output_prefix}_deconvolution.csv"
    if os.path.exists(deconv_csv_file):
        print(f"Found {deconv_csv_file}, we start plotting.")
        plot_deconvolution_results(deconv_csv_file, deconv_output_prefix, args.bar_height)
        # Clean up the temporary processed reference file
        if os.path.exists(processed_ref_expr_file):
            os.remove(processed_ref_expr_file)
        sys.exit(0)

    # Run ARIC deconvolution if CSV does not exist.
    ARIC(
        args.bulk_expr_file,
        processed_ref_expr_file,
        save_path=deconv_csv_file,
        marker_path=None,
        selected_marker=False,
        scale=0.1,
        delcol_factor=10,
        iter_num=10,
        confidence=0.75,
        w_thresh=10,
        unknown=False,
        is_methylation=False
    )

    # Check if ARIC produced the CSV
    if not os.path.exists(deconv_csv_file):
        print(f"Error: '{deconv_csv_file}' not found.")
        sys.exit(1)

    # Plot the results
    plot_deconvolution_results(deconv_csv_file, deconv_output_prefix, args.bar_height)

    # Remove the temporary processed reference file
    if os.path.exists(processed_ref_expr_file):
        os.remove(processed_ref_expr_file)

if __name__ == "__main__":
    main()

