#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RNA deconvolution using ARIC with optional output directory.
"""
import argparse
import os
import re
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from ARIC import ARIC  # ARIC deconvolution function

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42


def strip_ensembl_version(gene_id: str) -> str:
    """Remove version suffix from Ensembl ID."""
    return re.sub(r'\.\d+$', '', gene_id)


def is_gene_identifier_column(series: pd.Series) -> bool:
    """Check if a Series is non-numeric gene IDs."""
    numeric_count = pd.to_numeric(series, errors='coerce').notnull().sum()
    return numeric_count < len(series)


def determine_bulk_gene_identifier_type(bulk_file: str, ref_file: str) -> (str, float):
    """Auto-detect bulk gene ID type by overlap with reference."""
    bulk_df = pd.read_csv(bulk_file, usecols=[0])
    bulk_genes = bulk_df.iloc[:, 0].astype(str).apply(
        lambda x: strip_ensembl_version(x.strip()) if re.match(r'^ENSG\d+', x.strip()) else x.strip()
    )
    bulk_set = set(bulk_genes)

    ref_df = pd.read_csv(ref_file, usecols=[0,1])
    candidates = {}
    for col in ref_df.columns:
        col_data = ref_df[col].astype(str).str.strip()
        if is_gene_identifier_column(col_data):
            if re.match(r'^ENSG', col_data.iloc[0]):
                genes = col_data.apply(strip_ensembl_version)
            else:
                genes = col_data
            ref_set = set(genes)
            ov = len(bulk_set & ref_set)
            ratio = ov / len(bulk_set) if bulk_set else 0
            candidates[col] = (ov, ratio)
        else:
            candidates[col] = (0, 0)

    for col, (ov, ratio) in candidates.items():
        print(f"Ref column '{col}': {ov} overlaps ({ratio*100:.2f}%)")

    # choose best match
    valid = [c for c, (ov,_) in candidates.items() if ov > 0]
    chosen = valid[0] if len(valid)==1 else max(candidates, key=lambda k: candidates[k][1])
    first_val = str(ref_df[chosen].iloc[0]).strip()
    gene_type = "Ensembl_ID" if re.match(r'^ENSG', first_val) else "Gene_symbol"
    print(f"Chosen bulk gene id type: {gene_type} ({candidates[chosen][1]*100:.2f}%)")
    return gene_type, candidates[chosen][1]


def compute_overlap_ratio(bulk_file: str, ref_file: str, gene_type: str) -> float:
    """Compute overlap ratio for given gene type."""
    bulk_df = pd.read_csv(bulk_file, usecols=[0])
    bulk_genes = bulk_df.iloc[:,0].astype(str).apply(
        lambda x: strip_ensembl_version(x.strip()) if gene_type=='Ensembl_ID' else x.strip()
    )
    bulk_set = set(bulk_genes)

    ref_df = pd.read_csv(ref_file, usecols=[0,1])
    col = None
    for c in ref_df.columns:
        val = str(ref_df[c].iloc[0]).strip()
        if gene_type=='Ensembl_ID' and re.match(r'^ENSG', val):
            col = c; break
        if gene_type=='Gene_symbol' and not re.match(r'^ENSG', val):
            col = c; break
    if col is None:
        return 0.0
    ref_set = set(ref_df[col].astype(str).str.strip().apply(
        lambda x: strip_ensembl_version(x) if gene_type=='Ensembl_ID' else x
    ))
    ov = len(bulk_set & ref_set)
    return ov / len(bulk_set) if bulk_set else 0


def process_reference_expression_matrix(ref_file: str, gene_type: str, output_dir: str) -> str:
    """Filter ref matrix and save to output_dir."""
    df = pd.read_csv(ref_file)
    if df.shape[1] < 2:
        print("Error: Reference file needs at least two columns.")
        sys.exit(1)

    c0, c1 = df.columns[0], df.columns[1]
    id0, id1 = is_gene_identifier_column(df[c0]), is_gene_identifier_column(df[c1])
    if id0 and id1:
        if gene_type=='Ensembl_ID':
            chosen = c0 if re.match(r'^ENSG', str(df[c0].iloc[0]).strip()) else c1
        else:
            chosen = c0 if not re.match(r'^ENSG', str(df[c0].iloc[0]).strip()) else c1
    elif id0:
        chosen = c0
    elif id1:
        chosen = c1
    else:
        print("Error: No gene ID column in reference.")
        sys.exit(1)

    select_cols = [chosen]
    for col in df.columns:
        if col==chosen: continue
        ser = pd.to_numeric(df[col], errors='coerce')
        if ser.notnull().mean() >= 0.5:
            select_cols.append(col)
    proc = df[select_cols].copy()
    for col in proc.columns[1:]:
        proc[col] = pd.to_numeric(proc[col], errors='coerce')
    if gene_type=='Ensembl_ID':
        proc[chosen] = proc[chosen].astype(str).apply(strip_ensembl_version)

    out_file = os.path.join(output_dir, "processed_" + os.path.basename(ref_file))
    proc.to_csv(out_file, index=False)
    return out_file


def plot_deconvolution_results(deconv_csv: str, output_prefix: str, bar_height: float) -> None:
    """Plot and save deconvolution bar chart."""
    df = pd.read_csv(deconv_csv)
    ct_col, frac_col = df.columns[0], df.columns[1]
    df[frac_col] = df[frac_col] * 100
    df.sort_values(by=frac_col, ascending=False, inplace=True)

    sns.set_style("whitegrid")
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 6

    height = df.shape[0] * 0.12
    plt.figure(figsize=(3, height))
    ax = sns.barplot(x=frac_col, y=ct_col, data=df, palette=sns.color_palette("bright"))
    for p in ax.patches:
        p.set_height(bar_height)
    ax.set_xlabel("Percentage contribution%", fontsize=6)
    ax.set_ylabel("Cell/tissue types", fontsize=6)
    ax.set_title("RNA tissue/cellular source", fontsize=6)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", dpi=600)
    plt.savefig(f"{output_prefix}.pdf", dpi=600)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Perform RNA deconvolution using ARIC.")
    parser.add_argument("--bulk_expr_file", required=True,
                        help="CSV with 2 columns: gene ID and expression.")
    parser.add_argument("--ref_expr_file", required=True,
                        help="CSV with gene ID columns and expression values.")
    parser.add_argument("--bulk_gene_type", choices=["Gene_symbol","Ensembl_ID"],
                        default=None,
                        help="Specify bulk gene ID type." )
    parser.add_argument("--bar_height", type=float, default=1.0,
                        help="Bar height in plot.")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to save all outputs.")
    args = parser.parse_args()

    # ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # determine gene ID type
    if args.bulk_gene_type is None:
        gene_type, overlap = determine_bulk_gene_identifier_type(
            args.bulk_expr_file, args.ref_expr_file)
    else:
        gene_type = args.bulk_gene_type
        overlap = compute_overlap_ratio(
            args.bulk_expr_file, args.ref_expr_file, gene_type)
        print(f"Using provided gene type: {gene_type} ({overlap*100:.2f}%)")

    # process reference file into output_dir
    processed_ref = process_reference_expression_matrix(
        args.ref_expr_file, gene_type, args.output_dir)

    # build output prefix in output_dir
    bulk_base = os.path.splitext(os.path.basename(args.bulk_expr_file))[0]
    ref_base = os.path.splitext(os.path.basename(args.ref_expr_file))[0]
    ratio_str = f"{overlap*100:.2f}%"
    suffix = "ensembl_IDs" if gene_type=='Ensembl_ID' else "gene_symbols"
    prefix = os.path.join(
        args.output_dir,
        f"{bulk_base}_{ref_base}_deconvolution_fraction_{ratio_str}_overlaped_{suffix}"
    )
    print(f"Output prefix: {prefix}")

    csv_file = f"{prefix}_deconvolution.csv"
    if os.path.exists(csv_file):
        print("Found existing deconvolution CSV, plotting...")
        plot_deconvolution_results(csv_file, prefix, args.bar_height)
        os.remove(processed_ref)
        sys.exit(0)

    # run ARIC
    ARIC(
        args.bulk_expr_file,
        processed_ref,
        save_path=csv_file,
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

    if not os.path.exists(csv_file):
        print(f"Error: '{csv_file}' not found.")
        sys.exit(1)

    # plot and cleanup
    plot_deconvolution_results(csv_file, prefix, args.bar_height)
    os.remove(processed_ref)

if __name__ == "__main__":
    main()

