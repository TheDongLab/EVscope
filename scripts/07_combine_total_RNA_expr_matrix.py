#!/usr/bin/env python3
"""Merge RNA and circRNA expression data with full columns"""
import argparse
import csv

def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene_expr", required=True)
    parser.add_argument("--circRNA_expr", required=True)
    parser.add_argument("--out_matrix", required=True)
    return parser.parse_args()

def read_gene_data(file_path):
    """Read gene expression data"""
    data = []
    with open(file_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row["GeneType"] == "artifact":
                continue
            try:
                data.append({
                    "GeneID": row["GeneID"],
                    "GeneSymbol": row["GeneSymbol"],
                    "GeneType": row["GeneType"],
                    "ReadCounts": row["ReadCounts"],
                    "Norm_Expr": row["TPM"]
                })
            except KeyError:
                continue
    return data

def read_circ_data(file_path):
    """Read circRNA expression data and return as gene-like rows."""
    data = []
    with open(file_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                circ_id = row["circRNA_ID1"]
                circ_symbol = row["circRNA_ID3"]
                circ_counts = row["junction_read_counts"]
                circ_cpm = row["CPM"]
            except KeyError:
                continue
            data.append({
                "GeneID": circ_id,
                "GeneSymbol": circ_symbol,
                "GeneType": "circRNAs",
                "ReadCounts": circ_counts,
                "Norm_Expr": circ_cpm
            })
    return data

def main():
    args = parse_args()

    # Read data
    gene_data = read_gene_data(args.gene_expr)
    circ_data = read_circ_data(args.circRNA_expr)

    # Write merged matrix
    with open(args.out_matrix, 'w') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["GeneID","GeneSymbol","GeneType","ReadCounts","Norm_Expr"],
            delimiter='\t'
        )
        writer.writeheader()
        writer.writerows(gene_data + circ_data)

if __name__ == "__main__":
    main()

