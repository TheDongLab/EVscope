#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Convert RSEM output to expression matrix with gene annotations"
    )
    parser.add_argument(
        "--RSEM_out",
        required=True,
        help="Path to the RSEM genes results file (e.g., Ewa_DNAnexus_CTRL_SNCA_10000_RSEM.genes.results)"
    )
    parser.add_argument(
        "--GeneID_meta_table",
        required=True,
        help="Path to the gene annotation meta table (TSV format)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output expression matrix file"
    )
    args = parser.parse_args()

    # Read the gene annotation meta table using tab delimiter
    meta_df = pd.read_csv(args.GeneID_meta_table, sep="\t")
    # If the expected column is missing, try using whitespace delimiter
    if "GeneID" not in meta_df.columns:
        meta_df = pd.read_csv(args.GeneID_meta_table, sep='\s+')

    # Read the RSEM results file (tab-delimited)
    rsem_df = pd.read_csv(args.RSEM_out, sep="\t")

    # Merge meta table and RSEM results using GeneID from meta table and gene_id from RSEM file
    merged_df = pd.merge(meta_df, rsem_df, left_on="GeneID", right_on="gene_id", how="inner")

    # Select and rename columns:
    # - GeneID, GeneSymbol, GeneType from meta table
    # - ReadCounts (from expected_count) and TPM from RSEM results
    output_df = merged_df[["GeneID", "GeneSymbol", "GeneType", "expected_count", "TPM"]].copy()
    output_df.rename(columns={"expected_count": "ReadCounts"}, inplace=True)

    # Write the output expression matrix to a TSV file
    output_df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()

