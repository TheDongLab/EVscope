#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Correct featureCounts for gDNA contamination in strand-specific RNA-seq data.")
    parser.add_argument('--strand', required=True, choices=['forward', 'reverse'],
                        help="Strand specificity of RNA-seq library (forward or reverse).")
    parser.add_argument('--forward_featureCounts_table', required=True,
                        help="featureCounts output table from forward strand.")
    parser.add_argument('--reverse_featureCounts_table', required=True,
                        help="featureCounts output table from reverse strand.")
    parser.add_argument('--output', required=True,
                        help="Output file for gDNA-corrected read counts.")
    return parser.parse_args()

def read_featureCounts(file_path):
    """Read featureCounts file."""
    df = pd.read_csv(file_path, sep='\t', comment='#', header=0)
    return df

def correct_gdna_contamination(forward_df, reverse_df, strand):
    """Calculate gDNA-corrected counts based on strand specificity."""
    count_column = forward_df.columns[-1]

    merged_df = forward_df.copy()

    if strand == 'reverse':
        merged_df[count_column] = reverse_df[count_column] - forward_df[count_column]
    else:  # strand == 'forward'
        merged_df[count_column] = forward_df[count_column] - reverse_df[count_column]

    # Replace negative values with zero
    merged_df[count_column] = merged_df[count_column].clip(lower=0)

    return merged_df

def main():
    args = parse_arguments()

    # Read input files
    forward_counts_df = read_featureCounts(args.forward_featureCounts_table)
    reverse_counts_df = read_featureCounts(args.reverse_featureCounts_table)

    # Perform correction
    corrected_counts_df = correct_gdna_contamination(forward_counts_df, reverse_counts_df, args.strand)

    # Save results
    corrected_counts_df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()

