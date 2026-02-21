#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Example usage:
#     python get_umi_seq_motif_for_individual_fq.py -head SAMPLE_NAME -fq /path/to/read2.fastq.gz -n 14 -r 1000000 -o /path/to/motif_output
#     -n: Number of bases to extract from each read.
#     -r: Maximum number of reads to process.
#     -o: Output directory for CSV and motif logo images.


import os
import sys
import gzip
import argparse
import logging
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import logomaker
import subprocess
import tempfile
import numpy as np

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42

def setup_logging(output_dir='motif_output'):
    """Set up logging to console and file."""
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, 'processing.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_path)
        ]
    )

def parse_fastq_and_count(file_path, first_n_bases=14, max_reads=1000000):
    """
    Parse FASTQ file, extract the first N bases from each read,
    and count the frequency of bases at each position.
    """
    counts = [defaultdict(int) for _ in range(first_n_bases)]
    open_func = gzip.open if file_path.endswith('.gz') else open
    reads_processed = 0
    try:
        with open_func(file_path, 'rt') as f:
            while reads_processed < max_reads:
                try:
                    _header = next(f).strip()
                    seq = next(f).strip()[:first_n_bases]
                    _plus = next(f).strip()
                    _qual = next(f).strip()
                    if len(seq) == first_n_bases:
                        for i, base in enumerate(seq):
                            counts[i][base.upper()] += 1
                    reads_processed += 1
                except StopIteration:
                    break
    except Exception as e:
        logging.error(f"Error processing {file_path}: {e}")
        return None

    data = []
    for i, counter in enumerate(counts):
        total = sum(counter.values())
        freq = {b: (counter[b] / total if total else 0) for b in counter}
        freq['Position'] = i + 1
        data.append(freq)

    df = pd.DataFrame(data).fillna(0).set_index('Position').sort_index()
    # Ensure columns A, C, G, T exist
    for b in ['A', 'C', 'G', 'T']:
        if b not in df.columns:
            df[b] = 0
    return df[['A', 'C', 'G', 'T']]

def main():
    parser = argparse.ArgumentParser(
        description="Analyze the first N bases of a FASTQ file (Read2) and generate a motif logo.\n"
                    "Parameters:\n"
                    "  -head: Sample name header used for output file naming.\n"
                    "  -fq:   Path to the Read2 FASTQ file (supports .gz).\n"
                    "  -n:    Number of bases to extract from each read (default: 14).\n"
                    "  -r:    Maximum number of reads to process (default: 1000000).\n"
                    "  -o:    Output directory for CSV and motif logo images (default: motif_output)."
    )
    parser.add_argument('-head', required=True, help='Sample name header for output files.')
    parser.add_argument('-fq', required=True, help='Path to the Read2 FASTQ file.')
    parser.add_argument('-n', '--num_bases', type=int, default=14, help='Number of bases to extract from each read.')
    parser.add_argument('-r', '--num_reads', type=int, default=1000000, help='Maximum number of reads to process.')
    parser.add_argument('-o', '--output_dir_qc', default='motif_output', help='Output directory for CSV and motif logo images.')
    args = parser.parse_args()

    setup_logging(args.output_dir_qc)
    logging.info("Starting analysis...")

    os.makedirs(args.output_dir_qc, exist_ok=True)

    # Sample reads randomly using seqtk
    logging.info("Sampling reads with seqtk...")
    temp_fastq_path = None
    try:
        # Create a temporary file for sampled reads
        temp_fastq = tempfile.NamedTemporaryFile(delete=False, suffix=".fastq")
        temp_fastq_path = temp_fastq.name
        temp_fastq.close()
        cmd = ["seqtk", "sample", "-s", "42", args.fq, str(args.num_reads)]
        with open(temp_fastq_path, "w") as out_f:
            subprocess.run(cmd, stdout=out_f, check=True)
    except Exception as e:
        logging.error(f"Error sampling reads with seqtk: {e}")
        if temp_fastq_path and os.path.exists(temp_fastq_path):
            os.remove(temp_fastq_path)
        sys.exit(1)

    # Parse sampled FASTQ and count base frequencies
    try:
        df = parse_fastq_and_count(
            file_path=temp_fastq_path,
            first_n_bases=args.num_bases,
            max_reads=args.num_reads
        )
    finally:
        # Ensure temporary file is always cleaned up
        if temp_fastq_path and os.path.exists(temp_fastq_path):
            os.remove(temp_fastq_path)

    if df is None:
        logging.error("No data to plot. Exiting.")
        sys.exit(1)

    csv_path = os.path.join(args.output_dir_qc, f"{args.head}_base_distribution.csv")
    df.to_csv(csv_path)
    logging.info(f"Base distribution saved to {csv_path}")

    fig, ax = plt.subplots(figsize=(6, 2))
    try:
        logo = logomaker.Logo(df, ax=ax, shade_below=0.5, fade_below=0.5)
        logo.style_spines(visible=False)
        logo.style_spines(spines=['left', 'bottom'], visible=True)
        logo.style_xticks(anchor=0, spacing=1)
        ax.set_xlabel('Position')
        ax.set_ylabel('Frequency')

        combined_title = f"Frequency distribution of Read2 UMI linker region\n{args.head}"
        ax.set_title(combined_title, fontsize=10)

        # Set specific Y-axis ticks
        ax.set_yticks(np.arange(0, 1.01, 0.25))
        plt.tight_layout()

        # Save motif logo in three formats: PNG and PDF
        png_path = os.path.join(args.output_dir_qc, f"{args.head}_{args.num_bases}BP_motif_logo.png")
        pdf_path = os.path.join(args.output_dir_qc, f"{args.head}_{args.num_bases}BP_motif_logo.pdf")
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path, dpi=300)
        plt.close()
        logging.info(f"Motif logo saved as: {png_path}, {pdf_path}")
    except Exception as e:
        logging.error(f"Error plotting motif: {e}")
        plt.close()

if __name__ == "__main__":
    main()
