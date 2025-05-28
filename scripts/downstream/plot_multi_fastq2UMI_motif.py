#!/usr/bin/env python
import sys
import os
import argparse
import logging
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import math
import subprocess
import matplotlib as mpl  # for SVG settings

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('processing.log')
        ]
    )

def parse_fastq_and_count(file_path, first_n_bases=14, num_reads=1000000, sample_name=None):
    """
    Use seqtk to randomly sample FASTQ file and count base frequencies.
    
    Args:
        file_path (str): Path to FASTQ file (can be gzipped).
        first_n_bases (int): Number of bases to extract.
        num_reads (int): Number of reads to sample.
        sample_name (str): Optional sample name for visualization.
    
    Returns:
        tuple: (sample_name, frequency dataframe)
    """
    try:
        # Run seqtk sample command to randomly select reads
        result = subprocess.run(["seqtk", "sample", file_path, str(num_reads)],
                                capture_output=True, text=True, check=True)
        fastq_data = result.stdout.splitlines()
    except Exception as e:
        logging.error(f"Error running seqtk on {file_path}: {e}")
        return None

    # Initialize count list for each base position
    counts = [defaultdict(int) for _ in range(first_n_bases)]
    i = 0
    # Process FASTQ data in chunks of 4 lines
    while i < len(fastq_data):
        if i + 3 >= len(fastq_data):
            break
        header = fastq_data[i].strip()
        seq_line = fastq_data[i+1].strip()
        plus = fastq_data[i+2].strip()
        qual = fastq_data[i+3].strip()
        i += 4

        seq = seq_line[:first_n_bases]
        if len(seq) != first_n_bases:
            continue
        for j, base in enumerate(seq):
            counts[j][base.upper()] += 1

    # Build a dataframe with frequency for A, C, G, T at each position
    freq_data = []
    for pos, counter in enumerate(counts):
        total = sum(counter.values())
        freq = {base: (counter.get(base, 0) / total if total > 0 else 0) for base in ['A', 'C', 'G', 'T']}
        freq['Position'] = pos + 1
        freq_data.append(freq)

    df = pd.DataFrame(freq_data).set_index('Position')
    if sample_name is None:
        sample_name = os.path.basename(file_path).split('.')[0]
    return sample_name, df

def process_file(sample_name, file_path, output_dir, first_n_bases, num_reads):
    """
    Process one FASTQ file: sample reads, count base frequencies, and save CSV.
    """
    result = parse_fastq_and_count(file_path, first_n_bases, num_reads, sample_name=sample_name)
    if result is None:
        logging.warning(f"Skipping {file_path} due to errors.")
        return None
    sample_name, df = result
    logging.info(f"{sample_name}: Processed {num_reads} reads using seqtk.")
    csv_name = os.path.join(output_dir, f"{sample_name}_base_distribution.csv")
    try:
        df.to_csv(csv_name)
        logging.info(f"{sample_name}: Base distribution saved to {csv_name}")
    except Exception as e:
        logging.error(f"Error saving CSV for {sample_name}: {e}")
        return None
    return sample_name, df

def calculate_layout(num_samples, max_cols=5):
    """
    Calculate subplot layout based on number of samples.
    
    Args:
        num_samples (int): Number of samples.
        max_cols (int): Maximum number of columns.
    
    Returns:
        tuple: (nrows, ncols)
    """
    ncols = min(max_cols, num_samples)
    nrows = math.ceil(num_samples / ncols)
    return nrows, ncols

def main():
    parser = argparse.ArgumentParser(
        description="Extract random NUM_READS from FASTQ files using seqtk, analyze base distribution, and generate motif logos."
    )
    parser.add_argument('--fastq_files', required=True,
                        help='Tab-separated file with sample short name and FASTQ file path.')
    parser.add_argument('-n', '--num_bases', type=int, default=14,
                        help='Number of bases to extract (default: 14).')
    parser.add_argument('-r', '--num_reads', type=int, default=1000000,
                        help='Number of reads to extract from each FASTQ file (default: 1000000).')
    parser.add_argument('-o', '--output_dir', default='motif_output',
                        help='OUTPUT_DIR to save CSV files and motif logos (default: motif_output).')
    args = parser.parse_args()

    setup_logging()
    logging.info("Starting processing...")

    # Read TSV file with sample names and FASTQ file paths
    fastq_info = []
    try:
        with open(args.fastq_files, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    logging.warning(f"Skipping line due to insufficient columns: {line}")
                    continue
                sample_name = parts[0]
                file_path = parts[1]
                fastq_info.append((sample_name, file_path))
    except Exception as e:
        logging.error(f"Error reading FASTQ TSV file: {e}")
        sys.exit(1)

    if not fastq_info:
        logging.error("No FASTQ file information found.")
        sys.exit(1)

    first_n_bases = args.num_bases
    num_reads = args.num_reads
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    results = []
    for sample_name, file_path in fastq_info:
        res = process_file(sample_name, file_path, output_dir, first_n_bases, num_reads)
        if res is not None:
            results.append(res)

    if not results:
        logging.error("No samples were successfully processed.")
        sys.exit(1)

    # Plot motif logos for all samples
    num_samples = len(results)
    nrows, ncols = calculate_layout(num_samples, max_cols=5)
    width_per_subplot = 4
    height_per_subplot = 2
    fig_width = width_per_subplot * ncols
    fig_height = height_per_subplot * nrows
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width, fig_height))
    if num_samples == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for i, (sample_name, df) in enumerate(results):
        ax = axes[i]
        try:
            logo = logomaker.Logo(df, ax=ax, shade_below=0.5, fade_below=0.5, font_name='Arial')
            logo.style_spines(visible=False)
            logo.style_spines(spines=['left', 'bottom'], visible=True)
            logo.style_xticks(anchor=0, spacing=1)
            ax.set_xlabel('Position')
            ax.set_ylabel('Frequency')
            ax.set_title(sample_name)
        except Exception as e:
            logging.error(f"Error plotting logo for {sample_name}: {e}")

    # Remove unused subplots
    for j in range(len(results), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()

    # Save PNG image
    output_png = os.path.join(output_dir, f"{num_samples}_{first_n_bases}BP_motif_logos.png")
    try:
        plt.savefig(output_png, dpi=300)
        logging.info(f"Motif logos saved as {output_png}")
    except Exception as e:
        logging.error(f"Error saving motif logos PNG: {e}")

    # Save SVG image with editable vectorized text
    mpl.rcParams['svg.fonttype'] = 'none'
    output_svg = os.path.join(output_dir, f"{num_samples}_{first_n_bases}BP_motif_logos.svg")
    try:
        plt.savefig(output_svg, format='svg')
        logging.info(f"Motif logos saved as {output_svg}")
    except Exception as e:
        logging.error(f"Error saving motif logos SVG: {e}")
    plt.close()

if __name__ == "__main__":
    main()

