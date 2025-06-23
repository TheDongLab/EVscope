#!/usr/bin/env python3

import argparse
import gzip
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# Set publication-quality parameters
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 1.2
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.family'] = 'Arial'

def read_lengths_from_fastq(file):
    """Extract read lengths from a FASTQ or FASTQ.GZ file."""
    lengths = []
    opener = gzip.open if file.endswith('.gz') else open
    try:
        with opener(file, 'rt') as fq:
            for i, line in enumerate(fq):
                if i % 4 == 1:  # sequence line in fastq
                    lengths.append(len(line.strip()))
    except FileNotFoundError:
        raise ValueError(f"File not found: {file}")
    except Exception as e:
        raise ValueError(f"Error reading {file}: {str(e)}")
    return lengths

def plot_length_distribution(lengths_list, output_pdf, output_png, titles):
    """Generate read length distribution plot(s) with subplots for multiple files."""
    n_files = len(lengths_list)
    if n_files not in (1, 2):
        raise ValueError("Number of input files must be 1 or 2.")

    if n_files == 1:
        fig, ax = plt.subplots(figsize=(8, 6))
        axes = [ax]
    else:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
        axes = axes.flatten()

    # Determine global maximum read length for uniform x-axis
    max_length = max(max(lengths) for lengths in lengths_list) if lengths_list else 100
    x_ticks = np.arange(0, max_length + 10, 10)

    for ax, lengths, title in zip(axes, lengths_list, titles):
        if not lengths:
            raise ValueError(f"No reads found in {title}")
        counter = Counter(lengths)
        sorted_lengths = sorted(counter.items())
        x, y = zip(*sorted_lengths)
        x = np.array(x, dtype=int)
        y = np.log10(np.array(y, dtype=float))

        ax.bar(x, y, color='cornflowerblue', edgecolor='black', linewidth=0.8)
        ax.set_xlabel('Read Length (bp)', fontsize=12)
        if ax == axes[0]:
            ax.set_ylabel('Log10(Frequency)', fontsize=12)
        ax.set_title(title.split('.')[0], fontsize=14, weight='bold', pad=5)  # Use base filename without extension
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        ax.tick_params(axis='both', labelsize=10)
        ax.set_xticks(x_ticks)  # Uniform x-axis ticks
        ax.set_xlim(0, max_length + 10)  # Uniform x-axis range starting at 0

    # Set centered title above all subplots with one-line separation
    fig.suptitle('Read length distribution', fontsize=16, weight='bold', y=1.05)

    plt.tight_layout(pad=1.0)  # One-line separation
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate read length distribution from one or two FASTQ files (fastq or fastq.gz)')
    parser.add_argument('--input_fastqs', nargs='+', required=True, help='Input FASTQ or FASTQ.GZ file(s), up to 2 files (e.g., R1.fq.gz R2.fq.gz)')
    parser.add_argument('--output_pdf', required=True, help='Output PDF file for the distribution plot')
    parser.add_argument('--output_png', required=True, help='Output PNG file for the distribution plot')
    args = parser.parse_args()

    if len(args.input_fastqs) > 2:
        raise ValueError("Maximum of 2 FASTQ files are supported.")

    lengths_list = []
    titles = []
    for fastq_file in args.input_fastqs:
        lengths = read_lengths_from_fastq(fastq_file)
        lengths_list.append(lengths)
        titles.append(fastq_file.split('.')[0])  # Use base filename without extension

    plot_length_distribution(lengths_list, args.output_pdf, args.output_png, titles)
