#!/usr/bin/env python3

import argparse
import gzip
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

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
        with opener(file, 'rt', errors='ignore') as fq:
            for i, line in enumerate(fq):
                if i % 4 == 1:  # This is the sequence line in a FASTQ record
                    lengths.append(len(line.strip()))
    except FileNotFoundError:
        raise ValueError(f"File not found: {file}")
    except Exception as e:
        raise ValueError(f"Error reading {file}: {str(e)}")
    return lengths

def plot_length_distribution(lengths_list, output_pdf, output_png, titles):
    """
    Generate a vertically stacked read length distribution plot for three FASTQ files.
    """
    n_files = len(lengths_list)
    if n_files != 3:
        raise ValueError("This script is designed for exactly 3 input files.")

    # Create a figure with 3 vertically stacked subplots that share the same X-axis
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)

    # Determine the global maximum read length for a uniform x-axis across all plots
    # This ensures the scales are comparable.
    if not any(lengths_list):
        raise ValueError("All input files are empty.")
    max_length = max(max(lengths) for lengths in lengths_list if lengths)

    # Define x-axis ticks for clarity
    tick_step = 20 if max_length > 100 else 10
    x_ticks = np.arange(0, max_length + tick_step, tick_step)

    for i, (ax, lengths, title) in enumerate(zip(axes, lengths_list, titles)):
        if not lengths:
            print(f"Warning: No reads found in file for plot '{title}'. Skipping.")
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
            ax.set_title(title, fontsize=14, weight='bold')
            continue

        counter = Counter(lengths)
        sorted_lengths = sorted(counter.items())
        x, y = zip(*sorted_lengths)
        
        # Use log10 for frequency to handle large variations
        y_log = np.log10(np.array(y, dtype=float))

        ax.bar(x, y_log, color='cornflowerblue', edgecolor='black', linewidth=0.5)
        
        # Set individual subplot titles
        ax.set_title(title, fontsize=14, weight='bold', pad=10)
        
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        ax.tick_params(axis='both', which='major', labelsize=10, direction='in')
        ax.set_xticks(x_ticks)
        ax.set_xlim(0, max_length + 5)

    # Add a single, shared X-axis label to the bottom plot
    axes[-1].set_xlabel('Read Length (bp)', fontsize=12)

    # Add a single, centered Y-axis label for the entire figure
    fig.text(0.04, 0.5, 'Log10(Frequency)', va='center', rotation='vertical', fontsize=12)

    # Adjust layout to minimize whitespace and prevent titles/labels from overlapping
    plt.tight_layout(rect=[0.05, 0, 1, 1]) # rect=[left, bottom, right, top]
    fig.subplots_adjust(hspace=0.35) # Add a small amount of vertical space between plots

    # Save the figure in high-quality formats suitable for publication and web
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plots saved to {output_pdf} and {output_png}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate a publication-quality read length distribution plot from three FASTQ files.'
    )
    parser.add_argument(
        '--input_fastqs', 
        nargs=3, 
        required=True, 
        metavar='FILE',
        help='Three input FASTQ or FASTQ.GZ files (e.g., R1_trimmed.fq.gz R1_clean.fq.gz R2_clean.fq.gz).'
    )
    parser.add_argument(
        '--titles',
        nargs=3,
        required=True,
        metavar='TITLE',
        help='Three titles for the plots, in the same order as the input files. (e.g., "R1 Adapter-Trimmed" "R1 Clean Reads" "R2 Clean Reads")'
    )
    parser.add_argument('--output_pdf', required=True, help='Output PDF file path.')
    parser.add_argument('--output_png', required=True, help='Output PNG file path.')
    args = parser.parse_args()

    # Process each file to get read lengths
    lengths_list = []
    for fastq_file in args.input_fastqs:
        print(f"Processing {os.path.basename(fastq_file)}...")
        lengths = read_lengths_from_fastq(fastq_file)
        lengths_list.append(lengths)

    # Generate the plot
    plot_length_distribution(lengths_list, args.output_pdf, args.output_png, args.titles)

