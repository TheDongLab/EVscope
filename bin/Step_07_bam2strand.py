#!/usr/bin/env python3
"""
Infer strand specificity from a BAM file, generate a TSV summary,
and produce a pie chart.
The default test read number is 100,000,000.
"""

import subprocess
import sys
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set font properties for publication quality
plt.rcParams["font.family"] = "Arial"
mpl.rcParams['pdf.fonttype'] = 42


def run_infer_experiment(bam_file, refgene_bed, test_read_num, output_dir):
    """
    Runs RSeQC's infer_experiment.py, parses the output, and generates
    a TSV file and a formatted pie chart.
    """
    # Ensure the input BAM file exists.
    if not os.path.isfile(bam_file):
        print(f"Error: Input BAM file not found at {bam_file}")
        sys.exit(1)

    # Create the output directory if it doesn't exist.
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.basename(bam_file).replace('.bam', '')

    # --- Run infer_experiment.py and write to TSV ---
    tsv_path = os.path.join(output_dir, f"{base_name}_bam2strandness.tsv")
    frac_fwd, frac_rev, frac_failed, lib_type = "0", "0", "0", "Unknown"

    with open(tsv_path, "w") as f:
        f.write("BAM\tdata type\t'1++,1--,2+-,2-+ (forward)'\t'1+-,1-+,2++,2-- (reverse)'\t%reads failed to determine:\n")
        try:
            # Construct and run the command.
            cmd = ["infer_experiment.py", "-i", bam_file, "-r", refgene_bed, "-s", str(test_read_num)]
            proc = subprocess.run(cmd, check=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # Parse the output from stdout.
            lines = proc.stdout.strip().split('\n')
            lib_type = lines[0].strip() if lines else "Library type not found"
            frac_failed = lines[1].split(": ")[1] if len(lines) > 1 else "0"
            frac_fwd = lines[2].split(": ")[1] if len(lines) > 2 else "0"
            frac_rev = lines[3].split(": ")[1] if len(lines) > 3 else "0"
            
            f.write(f"{bam_file}\t{lib_type}\t{frac_fwd}\t{frac_rev}\t{frac_failed}\n")

        except subprocess.CalledProcessError as e:
            error_message = f"Error processing {bam_file}: {e.stderr}"
            f.write(f"{bam_file}\tError\tNA\tNA\tNA\nDetails: {error_message}\n")
            print(error_message)
            sys.exit(1)

    print(f"Strandness results written to {tsv_path}")

    # --- Generate Publication-Quality Pie Chart ---
    try:
        f_failed = float(frac_failed)
        f_forward = float(frac_fwd)
        f_reverse = float(frac_rev)
    except ValueError:
        print("Error: Could not convert strand fractions to numbers. Skipping chart generation.")
        return

    # Chart data and styling
    labels = ['Forward strand', 'Reverse strand', 'Undetermined strand']
    sizes = [f_forward, f_reverse, f_failed]
    colors = ['#FF4500', '#FFA500', '#87CEEB'] # Red-Orange, Orange, Sky Blue

    # Create figure and axes. Adjust figsize for legend placement.
    fig, ax = plt.subplots(figsize=(5, 3), dpi=300)

    # Draw the pie chart.
    wedges, texts = ax.pie(
        sizes,
        startangle=90,
        colors=colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1},
        radius=1.2 # Use a standard radius for the pie.
    )

    ax.axis('equal')  # Ensures the pie is circular.

    # Set the three-line title with optimized spacing.
    title_str = f"Percentage of strand specificity\n{lib_type}\n{base_name}"
    ax.set_title(title_str, fontsize=8, pad=10)

    # Create legend labels with percentages.
    legend_labels = [f"{lbl} ({sz*100:.1f}%)" for lbl, sz in zip(labels, sizes)]
    
    # Add the legend to the right of the pie chart, positioned as close as possible.
    ax.legend(
        wedges,
        legend_labels,
        loc='center left', # Anchor point of the legend box.
        bbox_to_anchor=(0.80, 0.5), # Position the legend extremely close to the pie.
        frameon=False,
        prop={'size': 7}
    )

    # Save the chart in multiple formats with no extra whitespace.
    pie_pdf = os.path.join(output_dir, f"{base_name}_bam2strandness_pie.pdf")
    pie_png = os.path.join(output_dir, f"{base_name}_bam2strandness_pie.png")
    
    # 'bbox_inches="tight"' is crucial for removing whitespace for HTML embedding.
    plt.savefig(pie_pdf, format='pdf', bbox_inches='tight')
    plt.savefig(pie_png, format='png', bbox_inches='tight')
    plt.close()

    print(f"Pie chart saved as: {pie_pdf}, {pie_png}")


def main():
    """Main function to parse arguments and run the analysis."""
    parser = argparse.ArgumentParser(
        description="Infer strand specificity and generate a TSV summary and pie chart.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--input_bam", required=True, help="Input BAM file.")
    parser.add_argument("--bed", required=True, help="Reference gene model in BED format.")
    parser.add_argument("--test_read_num", type=int, default=100000000,
                        help="Number of reads to sample for inference (default: 100,000,000).")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to save the output files.")
    args = parser.parse_args()

    run_infer_experiment(
        args.input_bam,
        args.bed,
        args.test_read_num,
        args.output_dir
    )


if __name__ == "__main__":
    main()
