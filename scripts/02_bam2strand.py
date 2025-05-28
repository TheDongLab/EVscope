#!/usr/bin/env python
"""
Infer strand specificity from a BAM file, generate a TSV summary,
and produce a pie chart showing percentages of forward, reverse,
and undetermined reads. The pie chart is saved in PDF, PNG, and SVG
formats. The default test read number is 100000000.
"""

import subprocess
import sys
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set SVG output to use vector fonts and set global font to Arial
mpl.rcParams['svg.fonttype'] = 'none'
plt.rcParams["font.family"] = "Arial"
# Set global font size to 6 for all text elements in the figure

def run_infer_experiment(bam_file, refgene_bed, test_read_num):
    # Check if the BAM file exists
    if not os.path.isfile(bam_file):
        print(f"Error: The file {bam_file} does not exist.")
        sys.exit(1)

    # Derive base name and output directory (use current working directory)
    base_name = os.path.basename(bam_file).replace('.bam', '')
    current_dir = os.getcwd()
    output_file = os.path.join(current_dir, f"{base_name}_bam2strandness.tsv")

    # Open TSV output file and write header
    with open(output_file, "w") as file:
        file.write("BAM\tdata type\t'1++,1--,2+-,2-+ (forward)'\t'1+-,1-+,2++,2-- (reverse)'\t%reads failed to determine:\n")

        try:
            # Build and run the infer_experiment command
            command = f"infer_experiment.py -i {bam_file} -r {refgene_bed} -s {test_read_num}"
            process = subprocess.run(command, shell=True, check=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            result_lines = process.stdout.strip().split('\n')

            # Parse command output
            library_type = result_lines[0].strip() if result_lines else "Library type not found"
            fraction_failed = result_lines[1].split(": ")[1] if len(result_lines) > 1 else "0"
            fraction_forward = result_lines[2].split(": ")[1] if len(result_lines) > 2 else "0"
            fraction_reverse = result_lines[3].split(": ")[1] if len(result_lines) > 3 else "0"

            # Write results to TSV file
            file.write(f"{bam_file}\t{library_type}\t{fraction_forward}\t{fraction_reverse}\t{fraction_failed}\n")

        except subprocess.CalledProcessError as e:
            file.write(f"{bam_file}\tError\tExit code: {e.returncode}, Error message: {e.stderr}\n")
            print(f"Error processing {bam_file}: {e.stderr}")
            sys.exit(1)

    print(f"Results have been written to {output_file}")

    # Convert fraction values to float for plotting
    try:
        f_failed = float(fraction_failed)
        f_forward = float(fraction_forward)
        f_reverse = float(fraction_reverse)
    except ValueError:
        print("Error: Unable to convert fraction values to float.")
        sys.exit(1)

    # Prepare pie chart data and labels
    labels = ['Forward strand', 'Reverse strand', 'Undetermined strand']
    sizes = [f_forward, f_reverse, f_failed]

    # Use only three colors: red, orange, skyblue
    colors = ['red', 'orange', 'skyblue']

    # Create pie chart with a larger figure size and improved aesthetics
    fig, ax = plt.subplots(figsize=(3, 3), dpi=300)
    # Draw pie chart without autopct
    wedges, texts = ax.pie(
        sizes,
        startangle=90,
        colors=colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1},
        textprops={'fontsize': 6},
        radius=0.2
    )

    ax.axis('equal')  # Ensures the pie is a circle
    plt.title(f"Percentage of read strand specificity\n{library_type}:{base_name}", fontsize=8)

    # Create legend labels with percentage appended to the right side of text
    legend_labels = [f"{label} ({size*100:.1f}%)" for label, size in zip(labels, sizes)]
    # Add legend at the bottom with reduced spacing between entries and closer to the pie chart
    ax.legend(wedges, legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.05),
              ncol=3, frameon=False, prop={'size': 6}, columnspacing=0.5)

    # Adjust layout to ensure legend is not clipped and is closer to the pie chart
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    # Define output filenames for pie chart in different formats (in current directory)
    pie_pdf = os.path.join(current_dir, f"{base_name}_bam2strandness_pie.pdf")
    pie_png = os.path.join(current_dir, f"{base_name}_bam2strandness_pie.png")
    pie_svg = os.path.join(current_dir, f"{base_name}_bam2strandness_pie.svg")

    # Save pie chart in PDF, PNG, and SVG formats with extra bounding
    plt.savefig(pie_pdf, format='pdf', bbox_inches='tight')
    plt.savefig(pie_png, format='png', bbox_inches='tight')
    plt.savefig(pie_svg, format='svg', bbox_inches='tight')
    plt.close()

    print(f"Pie chart saved as: {pie_pdf}, {pie_png}, {pie_svg}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Infer strand specificity from a BAM file and generate TSV and pie chart outputs."
    )
    parser.add_argument("--input_bam", required=True, help="Input BAM file")
    parser.add_argument("--bed", required=True, help="Reference gene model BED file")
    parser.add_argument("--test_read_num", type=int, default=100000000,
                        help="Number of test reads (default: 100000000)")
    args = parser.parse_args()

    run_infer_experiment(args.input_bam, args.bed, args.test_read_num)

if __name__ == "__main__":
    main()

