#!/usr/bin/env python
"""
Infer strand specificity from a BAM file, generate a TSV summary,
and produce a pie chart showing percentages of forward, reverse,
and undetermined reads. The pie chart is saved in PDF and PNG
formats. The default test read number is 100000000.
"""

import subprocess
import sys
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.rcParams["font.family"] = "Arial"
mpl.rcParams['pdf.fonttype'] = 42


def run_infer_experiment(bam_file, refgene_bed, test_read_num, output_dir):
    # Ensure BAM exists
    if not os.path.isfile(bam_file):
        print(f"Error: {bam_file} not found.")
        sys.exit(1)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.basename(bam_file).replace('.bam', '')

    # TSV output
    tsv_path = os.path.join(output_dir, f"{base_name}_bam2strandness.tsv")
    with open(tsv_path, "w") as f:
        f.write("BAM\tdata type\t'1++,1--,2+-,2-+ (forward)'\t'1+-,1-+,2++,2-- (reverse)'\t%reads failed to determine:\n")
        try:
            cmd = f"infer_experiment.py -i {bam_file} -r {refgene_bed} -s {test_read_num}"
            proc = subprocess.run(cmd, shell=True, check=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            lines = proc.stdout.strip().split('\n')
            lib_type = lines[0].strip() if lines else "Library type not found"
            frac_failed = lines[1].split(": ")[1] if len(lines) > 1 else "0"
            frac_fwd = lines[2].split(": ")[1] if len(lines) > 2 else "0"
            frac_rev = lines[3].split(": ")[1] if len(lines) > 3 else "0"
            f.write(f"{bam_file}\t{lib_type}\t{frac_fwd}\t{frac_rev}\t{frac_failed}\n")
        except subprocess.CalledProcessError as e:
            f.write(f"{bam_file}\tError\tExit code: {e.returncode}, Error message: {e.stderr}\n")
            print(f"Error processing {bam_file}: {e.stderr}")
            sys.exit(1)

    print(f"Results written to {tsv_path}")

    # Convert to float
    try:
        f_failed = float(frac_failed)
        f_forward = float(frac_fwd)
        f_reverse = float(frac_rev)
    except ValueError:
        print("Error: Unable to convert fractions.")
        sys.exit(1)

    # Pie chart data
    labels = ['Forward strand', 'Reverse strand', 'Undetermined strand']
    sizes = [f_forward, f_reverse, f_failed]
    colors = ['red', 'orange', 'skyblue']

    fig, ax = plt.subplots(figsize=(3, 3), dpi=300)
    wedges, _ = ax.pie(
        sizes,
        startangle=90,
        colors=colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1},
        textprops={'fontsize': 6},
        radius=0.2
    )

    ax.axis('equal')
    plt.title(f"Percentage of read strand specificity\n{lib_type}:{base_name}", fontsize=8)

    legend_labels = [f"{lbl} ({sz*100:.1f}%)" for lbl, sz in zip(labels, sizes)]
    ax.legend(wedges, legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.05),
              ncol=3, frameon=False, prop={'size': 6}, columnspacing=0.5)

    plt.tight_layout(rect=[0, 0.05, 1, 1])

    # Save pie charts
    pie_pdf = os.path.join(output_dir, f"{base_name}_bam2strandness_pie.pdf")
    pie_png = os.path.join(output_dir, f"{base_name}_bam2strandness_pie.png")
    plt.savefig(pie_pdf, format='pdf', bbox_inches='tight')
    plt.savefig(pie_png, format='png', bbox_inches='tight')
    plt.close()

    print(f"Pie chart saved as: {pie_pdf}, {pie_png}")


def main():
    parser = argparse.ArgumentParser(
        description="Infer strand specificity and generate outputs."
    )
    parser.add_argument("--input_bam", required=True, help="Input BAM file")
    parser.add_argument("--bed", required=True, help="Reference gene model BED file")
    parser.add_argument("--test_read_num", type=int, default=100000000,
                        help="Number of test reads")
    parser.add_argument("--output_dir", required=True,
                        help="Directory for outputs")
    args = parser.parse_args()

    run_infer_experiment(
        args.input_bam,
        args.bed,
        args.test_read_num,
        args.output_dir
    )


if __name__ == "__main__":
    main()

