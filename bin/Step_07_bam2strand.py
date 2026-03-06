#!/usr/bin/env python3
"""
Infer strand specificity from a BAM file, generate a TSV summary,
and produce a pie chart.

Also computes Splice reads per kilobase (splice/kb) as a DNA contamination
metric. Splice/kb cleanly separates DNase-treated from untreated samples
(DNase ~1.61 vs no-DNase ~0.19 in EV RNA-seq), unlike overall strand
specificity which is confounded by ultra-short insert sizes after DNase
treatment.

Since EVscope already incorporates a gDNA correction method (subtracting
opposite-strand read counts), splice/kb serves as an independent,
complementary QC metric for residual DNA contamination.

Reference: EVscope (https://www.biorxiv.org/content/10.1101/2025.06.24.660984v1)

The default test read number is 100,000,000.
"""

import subprocess
import sys
import os
import re
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set font properties for publication quality
plt.rcParams["font.family"] = "Arial"
mpl.rcParams['pdf.fonttype'] = 42


def compute_splice_per_kb(bam_file, star_log=None, sample_n=500000):
    """
    Compute splice reads per kilobase (splice/kb).

    Formula: splice/kb = (n_spliced_reads / n_total_reads) * (1000 / avg_mapped_len)

    This metric separates DNase-treated (~1.61 splice/kb) from
    untreated (~0.19 splice/kb) EV RNA-seq libraries.

    If star_log is provided, uses STAR-reported splice count, unique reads,
    and average mapped length. Otherwise, samples the BAM directly.

    Returns: (splice_per_kb, n_spliced, n_total, avg_len) or None on failure.
    """
    if star_log and os.path.isfile(star_log):
        # Parse from STAR Log.final.out
        metrics = {}
        try:
            with open(star_log) as f:
                for line in f:
                    if "|" in line:
                        k, v = [x.strip() for x in line.split("|")]
                        metrics[k] = v.replace('%', '').strip()
            n_splices   = int(metrics.get("Number of splices: Total", 0))
            n_unique    = int(metrics.get("Uniquely mapped reads number", 0))
            avg_len     = float(metrics.get("Average mapped length", 0))
            if n_unique > 0 and avg_len > 0:
                splice_per_kb = n_splices * 1000.0 / (n_unique * avg_len)
                return round(splice_per_kb, 4), n_splices, n_unique, round(avg_len, 1)
        except Exception as e:
            print(f"Warning: STAR log parse failed ({e}), falling back to BAM sampling.")

    # Fall back: sample reads directly from BAM
    try:
        cmd_total = ["samtools", "view", "-c", "-F", "256", bam_file]
        r_total = subprocess.run(cmd_total, capture_output=True, text=True, check=True)
        n_total = int(r_total.stdout.strip())

        # Sample up to sample_n reads for efficiency
        frac = min(1.0, sample_n / max(n_total, 1))
        sample_flag = ["-s", f"{frac:.6f}"] if frac < 1.0 else []
        cmd_view = ["samtools", "view", "-F", "256"] + sample_flag + [bam_file]
        r_view = subprocess.run(cmd_view, capture_output=True, text=True, check=True)

        lines = [l for l in r_view.stdout.strip().split("\n") if l]
        if not lines:
            return None

        n_sampled = len(lines)
        n_spliced_sampled = 0
        total_len = 0
        for l in lines:
            fields = l.split("\t")
            if len(fields) < 10:
                continue
            if re.search(r'\d+N', fields[5]):
                n_spliced_sampled += 1
            total_len += len(fields[9])

        avg_len = total_len / n_sampled if n_sampled > 0 else 0
        splice_rate = n_spliced_sampled / n_sampled if n_sampled > 0 else 0
        splice_per_kb = splice_rate * 1000.0 / avg_len if avg_len > 0 else 0

        # Scale n_spliced back to full BAM
        n_spliced_est = int(splice_rate * n_total)
        return round(splice_per_kb, 4), n_spliced_est, n_total, round(avg_len, 1)

    except Exception as e:
        print(f"Warning: BAM-based splice/kb computation failed: {e}")
        return None


def run_infer_experiment(bam_file, refgene_bed, test_read_num, output_dir, star_log=None):
    """
    Runs RSeQC's infer_experiment.py, parses the output, and generates
    a TSV file and a formatted pie chart. Also computes and saves splice/kb.
    """
    if not os.path.isfile(bam_file):
        print(f"Error: Input BAM file not found at {bam_file}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.basename(bam_file).replace('.bam', '')

    # --- Run infer_experiment.py and write to TSV ---
    tsv_path = os.path.join(output_dir, f"{base_name}_bam2strandness.tsv")
    frac_fwd, frac_rev, frac_failed, lib_type = "0", "0", "0", "Unknown"

    with open(tsv_path, "w") as f:
        f.write("BAM\tdata type\t'1++,1--,2+-,2-+ (forward)'\t'1+-,1-+,2++,2-- (reverse)'\t"
                "%reads failed to determine:\tSplice_per_kb\tn_spliced_reads\tn_total_reads\tavg_mapped_len_nt\n")
        try:
            cmd = ["infer_experiment.py", "-i", bam_file, "-r", refgene_bed, "-s", str(test_read_num)]
            proc = subprocess.run(cmd, check=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            lines = proc.stdout.strip().split('\n')
            lib_type    = lines[0].strip() if lines else "Library type not found"
            frac_failed = lines[1].split(": ")[1] if len(lines) > 1 else "0"
            frac_fwd    = lines[2].split(": ")[1] if len(lines) > 2 else "0"
            frac_rev    = lines[3].split(": ")[1] if len(lines) > 3 else "0"

            # Compute splice/kb
            splice_result = compute_splice_per_kb(bam_file, star_log=star_log)
            if splice_result:
                spkb, n_spl, n_tot, avg_l = splice_result
            else:
                spkb, n_spl, n_tot, avg_l = "NA", "NA", "NA", "NA"

            f.write(f"{bam_file}\t{lib_type}\t{frac_fwd}\t{frac_rev}\t{frac_failed}\t"
                    f"{spkb}\t{n_spl}\t{n_tot}\t{avg_l}\n")

        except subprocess.CalledProcessError as e:
            error_message = f"Error processing {bam_file}: {e.stderr}"
            f.write(f"{bam_file}\tError\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
                    f"Details: {error_message}\n")
            print(error_message)
            sys.exit(1)

    print(f"Strandness results written to {tsv_path}")
    print(f"Splice/kb = {spkb}")

    # --- Generate Publication-Quality Pie Chart ---
    try:
        f_failed  = float(frac_failed)
        f_forward = float(frac_fwd)
        f_reverse = float(frac_rev)
    except ValueError:
        print("Error: Could not convert strand fractions to numbers. Skipping chart generation.")
        return

    labels = ['Forward strand', 'Reverse strand', 'Undetermined strand']
    sizes  = [f_forward, f_reverse, f_failed]
    colors = ['#FF4500', '#FFA500', '#87CEEB']

    fig, ax = plt.subplots(figsize=(5, 3), dpi=300)
    wedges, texts = ax.pie(
        sizes, startangle=90, colors=colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1}, radius=1.2
    )
    ax.axis('equal')

    title_str = (f"Percentage of strand specificity\n{lib_type}\n{base_name}\n"
                 f"Splice/kb = {spkb}")
    ax.set_title(title_str, fontsize=8, pad=10)

    legend_labels = [f"{lbl} ({sz*100:.1f}%)" for lbl, sz in zip(labels, sizes)]
    ax.legend(wedges, legend_labels, loc='center left',
              bbox_to_anchor=(0.80, 0.5), frameon=False, prop={'size': 7})

    pie_pdf = os.path.join(output_dir, f"{base_name}_bam2strandness_pie.pdf")
    pie_png = os.path.join(output_dir, f"{base_name}_bam2strandness_pie.png")
    plt.savefig(pie_pdf, format='pdf', bbox_inches='tight')
    plt.savefig(pie_png, format='png', bbox_inches='tight')
    plt.close()

    print(f"Pie chart saved as: {pie_pdf}, {pie_png}")


def main():
    parser = argparse.ArgumentParser(
        description="Infer strand specificity and generate a TSV summary and pie chart. "
                    "Also computes splice reads per kilobase (splice/kb) as a DNA contamination metric.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--input_bam",    required=True,  help="Input BAM file.")
    parser.add_argument("--bed",          required=True,  help="Reference gene model in BED format (non-overlapping exon).")
    parser.add_argument("--test_read_num", type=int, default=100000000,
                        help="Number of reads to sample for strand inference (default: 100,000,000).")
    parser.add_argument("--output_dir",   required=True,  help="Directory to save the output files.")
    parser.add_argument("--star_log",     required=False, default=None,
                        help="STAR Log.final.out from refined alignment (Step 6). "
                             "Used for splice/kb computation. If omitted, BAM is sampled directly.")
    args = parser.parse_args()

    run_infer_experiment(
        args.input_bam,
        args.bed,
        args.test_read_num,
        args.output_dir,
        star_log=args.star_log
    )


if __name__ == "__main__":
    main()
