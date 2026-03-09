#!/usr/bin/env python3
"""
Infer strand specificity from a BAM file, generate a TSV summary,
and produce a pie chart.

Also computes Splices per Kilobase (splice/kb) as a DNA contamination
metric. Splice/kb cleanly separates DNase-treated from untreated samples
(DNase ~1.67 vs no-DNase ~0.13 in EV RNA-seq), unlike overall strand
specificity which is confounded by ultra-short insert sizes after DNase.

Metric definition:
  splice/kb = (total splice junction crossings) / (n_unique_reads x avg_mapped_len / 1000)

Both the STAR log path and the BAM fallback path count splice *junction crossings*
(not spliced reads). A read with CIGAR 30M500N30M200N30M crosses 2 junctions
and contributes 2, consistent with STAR "Number of splices: Total".

Since EVscope already incorporates gDNA correction (subtracting opposite-strand
read counts), splice/kb serves as an independent, complementary QC metric.

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

# Compiled regex for splice junction detection in CIGAR
_SPLICE_RE = re.compile(r'(\d+)N')
# Compiled regex for M (match/mismatch) operations in CIGAR
_MATCH_RE  = re.compile(r'(\d+)M')


def count_cigar_junctions(cigar):
    """Count the number of N operations (splice junctions) in a CIGAR string."""
    return len(_SPLICE_RE.findall(cigar))


def cigar_mapped_len(cigar):
    """Sum of M (match/mismatch) operations in CIGAR string.

    This matches STAR 'Average mapped length': only aligned M bases are counted.
    Soft-clipped (S) and intronic (N) bases are excluded, consistent with STAR
    reporting splice junction statistics from uniquely mapped reads only.

    Example: '19M170458N61M' -> 80 nt  (intron N excluded)
             '1S65M'         -> 65 nt  (soft-clip S excluded)
    """
    return sum(int(n) for n in _MATCH_RE.findall(cigar))


def compute_splice_per_kb(bam_file, star_log=None, sample_n=500000):
    """
    Compute splices per kilobase (splice/kb).

    Metric: splice/kb = n_splice_junctions / (n_unique_reads x avg_mapped_len / 1000)

    Both paths count splice *junction crossings* (consistent with STAR
    "Number of splices: Total"). A read spanning 2 introns = 2 junctions.

    Primary path: STAR Log.final.out (fast, exact).
    Fallback: samtools view with up to sample_n uniquely mapped reads (-q 255).

    Returns: (splice_per_kb, n_junctions, n_reads, avg_len) or None on failure.
    """
    if star_log and os.path.isfile(star_log):
        try:
            metrics = {}
            with open(star_log) as f:
                for line in f:
                    if "|" in line:
                        k, v = [x.strip() for x in line.split("|")]
                        metrics[k] = v.replace('%', '').strip()
            n_junctions = int(metrics.get("Number of splices: Total", 0))
            n_unique    = int(metrics.get("Uniquely mapped reads number", 0))
            avg_len     = float(metrics.get("Average mapped length", 0))
            if n_unique > 0 and avg_len > 0:
                splice_per_kb = n_junctions * 1000.0 / (n_unique * avg_len)
                return round(splice_per_kb, 4), n_junctions, n_unique, round(avg_len, 1)
        except Exception as e:
            print(f"Warning: STAR log parse failed ({e}), falling back to BAM sampling.")

    # Fallback: sample uniquely mapped reads from BAM (-q 255 = unique in STAR)
    try:
        cmd_total = [
            "samtools", "view", "-c", "-F", "256", "-q", "255", bam_file
        ]
        r_total = subprocess.run(cmd_total, capture_output=True, text=True, check=True)
        n_total = int(r_total.stdout.strip())
        if n_total == 0:
            return 0.0, 0, 0, 0.0

        frac = min(1.0, sample_n / n_total)
        sample_flag = ["-s", f"{frac:.6f}"] if frac < 1.0 else []
        # -q 255: uniquely mapped only (consistent with STAR path)
        cmd_view = ["samtools", "view", "-F", "256", "-q", "255"] + sample_flag + [bam_file]
        r_view = subprocess.run(cmd_view, capture_output=True, text=True, check=True)

        n_junctions_sampled = 0
        total_len = 0
        n_sampled = 0
        for line in r_view.stdout.splitlines():
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            n_sampled += 1
            n_junctions_sampled += count_cigar_junctions(fields[5])
            total_len += cigar_mapped_len(fields[5])  # M-ops only, matches STAR avg_mapped_len

        if n_sampled == 0:
            return 0.0, 0, n_total, 0.0

        avg_len = total_len / n_sampled
        junction_rate = n_junctions_sampled / n_sampled
        splice_per_kb = junction_rate * 1000.0 / avg_len if avg_len > 0 else 0.0
        n_junctions_est = int(junction_rate * n_total)
        return round(splice_per_kb, 4), n_junctions_est, n_total, round(avg_len, 1)

    except Exception as e:
        print(f"Warning: BAM-based splice/kb computation failed: {e}")
        return None


def run_infer_experiment(bam_file, refgene_bed, test_read_num, output_dir, star_log=None):
    """
    Runs RSeQC infer_experiment.py, parses the output, and generates a TSV
    and pie chart. Also computes and appends splice/kb to the TSV.
    """
    if not os.path.isfile(bam_file):
        print(f"Error: Input BAM file not found at {bam_file}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.basename(bam_file).replace('.bam', '')
    tsv_path  = os.path.join(output_dir, f"{base_name}_bam2strandness.tsv")

    frac_fwd, frac_rev, frac_failed, lib_type = "0", "0", "0", "Unknown"
    spkb, n_junc, n_tot, avg_l = "NA", "NA", "NA", "NA"

    with open(tsv_path, "w") as f:
        f.write(
            "BAM\tdata type\t'1++,1--,2+-,2-+ (forward)'\t'1+-,1-+,2++,2-- (reverse)'\t"
            "%reads failed to determine:\t"
            "Splice_per_kb\tn_splice_junctions\tn_unique_reads\tavg_mapped_len_nt\n"
        )
        try:
            cmd  = ["infer_experiment.py", "-i", bam_file, "-r", refgene_bed,
                    "-s", str(test_read_num)]
            proc = subprocess.run(cmd, check=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            lines       = proc.stdout.strip().split('\n')
            lib_type    = lines[0].strip() if lines else "Library type not found"
            frac_failed = lines[1].split(": ")[1] if len(lines) > 1 else "0"
            frac_fwd    = lines[2].split(": ")[1] if len(lines) > 2 else "0"
            frac_rev    = lines[3].split(": ")[1] if len(lines) > 3 else "0"

            result = compute_splice_per_kb(bam_file, star_log=star_log)
            if result:
                spkb, n_junc, n_tot, avg_l = result

            f.write(
                f"{bam_file}\t{lib_type}\t{frac_fwd}\t{frac_rev}\t{frac_failed}\t"
                f"{spkb}\t{n_junc}\t{n_tot}\t{avg_l}\n"
            )

        except subprocess.CalledProcessError as e:
            error_message = f"Error processing {bam_file}: {e.stderr}"
            f.write(f"{bam_file}\tError\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
                    f"Details: {error_message}\n")
            print(error_message)
            sys.exit(1)

    print(f"Strandness results written to {tsv_path}")
    print(f"Splice/kb = {spkb}  (n_junctions={n_junc}, n_reads={n_tot}, avg_len={avg_l}nt)")

    # --- Pie chart ---
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
    wedges, _ = ax.pie(
        sizes, startangle=90, colors=colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1}, radius=1.2
    )
    ax.axis('equal')
    ax.set_title(
        f"Percentage of strand specificity\n{lib_type}\n{base_name}\n"
        f"Splice/kb = {spkb}",
        fontsize=8, pad=10
    )
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
        description=(
            "Infer strand specificity and generate a TSV summary and pie chart. "
            "Also computes splices per kilobase (splice/kb) as a DNA contamination metric."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--input_bam",     required=True)
    parser.add_argument("--bed",           required=True,
                        help="Non-overlapping exon BED for infer_experiment.py.")
    parser.add_argument("--test_read_num", type=int, default=100000000)
    parser.add_argument("--output_dir",    required=True)
    parser.add_argument("--star_log",      required=False, default=None,
                        help="STAR Log.final.out from Step 6 refined alignment. "
                             "Primary source for splice/kb (fast, exact). "
                             "Falls back to BAM sampling if absent.")
    args = parser.parse_args()
    run_infer_experiment(
        args.input_bam, args.bed, args.test_read_num,
        args.output_dir, star_log=args.star_log
    )


if __name__ == "__main__":
    main()
