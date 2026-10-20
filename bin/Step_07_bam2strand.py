#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
Infer strand specificity from a BAM file, generate strandness summaries,
and produce publication-quality QC figures.

Step 07 keeps the legacy RSeQC infer_experiment output unchanged and adds:
  1) an RSeQC-compatible read-length-stratified strandness TSV;
  2) a read-length-stratified strandness line + point plot;
  3) an independent splice junction density plot.

Read-length stratification uses cumulative thresholds over RSeQC-style
alignment records. The X >= 0 point is anchored to the legacy
*_bam2strandness.tsv fractions so it exactly reproduces the original
RSeQC infer_experiment summary.

Reference: EVscope (https://www.biorxiv.org/content/10.1101/2025.06.24.660984v1)
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

import matplotlib as mpl
mpl.use("pdf")
import matplotlib.pyplot as plt
import pysam
from bx.intervals import Intersecter

# Publication-quality vector output defaults.
mpl.rcParams.update({
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 6,
    'axes.linewidth': 0.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 2,
    'ytick.major.size': 2,
    'lines.linewidth': 0.5,
})

_SPLICE_RE = re.compile(r'(\d+)N')
_MATCH_RE = re.compile(r'(\d+)[M=X]')

NO_DNASE_SPLICE_PER_KB = 0.13
DNASE_TREATED_SPLICE_PER_KB_VALUES = [1.67, 0.63, 1.64, 0.94]
MEAN_DNASE_TREATED_SPLICE_PER_KB = sum(DNASE_TREATED_SPLICE_PER_KB_VALUES) / len(DNASE_TREATED_SPLICE_PER_KB_VALUES)

# Keep the legacy strandness pie colors unchanged. The length-stratified
# line plot intentionally reuses the exact same colors.
STRAND_COLORS = {
    'Forward strand': '#FF4500',
    'Reverse strand': '#FFA500',
    'Undetermined': '#87CEEB',
}

DEFAULT_LENGTH_THRESHOLDS = list(range(0, 226, 25))
LENGTH_THRESHOLD_STEP = 25


def _safe_float(value):
    """Return float(value), or None for NA/non-finite plotting inputs."""
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if x != x or x in (float('inf'), float('-inf')):
        return None
    return x


def clean_sample_name_from_bam_basename(base_name):
    """Convert EVscope BAM basename to display sample name for plotting."""
    suffixes = [
        '_STAR_umi_dedup_Aligned.sortedByCoord.out',
        '_STAR_no_dedup_Aligned.sortedByCoord.out',
        '_Aligned.sortedByCoord.out',
        '_Aligned.sortedByCoord_umi_dedup.out',
        '_Aligned.sortedByCoord_no_dedup.out',
    ]
    sample_name = base_name
    for suffix in suffixes:
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]
            break
    return sample_name


def count_cigar_junctions(cigar):
    """Count the number of N operations (splice junctions) in a CIGAR string."""
    return len(_SPLICE_RE.findall(cigar or ''))


def cigar_mapped_len(cigar):
    """Sum aligned M/=/X operations in a CIGAR string."""
    return sum(int(n) for n in _MATCH_RE.findall(cigar or ''))


def plot_strandness_pie(output_dir, base_name, lib_type, frac_fwd, frac_rev, frac_failed):
    """Plot the legacy strand specificity pie chart as a pure pie chart."""
    f_failed = float(frac_failed)
    f_forward = float(frac_fwd)
    f_reverse = float(frac_rev)

    labels = ['Forward strand', 'Reverse strand', 'Undetermined']
    sizes = [f_forward, f_reverse, f_failed]
    colors = [STRAND_COLORS['Forward strand'], STRAND_COLORS['Reverse strand'], STRAND_COLORS['Undetermined']]

    fig, ax = plt.subplots(figsize=(3.2, 3.2), dpi=300)
    wedges, _ = ax.pie(
        sizes,
        startangle=90,
        colors=colors,
        wedgeprops={'edgecolor': 'white', 'linewidth': 0.8},
        radius=1.05,
    )
    ax.axis('equal')
    ax.set_title(f'Strand Specificity\n({lib_type})', fontsize=8, pad=8)
    legend_labels = [f'{label} ({size * 100:.1f}%)' for label, size in zip(labels, sizes)]
    ax.legend(
        wedges,
        legend_labels,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.18),
        frameon=False,
        fontsize=6,
        ncol=1,
        handlelength=1.0,
        handletextpad=0.4,
    )
    fig.subplots_adjust(left=0.06, right=0.94, top=0.86, bottom=0.24)

    pie_pdf = os.path.join(output_dir, f'{base_name}_bam2strandness_pie.pdf')
    pie_png = os.path.join(output_dir, f'{base_name}_bam2strandness_pie.png')
    fig.savefig(pie_pdf, format='pdf')
    fig.savefig(pie_png, format='png', dpi=300)
    plt.close(fig)
    print(f'Strandness pie chart saved as: {pie_pdf}, {pie_png}')


def plot_splice_junction_density(output_dir, base_name, spkb):
    """Plot splice/kb QC as an independent reference barplot."""
    sample_spkb = _safe_float(spkb)
    display_sample_name = clean_sample_name_from_bam_basename(base_name)
    bar_labels = [display_sample_name, 'DNase-treated reference', 'No-DNase reference']
    bar_values = [
        sample_spkb if sample_spkb is not None else 0.0,
        MEAN_DNASE_TREATED_SPLICE_PER_KB,
        NO_DNASE_SPLICE_PER_KB,
    ]
    bar_colors = ['#3B82B6', '#4DAF7C', '#D65F5F']

    fig, ax = plt.subplots(figsize=(5.9, 2.2), dpi=300)
    y_pos = [0.00, 0.42, 0.84]
    ax.barh(y_pos, bar_values, color=bar_colors, edgecolor='white', linewidth=0.35, height=0.22)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(bar_labels, fontsize=6.5)
    ax.tick_params(axis='y', pad=3, length=0)
    ax.set_xlabel('Splice junctions per kb aligned sequence', fontsize=7, labelpad=5)
    ax.set_title('Splice junction density (SJ/kb)', fontsize=8.5, pad=8)

    xmax = max(1.5, max(bar_values) * 1.24)
    dx = xmax * 0.018
    ax.set_xlim(0, xmax)
    ax.set_ylim(1.15, -0.25)
    ax.grid(axis='x', linestyle=':', linewidth=0.4, alpha=0.28)
    ax.set_axisbelow(True)

    for y, value, label, color in zip(y_pos, bar_values, bar_labels, bar_colors):
        if label == display_sample_name and sample_spkb is None:
            ax.text(dx, y, 'NA', va='center', ha='left', fontsize=6.3, color='#666666', fontweight='bold')
        else:
            ax.text(value + dx, y, f'{value:.2f}', va='center', ha='left', fontsize=6.3, color=color, fontweight='bold')

    fig.subplots_adjust(left=0.30, right=0.94, top=0.82, bottom=0.25)
    out_pdf = os.path.join(output_dir, f'{base_name}_splice_junction_density.pdf')
    out_png = os.path.join(output_dir, f'{base_name}_splice_junction_density.png')
    fig.savefig(out_pdf, format='pdf')
    fig.savefig(out_png, format='png', dpi=300)
    plt.close(fig)
    print(f'Splice junction density plot saved as: {out_pdf}, {out_png}')


def rseqc_qlen(aligned_read):
    """Return the RSeQC-style qlen used by infer_experiment.py.

    RSeQC uses aligned_read.qlen both for overlap and as the effective query
    length. In modern pysam this is closest to query_alignment_length and is an
    absolute nucleotide count for the alignment record.
    """
    qlen = getattr(aligned_read, 'qlen', None)
    if qlen is None:
        qlen = aligned_read.query_alignment_length
    if qlen is None:
        qlen = aligned_read.query_length
    return qlen


def _legacy_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _pct_from_fraction(frac):
    return round(float(frac) * 100.0, 4)


def _fraction_from_counts(count, total):
    return round(count / float(total), 4) if total > 0 else 'NA'


def _pct_from_counts(count, total):
    return round(count * 100.0 / float(total), 4) if total > 0 else 'NA'


def _format_fraction4(value):
    try:
        return f'{float(value):.4f}'
    except (TypeError, ValueError):
        return 'NA'


def build_length_thresholds(lengths):
    """Return cumulative thresholds with N > 0, always retaining 0 nt."""
    max_len = max(lengths) if lengths else 0
    thresholds = [t for t in DEFAULT_LENGTH_THRESHOLDS if t == 0 or t <= max_len]
    threshold = DEFAULT_LENGTH_THRESHOLDS[-1] + LENGTH_THRESHOLD_STEP
    while threshold <= max_len:
        thresholds.append(threshold)
        threshold += LENGTH_THRESHOLD_STEP
    return thresholds


def load_refgene_intervals(refbed):
    """Load a 12-column BED gene model into per-chromosome strand intervals."""
    gene_ranges = {}
    with open(refbed, 'r') as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.rstrip('\n').split()
            if len(fields) < 6:
                continue
            chrom = fields[0]
            try:
                tx_start = int(fields[1])
                tx_end = int(fields[2])
            except ValueError:
                continue
            strand = fields[5]
            if strand not in ('+', '-'):
                continue
            gene_ranges.setdefault(chrom, Intersecter()).insert(tx_start, tx_end, strand)
    return gene_ranges


def classify_strandness(aligned_read, gene_strands):
    """Classify one alignment record using RSeQC infer_experiment key semantics."""
    # RSeQC uses ':'.join(set(...)); recognized forward/reverse keys only occur
    # for single-strand overlaps, so sorted order does not affect informative calls.
    strand_from_gene = ':'.join(sorted(set(gene_strands)))
    map_strand = '-' if aligned_read.is_reverse else '+'

    if aligned_read.is_paired:
        if aligned_read.is_read1:
            read_id = '1'
        elif aligned_read.is_read2:
            read_id = '2'
        else:
            return 'Undetermined'
        key = read_id + map_strand + strand_from_gene
        if key in {'1++', '1--', '2+-', '2-+'}:
            return 'Forward strand'
        if key in {'1+-', '1-+', '2++', '2--'}:
            return 'Reverse strand'
        return 'Undetermined'

    key = map_strand + strand_from_gene
    if key in {'++', '--'}:
        return 'Forward strand'
    if key in {'+-', '-+'}:
        return 'Reverse strand'
    return 'Undetermined'


def collect_rseqc_usable_records(bam_file, refgene_bed, sample_size, min_mapq=30):
    """Collect the same sampled universe used by RSeQC infer_experiment.

    Unit is one read alignment record. Paired-end R1/R2 records are counted
    separately, matching the legacy *_bam2strandness.tsv semantics.
    """
    gene_ranges = load_refgene_intervals(refgene_bed)
    records = []
    metadata = {
        'n_total_records_seen': 0,
        'n_after_filter': 0,
        'n_missing_qlen': 0,
        'n_without_refgene_overlap': 0,
        'n_usable_reads_sampled': 0,
    }

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for aligned_read in bam:
            metadata['n_total_records_seen'] += 1
            if metadata['n_usable_reads_sampled'] >= sample_size:
                break
            if aligned_read.is_qcfail or aligned_read.is_duplicate:
                continue
            # Match RSeQC 5.x: skip unmapped and secondary, but do not add a
            # supplementary-alignment exclusion here.
            if aligned_read.is_unmapped or aligned_read.is_secondary:
                continue
            if aligned_read.mapping_quality < min_mapq:
                continue
            metadata['n_after_filter'] += 1

            qlen = rseqc_qlen(aligned_read)
            if qlen is None or qlen <= 0:
                metadata['n_missing_qlen'] += 1
                continue

            chrom = bam.get_reference_name(aligned_read.reference_id)
            if chrom not in gene_ranges:
                metadata['n_without_refgene_overlap'] += 1
                continue
            read_start = aligned_read.reference_start
            read_end = read_start + qlen
            gene_strands = set(gene_ranges[chrom].find(read_start, read_end))
            if not gene_strands:
                metadata['n_without_refgene_overlap'] += 1
                continue

            category = classify_strandness(aligned_read, gene_strands)
            records.append((int(qlen), category))
            metadata['n_usable_reads_sampled'] += 1

    return records, metadata


def summarize_length_subset(records, threshold_nt):
    subset = [(length, category) for length, category in records if length >= threshold_nt]
    n_forward = sum(1 for _, category in subset if category == 'Forward strand')
    n_reverse = sum(1 for _, category in subset if category == 'Reverse strand')
    n_undetermined = sum(1 for _, category in subset if category == 'Undetermined')
    n_evaluated = n_forward + n_reverse + n_undetermined
    n_informative = n_forward + n_reverse
    lengths = [length for length, _ in subset]
    mean_len = sum(lengths) / len(lengths) if lengths else 'NA'
    median_len = 'NA'
    if lengths:
        sorted_lengths = sorted(lengths)
        mid = len(sorted_lengths) // 2
        if len(sorted_lengths) % 2:
            median_len = sorted_lengths[mid]
        else:
            median_len = (sorted_lengths[mid - 1] + sorted_lengths[mid]) / 2
    return {
        'n_evaluated': n_evaluated,
        'n_forward': n_forward,
        'n_reverse': n_reverse,
        'n_undetermined': n_undetermined,
        'n_informative': n_informative,
        'mean_read_length_nt': round(mean_len, 2) if mean_len != 'NA' else 'NA',
        'median_read_length_nt': round(median_len, 2) if median_len != 'NA' else 'NA',
    }


def compute_length_stratified_strandness(
    bam_file,
    refgene_bed,
    sample_size,
    min_mapq=30,
    legacy_fractions=None,
):
    """Compute RSeQC-compatible cumulative read-length strandness rows.

    threshold_nt == 0 is anchored to legacy RSeQC fractions when supplied.
    Other thresholds are subsets of the same RSeQC-usable sampled records.
    """
    records, metadata = collect_rseqc_usable_records(bam_file, refgene_bed, sample_size, min_mapq=min_mapq)
    thresholds = build_length_thresholds([length for length, _ in records])
    legacy_fractions = legacy_fractions or {}

    all_counts = summarize_length_subset(records, 0)
    collector_total = all_counts['n_evaluated']
    collector_forward_fraction = _fraction_from_counts(all_counts['n_forward'], collector_total)
    collector_reverse_fraction = _fraction_from_counts(all_counts['n_reverse'], collector_total)
    collector_undetermined_fraction = _fraction_from_counts(all_counts['n_undetermined'], collector_total)

    legacy_forward = legacy_fractions.get('forward')
    legacy_reverse = legacy_fractions.get('reverse')
    legacy_undetermined = legacy_fractions.get('undetermined')
    has_legacy_anchor = legacy_forward is not None and legacy_reverse is not None and legacy_undetermined is not None

    concordance_status = 'not_checked'
    if has_legacy_anchor:
        concordance_status = 'pass' if (
            _format_fraction4(legacy_forward) == _format_fraction4(collector_forward_fraction)
            and _format_fraction4(legacy_reverse) == _format_fraction4(collector_reverse_fraction)
            and _format_fraction4(legacy_undetermined) == _format_fraction4(collector_undetermined_fraction)
        ) else 'fail'

    rows = []
    for threshold_nt in thresholds:
        stats = summarize_length_subset(records, threshold_nt)
        n_evaluated = stats['n_evaluated']
        n_forward = stats['n_forward']
        n_reverse = stats['n_reverse']
        n_undetermined = stats['n_undetermined']
        n_informative = stats['n_informative']

        fraction_forward = _fraction_from_counts(n_forward, n_evaluated)
        fraction_reverse = _fraction_from_counts(n_reverse, n_evaluated)
        fraction_undetermined = _fraction_from_counts(n_undetermined, n_evaluated)
        pct_forward_all = _pct_from_counts(n_forward, n_evaluated)
        pct_reverse_all = _pct_from_counts(n_reverse, n_evaluated)
        pct_undetermined_all = _pct_from_counts(n_undetermined, n_evaluated)
        source = 'evscope_rseqc_compatible_collector'
        legacy_anchor = 'no'

        if threshold_nt == 0 and has_legacy_anchor:
            fraction_forward = round(float(legacy_forward), 4)
            fraction_reverse = round(float(legacy_reverse), 4)
            fraction_undetermined = round(float(legacy_undetermined), 4)
            pct_forward_all = _pct_from_fraction(legacy_forward)
            pct_reverse_all = _pct_from_fraction(legacy_reverse)
            pct_undetermined_all = _pct_from_fraction(legacy_undetermined)
            source = 'legacy_rseqc_bam2strandness_tsv'
            legacy_anchor = 'yes'

        if n_evaluated > 0:
            informative_rate = round(n_informative * 100.0 / n_evaluated, 4)
        else:
            informative_rate = 'NA'
        if n_informative > 0:
            pct_forward_informative = round(n_forward * 100.0 / n_informative, 4)
            pct_reverse_informative = round(n_reverse * 100.0 / n_informative, 4)
            dominant_strand = 'Forward strand' if n_forward >= n_reverse else 'Reverse strand'
            dominant_pct_informative = max(pct_forward_informative, pct_reverse_informative)
        else:
            pct_forward_informative = pct_reverse_informative = dominant_pct_informative = 'NA'
            dominant_strand = 'NA'

        label = f'Aligned query length >= {threshold_nt} nt'
        rows.append({
            'threshold_nt': threshold_nt,
            'threshold_label': f'{threshold_nt} nt',
            'length_label': label,
            'threshold_type': 'ge',
            'n_evaluated': n_evaluated,
            'n_forward': n_forward,
            'n_reverse': n_reverse,
            'n_undetermined': n_undetermined,
            'n_informative': n_informative,
            'fraction_forward': fraction_forward,
            'fraction_reverse': fraction_reverse,
            'fraction_undetermined': fraction_undetermined,
            'pct_forward_all': pct_forward_all,
            'pct_reverse_all': pct_reverse_all,
            'pct_undetermined_all': pct_undetermined_all,
            'pct_forward_informative': pct_forward_informative,
            'pct_reverse_informative': pct_reverse_informative,
            'informative_rate': informative_rate,
            'dominant_strand': dominant_strand,
            'dominant_pct_informative': dominant_pct_informative,
            'mean_read_length_nt': stats['mean_read_length_nt'],
            'median_read_length_nt': stats['median_read_length_nt'],
            'source': source,
            'legacy_anchor': legacy_anchor,
            'concordance_status': concordance_status,
            'collector_fraction_forward': collector_forward_fraction,
            'collector_fraction_reverse': collector_reverse_fraction,
            'collector_fraction_undetermined': collector_undetermined_fraction,
            'min_mapq': min_mapq,
            'sample_size': sample_size,
            'exclude_qcfail': 'yes',
            'exclude_duplicate': 'yes',
            'exclude_secondary': 'yes',
            'exclude_supplementary': 'no_rseqc_compatible',
            'length_metric': 'rseqc_qlen_nt',
            'unit': 'read_alignment_record',
            'rseqc_compatible': 'yes',
            **metadata,
        })
    return rows


def write_length_stratified_tsv(rows, output_path):
    """Write read-length-stratified strandness rows to TSV."""
    fieldnames = [
        'threshold_nt', 'threshold_label', 'length_label', 'threshold_type',
        'n_evaluated', 'n_forward', 'n_reverse', 'n_undetermined', 'n_informative',
        'fraction_forward', 'fraction_reverse', 'fraction_undetermined',
        'pct_forward_all', 'pct_reverse_all', 'pct_undetermined_all',
        'pct_forward_informative', 'pct_reverse_informative', 'informative_rate',
        'dominant_strand', 'dominant_pct_informative',
        'mean_read_length_nt', 'median_read_length_nt',
        'source', 'legacy_anchor', 'concordance_status',
        'collector_fraction_forward', 'collector_fraction_reverse', 'collector_fraction_undetermined',
        'min_mapq', 'sample_size',
        'exclude_qcfail', 'exclude_duplicate', 'exclude_secondary', 'exclude_supplementary',
        'length_metric', 'unit', 'rseqc_compatible',
        'n_total_records_seen', 'n_after_filter', 'n_missing_qlen',
        'n_without_refgene_overlap', 'n_usable_reads_sampled',
    ]
    with open(output_path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _num_or_none(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def plot_length_stratified_strandness(rows, output_dir, base_name):
    """Plot read-length-stratified strandness as line + point plot."""
    thresholds = [int(row['threshold_nt']) for row in rows]
    labels = [str(t) for t in thresholds]
    forward = [_num_or_none(row['pct_forward_all']) for row in rows]
    reverse = [_num_or_none(row['pct_reverse_all']) for row in rows]
    undetermined = [_num_or_none(row['pct_undetermined_all']) for row in rows]
    x = list(range(len(thresholds)))
    fig_width = max(5.2, 0.34 * len(thresholds) + 2.2)
    fig, ax = plt.subplots(figsize=(fig_width, 3.3), dpi=300)
    series = [
        ('Forward strand', forward),
        ('Reverse strand', reverse),
        ('Undetermined', undetermined),
    ]
    for label, values in series:
        ax.plot(
            x,
            values,
            marker='o',
            markersize=4.4,
            linewidth=1.7,
            color=STRAND_COLORS[label],
            label=label,
        )

    ax.axvline(0, color='#D1D5DB', linewidth=0.8, linestyle='--', zorder=0)
    ax.set_ylim(0, 105)
    ax.set_ylabel('Assigned alignment records (%)', fontsize=8)
    ax.set_xlabel('Minimum aligned query length cutoff (nt)', fontsize=8)
    ax.set_title('Strand assignment by BAM aligned query length', fontsize=9.5, pad=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=5.8)
    ax.tick_params(axis='both', labelsize=6.5)
    ax.grid(axis='y', linestyle=':', linewidth=0.45, alpha=0.32)
    ax.set_axisbelow(True)

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.31), ncol=3, frameon=False, fontsize=6.5, handlelength=1.6)
    fig.text(
        0.5,
        0.012,
        'Records with aligned query length >= cutoff; 0 nt = no length filter.',
        ha='center',
        fontsize=5.7,
        color='#555555',
    )
    fig.subplots_adjust(left=0.11, right=0.98, top=0.84, bottom=0.34)

    out_pdf = os.path.join(output_dir, f'{base_name}_read_length_stratified_strandness_lineplot.pdf')
    out_png = os.path.join(output_dir, f'{base_name}_read_length_stratified_strandness_lineplot.png')
    fig.savefig(out_pdf, format='pdf')
    fig.savefig(out_png, format='png', dpi=300)
    plt.close(fig)
    print(f'Read-length stratified strandness line plot saved as: {out_pdf}, {out_png}')

def compute_splice_per_kb(bam_file, star_log=None, sample_n=500000):
    """Compute splices per kilobase (splice/kb), preferring STAR Log.final.out."""
    if star_log and os.path.isfile(star_log):
        try:
            metrics = {}
            with open(star_log) as f:
                for line in f:
                    if '|' in line:
                        k, v = [x.strip() for x in line.split('|')]
                        metrics[k] = v.replace('%', '').strip()
            n_junctions = int(metrics.get('Number of splices: Total', 0))
            n_unique = int(metrics.get('Uniquely mapped reads number', 0))
            avg_len = float(metrics.get('Average mapped length', 0))
            if n_unique > 0 and avg_len > 0:
                splice_per_kb = n_junctions * 1000.0 / (n_unique * avg_len)
                return round(splice_per_kb, 4), n_junctions, n_unique, round(avg_len, 1), f'STAR refined Log.final.out: {star_log}'
        except Exception as e:
            print(f'Warning: STAR log parse failed ({e}), falling back to BAM sampling.')

    try:
        cmd_total = ['samtools', 'view', '-c', '-F', '2304', '-q', '255', bam_file]
        r_total = subprocess.run(cmd_total, capture_output=True, text=True, check=True)
        n_total = int(r_total.stdout.strip())
        if n_total == 0:
            return 0.0, 0, 0, 0.0, 'BAM_alignment_record_estimate'

        frac = min(1.0, sample_n / n_total)
        sample_flag = ['-s', f'{frac:.6f}'] if frac < 1.0 else []
        cmd_view = ['samtools', 'view', '-F', '2304', '-q', '255'] + sample_flag + [bam_file]
        r_view = subprocess.run(cmd_view, capture_output=True, text=True, check=True)

        n_junctions_sampled = 0
        total_len = 0
        n_sampled = 0
        for line in r_view.stdout.splitlines():
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 10:
                continue
            n_sampled += 1
            n_junctions_sampled += count_cigar_junctions(fields[5])
            total_len += cigar_mapped_len(fields[5])

        if n_sampled == 0:
            return 0.0, 0, n_total, 0.0, 'BAM_alignment_record_estimate'

        avg_len = total_len / n_sampled
        junction_rate = n_junctions_sampled / n_sampled
        splice_per_kb = junction_rate * 1000.0 / avg_len if avg_len > 0 else 0.0
        n_junctions_est = int(junction_rate * n_total)
        return round(splice_per_kb, 4), n_junctions_est, n_total, round(avg_len, 1), 'BAM_alignment_record_estimate'

    except Exception as e:
        print(f'Warning: BAM-based splice/kb computation failed: {e}')
        return None


def run_infer_experiment(bam_file, refgene_bed, test_read_num, output_dir, star_log=None, mapq=30):
    """Run RSeQC infer_experiment, then generate Step 07 TSV and figures."""
    if not os.path.isfile(bam_file):
        print(f'Error: Input BAM file not found at {bam_file}')
        sys.exit(1)
    if not os.path.isfile(refgene_bed):
        print(f'Error: Reference BED file not found at {refgene_bed}')
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.basename(bam_file).replace('.bam', '')
    tsv_path = os.path.join(output_dir, f'{base_name}_bam2strandness.tsv')
    length_tsv_path = os.path.join(output_dir, f'{base_name}_read_length_stratified_strandness.tsv')

    frac_fwd, frac_rev, frac_failed, lib_type = '0', '0', '0', 'Unknown'
    spkb, n_junc, n_tot, avg_l, spkb_source = 'NA', 'NA', 'NA', 'NA', 'NA'

    with open(tsv_path, 'w') as f:
        f.write(
            "BAM\tdata type\t'1++,1--,2+-,2-+ (forward)'\t'1+-,1-+,2++,2-- (reverse)'\t"
            "%reads failed to determine:\t"
            "Splice_per_kb\tn_splice_junctions\tn_unique_reads\tavg_mapped_len_nt\tsplice_per_kb_source\n"
        )
        try:
            cmd = [
                'infer_experiment.py', '-i', bam_file, '-r', refgene_bed,
                '-s', str(test_read_num), '-q', str(mapq),
            ]
            proc = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            lines = proc.stdout.strip().split('\n')
            lib_type = lines[0].strip() if lines else 'Library type not found'
            frac_failed = lines[1].split(': ')[1] if len(lines) > 1 else '0'
            frac_fwd = lines[2].split(': ')[1] if len(lines) > 2 else '0'
            frac_rev = lines[3].split(': ')[1] if len(lines) > 3 else '0'

            result = compute_splice_per_kb(bam_file, star_log=star_log)
            if result:
                spkb, n_junc, n_tot, avg_l, spkb_source = result

            f.write(
                f'{bam_file}\t{lib_type}\t{frac_fwd}\t{frac_rev}\t{frac_failed}\t'
                f'{spkb}\t{n_junc}\t{n_tot}\t{avg_l}\t{spkb_source}\n'
            )
        except subprocess.CalledProcessError as e:
            error_message = f'Error processing {bam_file}: {e.stderr}'
            f.write(f'{bam_file}\tError\tNA\tNA\tNA\tNA\tNA\tNA\tNA\terror\n')
            print(error_message)
            sys.exit(1)

    print(f'Strandness results written to {tsv_path}')
    print(f'Splice/kb = {spkb}  (n_junctions={n_junc}, n_reads={n_tot}, avg_len={avg_l}nt, source={spkb_source})')

    # Legacy-compatible pure strandness pie plus independent splice-density plot.
    try:
        plot_strandness_pie(output_dir, base_name, lib_type, frac_fwd, frac_rev, frac_failed)
        plot_splice_junction_density(output_dir, base_name, spkb)
    except ValueError:
        print('Error: Could not convert strand fractions to numbers. Skipping strand/splice chart generation.')

    # New read-length-stratified strandness table and 100% stacked barplot.
    try:
        length_rows = compute_length_stratified_strandness(
            bam_file,
            refgene_bed,
            sample_size=test_read_num,
            min_mapq=mapq,
            legacy_fractions={
                'forward': frac_fwd,
                'reverse': frac_rev,
                'undetermined': frac_failed,
            },
        )
        write_length_stratified_tsv(length_rows, length_tsv_path)
        plot_length_stratified_strandness(length_rows, output_dir, base_name)
        print(f'Read-length stratified strandness TSV written to {length_tsv_path}')
    except Exception as e:
        print(f'Warning: read-length-stratified strandness generation failed: {e}', file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Infer strand specificity and generate a TSV summary and figures. '
            'Also computes splices per kilobase (splice/kb) as a complementary '
            'gDNA-contamination QC proxy.'
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument('--input_bam', required=True)
    parser.add_argument('--bed', required=True, help='Non-overlapping exon BED for infer_experiment.py.')
    parser.add_argument('--test_read_num', type=int, default=100000000)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--star_log', required=False, default=None, help='STAR Log.final.out from Step 6 refined alignment.')
    parser.add_argument('--mapq', type=int, default=30, help='Minimum MAPQ for strand inference; default matches RSeQC infer_experiment.py.')
    args = parser.parse_args()
    run_infer_experiment(
        args.input_bam,
        args.bed,
        args.test_read_num,
        args.output_dir,
        star_log=args.star_log,
        mapq=args.mapq,
    )


if __name__ == '__main__':
    main()
