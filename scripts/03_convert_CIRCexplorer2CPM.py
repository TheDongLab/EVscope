#!/usr/bin/env python3
"""
Usage: python 03_convert_CIRCexplorer2CPM.py --CIRCexplorer2_result <result.txt> --input_bam <dedup.bam> --output <out.tsv> --GeneID_meta_table <meta_table.tsv>

This script refines CIRCexplorer2 output by updating junction read counts and calculating CPM 
based on unique fragments from the input BAM. It also maps GeneID to gene symbol using a provided meta table.
Output columns (tab-delimited):
  1. circRNA_ID1: "chrom:start|end|GeneID"
  2. circRNA_ID2: "chrom:start|end|GeneSymbol"
  3. junction_read_counts
  4. CPM
  5. circRNA_type
  6. Strandness
Rows with zero junction read counts are omitted.
"""

import argparse
import pysam

def parse_args():
    parser = argparse.ArgumentParser(
        description="Refine CIRCexplorer2 output and calculate CPM with gene symbol mapping."
    )
    parser.add_argument("--CIRCexplorer2_result", required=True,
                        help="Path to the CIRCexplorer2 result file (TSV, no header)")
    parser.add_argument("--input_bam", required=True,
                        help="Path to deduplicated sorted BAM file")
    parser.add_argument("--output", required=True,
                        help="Path for the output TSV file")
    parser.add_argument("--GeneID_meta_table", required=True,
                        help="Path to gene meta table (TSV with header: GeneID, GeneSymbol, GeneType)")
    return parser.parse_args()

def load_bam_unique_read_ids(bam_file):
    """Load unique query names from BAM; each unique name represents one fragment."""
    unique_reads = set()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if not aln.is_unmapped:
                unique_reads.add(aln.query_name)
    return unique_reads

def load_gene_meta_table(meta_file):
    """Load gene meta table; return dict mapping GeneID to GeneSymbol."""
    gene_meta = {}
    with open(meta_file, "r") as fin:
        header = fin.readline()  # skip header
        for line in fin:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene_id = parts[0].strip()
            gene_symbol = parts[1].strip()
            gene_meta[gene_id] = gene_symbol
    return gene_meta

def update_circexplorer2_results(result_file, unique_fragments, gene_meta, output_file):
    """
    Process the CIRCexplorer2 file and write output with 6 columns:
      circRNA_ID1, circRNA_ID2, junction_read_counts, CPM, circRNA_type, Strandness.
    circRNA_ID1: "chrom:start|end|GeneID"
    circRNA_ID2: "chrom:start|end|GeneSymbol"
    """
    with open(result_file, "r") as fin, open(output_file, "w") as fout:
        # Write header
        fout.write("\t".join(["circRNA_ID1", "circRNA_ID2", "junction_read_counts", "CPM", "circRNA_type", "Strandness"]) + "\n")
        for line in fin:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            # Expecting at least 18 fields
            if len(fields) < 18:
                continue
            chrom = fields[0]
            start = fields[1]
            end = fields[2]
            strand = fields[5]
            try:
                junction_reads = int(fields[12])
            except ValueError:
                continue
            if junction_reads == 0:
                continue
            circRNA_type = fields[13]
            gene_id = fields[14].strip()
            if not gene_id or gene_id.lower() == "n/a":
                gene_id = "novel"
            # Lookup gene symbol using meta table; default to "novel" if not found.
            gene_symbol = gene_meta.get(gene_id, "novel")
            # Build circRNA_ID1 and circRNA_ID2
            circRNA_ID1 = f"{chrom}:{start}|{end}|{gene_id}"
            circRNA_ID2 = f"{chrom}:{start}|{end}|{gene_symbol}"
            # Calculate CPM
            cpm = (junction_reads / unique_fragments) * 1e6 if unique_fragments > 0 else 0
            out_fields = [circRNA_ID1, circRNA_ID2, str(junction_reads), f"{cpm:.2f}", circRNA_type, strand]
            fout.write("\t".join(out_fields) + "\n")

def main():
    args = parse_args()
    unique_reads = load_bam_unique_read_ids(args.input_bam)
    unique_fragments = len(unique_reads)
    gene_meta = load_gene_meta_table(args.GeneID_meta_table)
    update_circexplorer2_results(args.CIRCexplorer2_result, unique_fragments, gene_meta, args.output)
    print(f"Updated CIRCexplorer2 results written to: {args.output}")

if __name__ == "__main__":
    main()

