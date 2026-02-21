#!/usr/bin/env python3
"""
Usage: python 04_convert_CIRI2CPM.py --CIRI2_result <CIRI2_result.tsv> --input_sam <UMI_dedup.sorted.sam> --output <CIRI2_dedup_out.tsv> --GeneID_meta_table <geneID_meta_table.tsv>

This script refines CIRI2 output using a UMI-deduplicated SAM file to calculate CPM.
It outputs 6 columns:
  1. circRNA_ID1: original_id|gene_id (e.g. chr15:27326808|27328888|ENSG00000182256.13)
  2. circRNA_ID2: original_id|gene_symbol (from meta table, or 'novel' if missing)
  3. junction_read_counts: count of junction reads
  4. CPM: (junction_read_counts / total_unique_fragments) * 1e6
  5. circRNA_type: as in CIRI2 output
  6. Strandness: strand information
Rows with junction_read_counts equal to zero are omitted.
"""

import argparse
import pysam

def parse_args():
    parser = argparse.ArgumentParser(
        description="Refine CIRI2 output using a UMI-deduplicated SAM file to calculate CPM."
    )
    parser.add_argument("--CIRI2_result", dest="ciri2_file", required=True,
                        help="Path to the CIRI2 result TSV file")
    parser.add_argument("--input_sam", dest="sam_file", required=True,
                        help="Path to the UMI-deduplicated sorted SAM file")
    parser.add_argument("--output", required=True,
                        help="Path for the output TSV file")
    parser.add_argument("--GeneID_meta_table", dest="meta_table", required=True,
                        help="Path to the GeneID meta table TSV file")
    return parser.parse_args()

def load_sam_unique_read_ids(sam_file):
    """Load unique read names from the SAM file."""
    unique_reads = set()
    with pysam.AlignmentFile(sam_file, "r") as sam:
        for alignment in sam.fetch(until_eof=True):
            if not alignment.is_unmapped:
                unique_reads.add(alignment.query_name)
    return unique_reads

def load_gene_meta_table(meta_table_file):
    """Load gene meta table and return a dict mapping GeneID to GeneSymbol."""
    gene_meta = {}
    with open(meta_table_file, "r") as fin:
        header = fin.readline().strip().split("\t")
        # Assume first column is GeneID and second column is GeneSymbol
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            gene_id = parts[0].strip()
            gene_symbol = parts[1].strip()
            gene_meta[gene_id] = gene_symbol
    return gene_meta

def update_ciri2_results(ciri2_file, mapped_unique_fragments, output_file, gene_meta):
    """Update CIRI2 results and write output TSV with 6 columns."""
    with open(ciri2_file, "r") as fin, open(output_file, "w") as fout:
        # Write header
        fout.write("\t".join(["circRNA_ID1", "circRNA_ID2", "junction_read_counts", "CPM", "circRNA_type", "Strandness"]) + "\n")
        # Skip header line of the CIRI2 file
        header = fin.readline()
        for line in fin:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            # Get original id from column 0 (e.g. "chr15:27326808|27328888")
            chrID=fields[1]
            start_posi=int(fields[2])-1
            end_posi=fields[3]
            original_id=chrID+":"+str(start_posi)+"|"+str(end_posi)
            # Get circRNA type from column 8
            circRNA_type = fields[8]
            # Get gene id from column 9 and remove trailing comma
            gene_id_raw = fields[9].strip()
            gene_id_val = gene_id_raw.rstrip(",")
            if not gene_id_val or gene_id_val.lower() == "n/a":
                gene_id_val = "novel"
                gene_symbol = "novel"
            else:
                gene_symbol = gene_meta.get(gene_id_val, "novel")
            # Get strand from column 10
            strand = fields[10]
            # Get junction read IDs from column 11 and split by comma
            junction_reads_field = fields[11].strip()
            if junction_reads_field:
                read_ids_list = [rid for rid in junction_reads_field.split(",") if rid]
            else:
                read_ids_list = []
            updated_count = len(read_ids_list)
            if updated_count == 0:
                continue
            # Compute CPM
            junction_reads_CPM = (updated_count / mapped_unique_fragments) * 1e6 if mapped_unique_fragments > 0 else 0
            # Build circRNA_ID1 and circRNA_ID2
            circRNA_ID1 = f"{original_id}|{gene_id_val}"
            circRNA_ID2 = f"{original_id}|{gene_symbol}"
            # Write output row
            fout.write("\t".join([circRNA_ID1, circRNA_ID2, str(updated_count),
                                  f"{junction_reads_CPM:.2f}", circRNA_type, strand]) + "\n")

def main():
    args = parse_args()
    # Load unique read IDs from the SAM file
    unique_read_ids = load_sam_unique_read_ids(args.sam_file)
    mapped_unique_fragments = len(unique_read_ids)
    # Load gene meta table mapping
    gene_meta = load_gene_meta_table(args.meta_table)
    # Update CIRI2 results and write output
    update_ciri2_results(args.ciri2_file, mapped_unique_fragments, args.output, gene_meta)
    print(f"Updated CIRI2 results written to: {args.output}")

if __name__ == "__main__":
    main()

