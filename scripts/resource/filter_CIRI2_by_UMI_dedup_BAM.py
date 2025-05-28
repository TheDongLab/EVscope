#!/usr/bin/env python3
"""
#usage: python /home/yz2474/scripts/donglab/EV-RNA-Profiler/filter_CIRI2_by_UMI_dedup_BAM.py  --input_bam /home/yz2474/yiyong_2023/Research/EV-RNA-seq/R01_PD_IMAP/Ewa/run_SMARTer_pipe/Ewa_DNAnexus_PD_SNCA_d42d44/07_UMI_dedup_bam/Ewa_DNAnexus_PD_SNCA_d42d44_UMI_dedup.sorted.bam --CIRI2_result Ewa_DNAnexus_PD_SNCA_d42d44_CIRI2_out.tsv --output CIRI2_dedup_out.tsv

This script refines CIRI2 output by filtering the circular junction read IDs 
using a UMI-deduplicated BAM file. It updates the junction read count and adds 
a CPM column based on the total number of unique fragments from the BAM.
Output columns (tab-delimited):
  1. circRNA_ID (merged as "originalID|gene_id" or "originalID|novel" if gene_id is missing or "n/a")
  2. #junction_reads (updated count)
  3. junction_reads_CPM (calculated as (updated_count/total_unique_fragments)*1e6)
  4. circRNA_type
  5. strand
  6. junction_reads_ID (filtered, comma-separated)
Rows with an updated junction read count of zero are omitted.
"""

import sys
import argparse
import pysam

def parse_args():
    parser = argparse.ArgumentParser(
        description="Refine CIRI2 output using a UMI-deduplicated BAM file."
    )
    parser.add_argument("--CIRI2_result", dest="ciri2_file", required=True,
                        help="Path to the CIRI2 result TSV file")
    parser.add_argument("--input_bam", dest="bam_file", required=True,
                        help="Path to the UMI-deduplicated sorted BAM file")
    parser.add_argument("--output", required=True,
                        help="Path for the output deduplicated TSV file")
    return parser.parse_args()

def load_bam_read_ids(bam_file):
    """
    Load unique read names (query names) from the BAM file into a set.
    For paired-end data, each unique query name represents one fragment.
    """
    read_ids = set()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for alignment in bam.fetch(until_eof=True):
            read_ids.add(alignment.query_name)
    return read_ids

def dedup_ciri2_results(ciri2_file, bam_read_ids, total_unique_fragments, output_file):
    """
    Reads the CIRI2 TSV file, filters the junction read IDs (column 12) to retain only 
    those present in bam_read_ids, updates the junction read count, computes junction_reads_CPM, 
    and writes a new TSV file with selected columns.
    
    Output columns:
      1. circRNA_ID (merged as "originalID|gene_id" or "originalID|novel" if gene_id is missing or "n/a")
      2. junction_reads (updated count)
      3. junction_reads_CPM (calculated as (updated_count/total_unique_fragments)*1e6)
      4. circRNA_type
      5. strand
      6. junction_reads_ID (filtered, comma-separated)
    """
    with open(ciri2_file, "r") as fin, open(output_file, "w") as fout:
        header = fin.readline().strip()
        # Write new header (omit separate gene_id column as it is merged)
        fout.write("\t".join(["circRNA_ID", "junction_reads", "junction_reads_CPM", 
                              "circRNA_type", "strand", "junction_reads_ID"]) + "\n")
        
        for line in fin:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            # Fields:
            # 0: circRNA_ID
            # 8: circRNA_type
            # 9: gene_id (may have trailing comma or "n/a")
            # 10: strand
            # 11: junction_reads_ID (comma-separated)
            original_id = fields[0]
            gene_id_raw = fields[9].strip()
            gene_id_val = gene_id_raw.rstrip(",")
            # Replace empty or "n/a" (case-insensitive) with "novel"
            if not gene_id_val or gene_id_val.lower() == "n/a":
                gene_id_val = "novel"
            # Merge original circRNA_ID with gene_id_val
            merged_id = f"{original_id}|{gene_id_val}"
            
            circRNA_type = fields[8]
            strand = fields[10]
            junction_reads_field = fields[11].strip()
            if junction_reads_field:
                read_ids_list = junction_reads_field.split(",")
            else:
                read_ids_list = []
            # Filter read IDs that are present in the BAM file
            filtered_ids = [rid for rid in read_ids_list if rid in bam_read_ids]
            updated_count = len(filtered_ids)
            # Skip rows with zero updated junction reads
            if updated_count == 0:
                continue
            # Calculate junction_reads_CPM = (updated_count / total_unique_fragments) * 1e6
            junction_reads_CPM = (updated_count / total_unique_fragments) * 1e6 if total_unique_fragments > 0 else 0
            filtered_ids_str = ",".join(filtered_ids)
            out_fields = [merged_id, str(updated_count), f"{junction_reads_CPM:.2f}",
                          circRNA_type, strand, filtered_ids_str]
            fout.write("\t".join(out_fields) + "\n")

def main():
    args = parse_args()
    bam_read_ids = load_bam_read_ids(args.bam_file)
    total_unique_fragments = len(bam_read_ids)
    
    # Use the provided output file path.
    output_file = args.output
    
    dedup_ciri2_results(args.ciri2_file, bam_read_ids, total_unique_fragments, output_file)
    print(f"Deduplicated CIRI2 results written to: {output_file}")

if __name__ == "__main__":
    main()

