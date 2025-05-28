#!/usr/bin/env python3
import argparse
import csv

def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Convert featureCounts output to a TPM expression matrix."
    )
    # Note: using --featureCounts_out (spelling corrected)
    parser.add_argument("--featureCounts_out", required=True,
                        help="Path to the featureCounts output file")
    parser.add_argument("--GeneID_meta_table", required=True,
                        help="Path to the gene meta table (mapping GeneID to GeneSymbol and GeneType)")
    parser.add_argument("--output", required=True,
                        help="Output file path for the TPM expression matrix")
    return parser.parse_args()

def read_meta_table(meta_file):
    """
    Reads the gene meta table and returns a dictionary mapping:
       GeneID -> (GeneSymbol, GeneType)
    """
    meta_dict = {}
    with open(meta_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)  # skip header
        for row in reader:
            if len(row) < 3:
                continue
            gene_id = row[0].strip()
            gene_symbol = row[1].strip() if row[1].strip() else "NA"
            gene_type = row[2].strip() if row[2].strip() else "NA"
            meta_dict[gene_id] = (gene_symbol, gene_type)
    return meta_dict

def read_featureCounts(featurecounts_file):
    """
    Reads the featureCounts output file.
    Skips comment lines (starting with '#') and the header line.
    For each gene, extracts GeneID (column 1), gene length (column 6),
    and read counts (column 7). Genes with zero length or zero read counts
    are skipped.
    Calculates RPK using:
         RPK = read_count / (length/1000) = read_count * 1000 / length
    Returns:
         genes_data: list of dictionaries with keys: GeneID, ReadCounts, Length, RPK
         total_rpk: sum of all RPK values
    """
    genes_data = []
    total_rpk = 0.0
    header_found = False
    with open(featurecounts_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Skip comment lines starting with '#'
            if line.startswith("#"):
                continue
            # The first non-comment line is the header; skip it.
            if not header_found:
                header_found = True
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            gene_id = parts[0].strip()
            try:
                length = float(parts[5])
                read_count = float(parts[6])
            except ValueError:
                continue
            # Skip genes with zero length or zero read counts
            if length <= 0 or read_count <= 0:
                continue
            # Calculate RPK (reads per kilobase)
            rpk = read_count * 1000.0 / length
            total_rpk += rpk
            genes_data.append({
                "GeneID": gene_id,
                "ReadCounts": read_count,
                "Length": length,
                "RPK": rpk
            })
    return genes_data, total_rpk

def compute_tpm(genes_data, total_rpk):
    """
    Computes TPM for each gene using:
         TPM = (gene_RPK / total_RPK) * 1e6
    Adds the TPM value to each gene dictionary.
    """
    for gene in genes_data:
        gene["TPM"] = gene["RPK"] / total_rpk * 1e6 if total_rpk > 0 else 0
    return genes_data

def write_output(genes_data, meta_dict, output_file):
    """
    Writes the final TPM expression matrix to the output file.
    The output file contains the following columns:
         GeneID, GeneSymbol, GeneType, ReadCounts, TPM
    If a GeneID is not found in the meta table, 'NA' is used for GeneSymbol and GeneType.
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["GeneID", "GeneSymbol", "GeneType", "ReadCounts", "TPM"])
        for gene in genes_data:
            gene_id = gene["GeneID"]
            read_counts = gene["ReadCounts"]
            tpm = gene["TPM"]
            gene_symbol, gene_type = meta_dict.get(gene_id, ("NA", "NA"))
            writer.writerow([gene_id, gene_symbol, gene_type, read_counts, tpm])

def main():
    args = parse_args()
    # Read gene meta table for annotation mapping.
    meta_dict = read_meta_table(args.GeneID_meta_table)
    # Read featureCounts output to get read counts and gene lengths.
    genes_data, total_rpk = read_featureCounts(args.featureCounts_out)
    if total_rpk <= 0:
        print("Warning: total RPK is zero. No expressed genes found.")
    # Compute TPM values.
    genes_data = compute_tpm(genes_data, total_rpk)
    # Write the TPM expression matrix to the output file.
    write_output(genes_data, meta_dict, args.output)

if __name__ == "__main__":
    main()

