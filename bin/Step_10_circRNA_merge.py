#!/usr/bin/env python3
import argparse
import csv
import os
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Combine circRNA outputs from CIRCexplorer2 and CIRI2 and generate Venn diagrams."
    )
    parser.add_argument("--CIRCexplorer2", required=True, help="Input TSV from CIRCexplorer2")
    parser.add_argument("--CIRI2", required=True, help="Input TSV from CIRI2")
    parser.add_argument("--output_matrix", required=True, help="Output combined TSV file")
    parser.add_argument("--out_venn", required=True, help="Output PDF file for Venn diagram")
    return parser.parse_args()

def parse_circRNA_id1(cid):
    # Parse circRNA_ID1: expected format 'chr15:27326808|27328888|...'
    chrom, rest = cid.split(":", 1)
    parts = rest.split("|")
    start = int(parts[0])
    end = int(parts[1])
    return chrom, start, end

def extract_gene_id(cid2):
    # Extract geneID from circRNA_ID2; return None if "novel"
    parts = cid2.split("|")
    if len(parts) >= 3:
        gene = parts[2]
        if gene.lower() == "novel":
            return None
        else:
            return gene
    return None

def get_coordinate(cid2):
    # Return coordinate portion (first two fields) of circRNA_ID2
    parts = cid2.split("|")
    if len(parts) >= 2:
        return "|".join(parts[:2])
    return cid2

def read_file(filename, method):
    # Read TSV file and add extra fields.
    rows = []
    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            row['method'] = method
            row['junction_read_counts'] = int(row['junction_read_counts']) if row['junction_read_counts'] else 0
            row['CPM'] = float(row['CPM']) if row['CPM'] else 0.0
            chrom, start, end = parse_circRNA_id1(row['circRNA_ID1'])
            row['chrom'] = chrom
            row['start'] = start
            row['end'] = end
            row['gene_id'] = extract_gene_id(row['circRNA_ID2'])
            row['coordinate'] = get_coordinate(row['circRNA_ID2'])
            rows.append(row)
    return rows

def merge_entries(rows):
    # Merge circRNAs based on criteria.
    groups = defaultdict(list)
    for row in rows:
        key = (row['chrom'], row['end'], row['Strandness'])
        groups[key].append(row)

    merged = []
    for key, group in groups.items():
        group.sort(key=lambda r: r['start'])
        clusters = []
        for r in group:
            if not clusters:
                clusters.append([r])
            else:
                current_cluster = clusters[-1]
                if abs(r['start'] - current_cluster[0]['start']) < 2:
                    current_cluster.append(r)
                else:
                    clusters.append([r])
        for cluster in clusters:
            rep = max(cluster, key=lambda r: r['junction_read_counts'])
            methods = {r['method'] for r in cluster}
            if len(methods) > 1:
                rep['Source'] = "CIRI2 (overlapped)" if rep['method'] == "CIRI2" else "CIRCexplorer2 (overlapped)"
            else:
                rep['Source'] = rep['method']
            merged.append(rep)
    return merged

def assign_circRNA_ID3(rows):
    # Assign new circRNA_ID3 names.
    groups = defaultdict(list)
    for row in rows:
        key = row['gene_id'] if row['gene_id'] is not None else "de_novo"
        groups[key].append(row)

    for gene, group in groups.items():
        if gene == "de_novo":
            for r in group:
                r['circRNA_ID3'] = r['coordinate']
        else:
            group.sort(key=lambda r: r['start'])
            for idx, r in enumerate(group, start=1):
                r['circRNA_ID3'] = f"circRNA{idx}|{gene}"
    return rows

def write_output(rows, output_file):
    # Write output TSV with selected columns.
    fieldnames = ['circRNA_ID1', 'circRNA_ID2', 'circRNA_ID3', 'Source',
                  'junction_read_counts', 'CPM', 'circRNA_type', 'Strandness']
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for r in rows:
            writer.writerow({
                'circRNA_ID1': r['circRNA_ID1'],
                'circRNA_ID2': r['circRNA_ID2'],
                'circRNA_ID3': r.get('circRNA_ID3', ''),
                'Source': r['Source'],
                'junction_read_counts': r['junction_read_counts'],
                'CPM': r['CPM'],
                'circRNA_type': r['circRNA_type'],
                'Strandness': r['Strandness']
            })

def plot_venn_id1(ce2_rows, c2_rows, output_path):
    # Generate Venn diagram and save as PDF.
    set_ce2 = set(r["circRNA_ID1"] for r in ce2_rows)
    set_c2 = set(r["circRNA_ID1"] for r in c2_rows)

    plt.figure(figsize=(6,6))
    venn2(subsets=(len(set_ce2 - set_c2), len(set_c2 - set_ce2), len(set_ce2 & set_c2)),
          set_labels=("CIRCexplorer2", "CIRI2"))
    plt.title("Venn diagram of circRNAs identified between CIRCexplorer2 and CIRI2")
    # Save as both PNG and PDF
    base, ext = os.path.splitext(output_path)
    png_path = base + ".png" if ext.lower() != ".png" else output_path
    pdf_path = base + ".pdf" if ext.lower() != ".pdf" else output_path
    plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.close()

def main():
    args = parse_args()
    # Read input files
    ce2_rows = read_file(args.CIRCexplorer2, "CIRCexplorer2")
    c2_rows = read_file(args.CIRI2, "CIRI2")
    all_rows = ce2_rows + c2_rows

    # Merge entries and assign new IDs
    merged_rows = merge_entries(all_rows)
    final_rows = assign_circRNA_ID3(merged_rows)
    
    # Write combined TSV output
    write_output(final_rows, args.output_matrix)
    
    # Plot Venn diagram and save as specified PDF
    plot_venn_id1(ce2_rows, c2_rows, args.out_venn)
    
    print(f"Combined TSV written to {args.output_matrix}")
    print(f"Venn diagram saved as {args.out_venn}")

if __name__ == "__main__":
    main()

