#!/usr/bin/env bash
# Input bigWig file and set output directory.
bw_file="$1"
output_dir="$(dirname "$bw_file")"

# List of bed files.
bed_files=(
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_protein_coding.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_tRNAs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_rRNAs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_Y_RNAs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_snoRNAs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_snRNAs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_scaRNAs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_vault_RNAs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_lncRNAs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_TEC_protein_coding.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_IG_genes.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_TR_genes.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_misc-sncRNAs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_miRNAs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_pseudogenes.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_ERVs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_LINEs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_SINEs.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_piRNAs_gold_standard.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_RNAtype_bed/HG38_artifact.bed"
)

# Define output prefix based on bigWig file name.
output_prefix=$(basename "$bw_file" .bw)
num_colors=${#bed_files[@]}

# Generate color list from matplotlib tab10 colormap.
python_code=$(cat <<EOF
import matplotlib.pyplot as plt
import sys
num_colors = int(sys.argv[1])
palette = plt.get_cmap("tab10").colors
colors = []
for color in palette[:num_colors]:
    r, g, b = color[:3]
    colors.append("#%02x%02x%02x" % (int(r*255), int(g*255), int(b*255)))
print(" ".join(colors))
EOF
)
color_list=$(python -c "$python_code" "$num_colors")
IFS=' ' read -r -a colors <<< "$color_list"

# Check if computeMatrix output exists.
matrix_file="$output_dir/${output_prefix}_bed_stacked_matrix_RNA_types.gz"
sorted_regions_file="$output_dir/${output_prefix}_bed_stacked_sorted_regions_RNA_types.bed"

if [ ! -f "$matrix_file" ]; then
    echo "computeMatrix output not found, running computeMatrix..."
    computeMatrix scale-regions \
        -S "$bw_file" \
        -R "${bed_files[@]}" \
        --beforeRegionStartLength 0 \
        --regionBodyLength 1000 \
        --afterRegionStartLength 0 \
        -o "$matrix_file" \
        --binSize 100 \
        -p 20 \
        --blackListFileName /home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/Encode_hg38-blacklist.v2.bed \
        --outFileSortedRegions "$sorted_regions_file"
else
    echo "computeMatrix output exists, skipping computeMatrix..."
fi

# Create region labels with counts from sorted regions.
labels=("protein_coding" "tRNAs" "rRNAs" "Y_RNAs" "snoRNAs" "snRNAs" "scaRNAs" "vault_RNAs" "lncRNAs" "TEC_protein_coding" "IG_genes" "TR_genes" "misc-sncRNAs" "miRNAs" "pseudogenes" "ERVs" "LINEs" "SINEs" "piRNAs" "artifact")
regions_label=()
for i in "${!bed_files[@]}"; do
    count=$(grep -c "$(basename "${bed_files[i]}")" "$sorted_regions_file")
    regions_label+=("${labels[i]} (${count})")
done

# Run plotProfile to generate SVG output.
echo "Running plotProfile for SVG..."
plotProfile -m "$matrix_file" \
    -out "$output_dir/${output_prefix}_bed_stacked_profile_RNA_types.svg" \
    --colors "${colors[@]}" \
    --legendLocation upper-left \
    --startLabel "Start (5')" \
    --endLabel "End (3')" \
    --plotWidth 12 --plotHeight 15 --dpi 300 \
    --yAxisLabel "Signal intensity (BPM) over scaled regions" \
    --yMin 0 --yMax 10 \
    --plotType lines \
    --numPlotsPerRow 1 \
    --regionsLabel "${regions_label[@]}"

echo "SVG profile generated at ${output_prefix}_bed_stacked_profile_RNA_types.svg"

# Run plotProfile to generate PNG output.
echo "Running plotProfile for PNG..."
plotProfile -m "$matrix_file" \
    -out "$output_dir/${output_prefix}_bed_stacked_profile_RNA_types.png" \
    --colors "${colors[@]}" \
    --legendLocation upper-left \
    --startLabel "Start (5')" \
    --endLabel "End (3')" \
    --plotWidth 12 --plotHeight 15 --dpi 300 \
    --yAxisLabel "Signal intensity (BPM) over scaled regions" \
    --yMin 0 --yMax 10 \
    --plotType lines \
    --numPlotsPerRow 1 \
    --regionsLabel "${regions_label[@]}"

echo "PNG profile generated at ${output_prefix}_bed_stacked_profile_RNA_types.png"

