#!/usr/bin/env bash
# run_deeptools_plotprfile_bed_stackd_202412_ZYY.sh
# This script generates a stacked profile plot from a bigWig file using computeMatrix and plotProfile.
# The Y-axis scale is fixed (--yMin 0 --yMax 1) so that the scale remains consistent across samples.
# The output plot is saved as ${output_prefix}_bed_stacked_profile.pdf.

bw_file="$1"
output_dir="$(dirname "$bw_file")"
bed_files=(
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_promoter_1500_500bp_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_5UTR_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_exon_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intron_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_3UTR_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_downstream_2kb_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intergenic_noOverlap.bed"
)
output_prefix=$(basename "$bw_file" .bw)
num_colors=${#bed_files[@]}

# Get a color list from tab10.
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

# Run computeMatrix to generate the matrix and sorted regions.
echo "Running computeMatrix..."
#computeMatrix scale-regions \
    -S "$bw_file" \
    -R "${bed_files[@]}" \
    --beforeRegionStartLength 0 \
    --regionBodyLength 1000 \
    --afterRegionStartLength 0 \
    -o "$output_dir/${output_prefix}_bed_stacked_matrix_meta_gene.gz" \
    -p 20 \
    --blackListFileName /home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/Encode_hg38-blacklist.v2.bed \
    --outFileSortedRegions "$output_dir/${output_prefix}_bed_stacked_sorted_regions_meta_gene.bed"

# Create labels for each region based on sorted regions.
labels=("Promoters" "5'UTR" "Exon" "Intron" "3'UTR" "Downstream" "Intergenic")
regions_label=()
for i in "${!bed_files[@]}"; do
    count=$(grep -c "$(basename "${bed_files[i]}")" "$output_dir/${output_prefix}_bed_stacked_sorted_regions_meta_gene.bed")
    regions_label+=("${labels[i]} (${count})")
done

echo "Running plotProfile..."
plotProfile -m "$output_dir/${output_prefix}_bed_stacked_matrix_meta_gene.gz" \
    -out "$output_dir/${output_prefix}_bed_stacked_profile_meta_gene.png" \
    --colors "${colors[@]}" \
    --legendLocation upper-left \
    --startLabel "Start (5')" \
    --endLabel "End (3')" \
    --plotWidth 12 --plotHeight 15 --dpi 300 \
    --yAxisLabel "Signal intensity (BPM) over scaled regions" \
    --yMin 0 --yMax 5 \
    --plotType lines \
    --numPlotsPerRow 1 \
    --regionsLabel "${regions_label[@]}"

plotProfile -m "$output_dir/${output_prefix}_bed_stacked_matrix_meta_gene.gz" \
    -out "$output_dir/${output_prefix}_bed_stacked_profile_meta_gene.svg" \
    --colors "${colors[@]}" \
    --legendLocation upper-left \
    --startLabel "Start (5')" \
    --endLabel "End (3')" \
    --plotWidth 12 --plotHeight 15 --dpi 300 \
    --yAxisLabel "Signal intensity (BPM) over scaled regions" \
    --yMin 0 --yMax 5 \
    --plotType lines \
    --numPlotsPerRow 1 \
    --regionsLabel "${regions_label[@]}"
echo "Plot profile generated at ${output_prefix}_bed_stacked_profile_meta_gene.pdf"

