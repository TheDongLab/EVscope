#!/usr/bin/env bash
bw_file="$1"
output_dir="$(dirname "$bw_file")"
bed_files=(
    "/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_promoter_1500_500bp_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_5UTR_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_exon_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intron_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_3UTR_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_downstream_2kb_noOverlap.bed"
    "/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intergenic_noOverlap.bed"
)
output_prefix=$(basename "$bw_file" .bw)
num_colors=${#bed_files[@]}
python_code=$(cat <<EOF
import matplotlib.pyplot as plt
import sys

# 从命令行读取 num_colors 变量
num_colors = int(sys.argv[1])
palette = plt.colormaps["tab10"]  # 使用新的方式访问 colormap
colors = []
for color in palette.colors[:num_colors]:
    r, g, b = color[:3]  # 确保只解包 RGB 部分
    colors.append("#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255)))
print(" ".join(colors))
EOF
)
color_list=$(python -c "$python_code" "$num_colors")
IFS=' ' read -r -a colors <<< "$color_list"
# 运行 computeMatrix 生成矩阵文件和排序后的区域文件
echo "Running computeMatrix..."
computeMatrix scale-regions \
    -S "$bw_file" \
    -R "${bed_files[@]}" \
    --beforeRegionStartLength 1000 \
    --regionBodyLength 1000 \
    --afterRegionStartLength 1000 \
    -o "$output_dir/${output_prefix}_bed_stacked_matrix.gz" \
    -p 20 \
    --blackListFileName /home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/Encode_Blacklist/lists/blacklist.v2.bed \
    --outFileSortedRegions "$output_dir/${output_prefix}_bed_stacked_sorted_regions.bed"
# 统计每个 BED 文件的有效区域数量
labels=("Promoters" "5'UTR" "Exon" "Intron" "3'UTR" "Downstream" "Intergenic")
regions_label=()
for i in "${!bed_files[@]}"; do
    # 使用 grep 统计排序文件中包含该 BED 文件名称的行数
    count=$(grep -c "$(basename "${bed_files[i]}")" "$output_dir/${output_prefix}_bed_stacked_sorted_regions.bed")
    regions_label+=("${labels[i]} (${count})")
done
echo "Running plotProfile..."
plotProfile -m "$output_dir/${output_prefix}_bed_stacked_matrix.gz" \
    -out "$output_dir/${output_prefix}_bed_stacked_density_profile.pdf" \
    --colors "${colors[@]}" \
    --legendLocation upper-left \
    --startLabel "Start (5')" \
    --endLabel "End (3')" \
    --plotWidth 12 --plotHeight 15 --dpi 300 \
    --yAxisLabel "Signal intensity (BPM) over scaled regions" \
    --plotType lines \
    --numPlotsPerRow 1 \
    --regionsLabel "${regions_label[@]}"

echo "Plot profile generated at ${output_prefix}_bed_stacked_profile.pdf"

