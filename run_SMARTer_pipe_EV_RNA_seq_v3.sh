#!/usr/bin/env bash
#Usage: bash run_SMARTer_pipe_v3_ZYY.sh B4_TDP43_KD_caRNA_rep3 /path/to/R1.fq.gz /path/to/R2.fq.gz

#######################################
# User-defined Variables
#######################################
set -eo pipefail
Header=$1
R1_path=$2
R2_path=$3
outputDir="${Header}_processed"

genomeDir="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_STAR-2.7.10b_sjdbOverhang100_Gene"
Gene_GTF="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/3_gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.gtf"
refflat="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/4_gencode.v45.chr_patch_hapl_scaff.annotation_UCSC_RepeatMasker_v4.0.7_Dfam2.0_rmskJoinedCurrent.bed12_converted.LINE.refflat"
genome_fa="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Sequence/hg38p14/hg38.p14.whole.genome.fa"
rRNA_fa="/home/yz2474/yiyong_2023/soft/sortmerna/index_build/smr_v4.3_default_db.fasta"

# function: run a step only if 'step.done' doesn't exist
run_step() {
    local step_dir=$1
    echo "Running in: ${step_dir}"
    local step_cmds=$2
    mkdir -p "$step_dir"
    if [ ! -f "$step_dir/step.done" ]; then
        eval "$step_cmds"
        touch "$step_dir/step.done"
    fi
}

#######################################
# Main Pipeline
#######################################
mkdir -p "${outputDir}"

# Step 1: raw data QC + 14-bp motif for R2
run_step "${outputDir}/1_raw_QC_and_motif" "
    # 1. QC
    fastqc -o ${outputDir}/1_raw_QC_and_motif -t 20 $R1_path $R2_path
    multiqc ${outputDir}/1_raw_QC_and_motif -o ${outputDir}/1_raw_QC_and_motif

    # 2. Generate a small 'list' file for the motif script (only R2)
    echo -e '${R2_path}\t${Header}' > ${outputDir}/1_raw_QC_and_motif/motif_fastq.list

    # 3. Run get_SMARTer_Read2_14BP_motif.py
    #    We assume you have it in PATH or specify the full path below:
    python /home/yz2474/scripts/donglab/get_SMARTer_Read2_14BP_motif/get_SMARTer_Read2_14BP_motif.py \
      ${outputDir}/1_raw_QC_and_motif/motif_fastq.list \
      -n 14 -r 100000 \
      -o ${outputDir}/1_raw_QC_and_motif/motif_output \
      -p 10
"

# Step 2: UMI extraction and adapter trimming
run_step "${outputDir}/2_UMI_trim" "
    umi_tools extract --bc-pattern='NNNNNNNN' \
      --stdin=$R2_path --stdout=${outputDir}/2_UMI_trim/R2_umi.fq.gz \
      --read2-in=$R1_path --read2-out=${outputDir}/2_UMI_trim/R1_umi.fq.gz \
      --log=${outputDir}/2_UMI_trim/umi_extract.log --umi-separator='_'

    cutadapt -u 14 -o ${outputDir}/2_UMI_trim/R2_umi_trim.fq.gz ${outputDir}/2_UMI_trim/R2_umi.fq.gz

    trim_galore --phred33 --fastqc --illumina --gzip --stringency 3 --cores 20 --paired --length 10 \
      ${outputDir}/2_UMI_trim/R1_umi.fq.gz ${outputDir}/2_UMI_trim/R2_umi_trim.fq.gz \
      -o ${outputDir}/2_UMI_trim

    fastqc -o ${outputDir}/2_UMI_trim -t 20 \
      ${outputDir}/2_UMI_trim/R1_umi_val_1.fq.gz ${outputDir}/2_UMI_trim/R2_umi_trim_val_2.fq.gz

    multiqc ${outputDir}/2_UMI_trim -o ${outputDir}/2_UMI_trim
"

# Step 3: rRNA detection from trimmed data
run_step "${outputDir}/3_ribodetector" "
    ribodetector_cpu -t 20 \
      -l 100 \
      -i ${outputDir}/2_UMI_trim/R1_umi_val_1.fq.gz ${outputDir}/2_UMI_trim/R2_umi_trim_val_2.fq.gz \
      -e rrna \
      --chunk_size 256 \
      -r ${outputDir}/3_ribodetector/R1_umi_rRNA.fq.gz ${outputDir}/3_ribodetector/R2_umi_rRNA.fq.gz \
      -o /dev/null /dev/null
"

# Step 4: Kraken downsampling and annotation
run_step "${outputDir}/4_kraken" "
    seqtk sample -s100 ${outputDir}/2_UMI_trim/R1_umi_val_1.fq.gz 100000 | gzip > ${outputDir}/4_kraken/${Header}_R1_downsampled.fq.gz
    seqtk sample -s100 ${outputDir}/2_UMI_trim/R2_umi_trim_val_2.fq.gz 100000 | gzip > ${outputDir}/4_kraken/${Header}_R2_downsampled.fq.gz

    source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
    conda activate kraken2_env

    kraken2 --db /home/yz2474/yiyong_2023/soft/kraken2/standard_krakendb --threads 20 \
      --report ${outputDir}/4_kraken/${Header}_report.txt --report-minimizer-data \
      --paired --gzip-compressed \
      ${outputDir}/4_kraken/${Header}_R1_downsampled.fq.gz \
      ${outputDir}/4_kraken/${Header}_R2_downsampled.fq.gz

    python /home/yz2474/yiyong_2023/soft/KrakenTools/kreport2krona.py \
      -r ${outputDir}/4_kraken/${Header}_report.txt \
      -o ${outputDir}/4_kraken/${Header}_krona_input.txt

    ktImportText ${outputDir}/4_kraken/${Header}_krona_input.txt -o ${outputDir}/4_kraken/${Header}_krona.html

    conda deactivate
"

# Step 5: STAR alignment
run_step "${outputDir}/5_STAR" "
    ulimit -n 65535
    STAR --genomeDir $genomeDir \
         --readFilesIn ${outputDir}/2_UMI_trim/R1_umi_val_1.fq.gz ${outputDir}/2_UMI_trim/R2_umi_trim_val_2.fq.gz \
         --outFileNamePrefix ${outputDir}/5_STAR/${Header}_ \
         --runThreadN 20 --twopassMode Basic --runMode alignReads --quantMode GeneCounts \
         --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
         --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \
         --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN 20
    samtools index -@ 10 ${outputDir}/5_STAR/${Header}_Aligned.sortedByCoord.out.bam
    python /home/yz2474/scripts/donglab/run_infer_experiment_ZYY.py ${outputDir}/5_STAR/${Header}_Aligned.sortedByCoord.out.bam
"

# Step 6: UMI deduplication
run_step "${outputDir}/6_umi_dedup_bam" "
    umi_tools dedup -I ${outputDir}/5_STAR/${Header}_Aligned.sortedByCoord.out.bam \
      -S ${outputDir}/6_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam \
      --log=${outputDir}/6_umi_dedup_bam/${Header}_umi_dedup.log \
      --extract-umi-method=read_id --paired

    samtools index -@ 20 ${outputDir}/6_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam
    samtools stats ${outputDir}/6_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam \
      > ${outputDir}/6_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup_samtools_stats.txt
    multiqc ${outputDir}/6_umi_dedup_bam -o ${outputDir}/6_umi_dedup_bam
"

# Step 8: FeatureCounts
run_step "${outputDir}/8_featurecounts" "
    featureCounts -a $Gene_GTF \
                  -o ${outputDir}/8_featurecounts/${Header}_featureCounts.txt \
                  -p -T 20 -s 2 -g gene_id -t exon -B -C -M --fraction --countReadPairs \
                  ${outputDir}/6_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam
"

# Step 9: CIRCexplorer2
run_step "${outputDir}/9_CIRCexplorer2" "
    CIRCexplorer2 parse -t STAR -b ${outputDir}/9_CIRCexplorer2/${Header}_back_spliced_junction.bed \
      ${outputDir}/5_STAR/${Header}_Chimeric.out.junction

    CIRCexplorer2 annotate -r $refflat -g $genome_fa \
      -b ${outputDir}/9_CIRCexplorer2/${Header}_back_spliced_junction.bed \
      -o ${outputDir}/9_CIRCexplorer2/${Header}_circularRNA_known.txt
"

# Step 7: bigWig density
# We'll define python_code *outside* the run_step command to avoid nested quotes issues.

python_code=$(cat <<'EOF_PY'
import matplotlib.pyplot as plt
import sys

num_colors = int(sys.argv[1])
palette = plt.colormaps['tab10']
colors = []
for color in palette.colors[:num_colors]:
    r, g, b = color[:3]
    colors.append("#%02x%02x%02x" % (int(r*255), int(g*255), int(b*255)))
print(" ".join(colors))
EOF_PY
)

run_step "${outputDir}/7_bw_density" "
    bamCoverage -b ${outputDir}/6_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam \
      -o ${outputDir}/7_bw_density/${Header}_umi_dedup_unstrand.bw \
      --normalizeUsing BPM --numberOfProcessors 20

    bw_file=${outputDir}/7_bw_density/${Header}_umi_dedup_unstrand.bw
    output_dir=\"\$(dirname \"\$bw_file\")\"

    bed_files=(
\"/mnt/data/projects/donglab/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/HG38_promoter_1500_500bp_noOverlap.bed\"
\"/mnt/data/projects/donglab/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/HG38_5UTR_noOverlap.bed\"
\"/mnt/data/projects/donglab/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/HG38_exon_noOverlap.bed\"
\"/mnt/data/projects/donglab/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/HG38_intron_noOverlap.bed\"
\"/mnt/data/projects/donglab/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/HG38_3UTR_noOverlap.bed\"
\"/mnt/data/projects/donglab/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/HG38_downstream_2kb_noOverlap.bed\"
\"/mnt/data/projects/donglab/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/HG38_intergenic_noOverlap.bed\"
    )

    output_prefix=\$(basename \"\$bw_file\" .bw)
    num_colors=\${#bed_files[@]}

    # We'll just run python_code from the variable
    color_list=\$(python -c \"\$python_code\" \"\$num_colors\")
    IFS=' ' read -r -a colors <<< \"\$color_list\"

    echo \"Running computeMatrix...\"
    computeMatrix scale-regions \
        -S \"\$bw_file\" \
        -R \"\${bed_files[@]}\" \
        --beforeRegionStartLength 0 \
        --regionBodyLength 1000 \
        --afterRegionStartLength 0 \
        -o \"\$output_dir/\${output_prefix}_bed_stacked_matrix.gz\" \
        -p 20 \
        --outFileSortedRegions \"\$output_dir/\${output_prefix}_bed_stacked_sorted_regions.bed\"

    # Updated labels so they don't have single quotes:
    labels=(\"Promoters\" \"5_UTR\" \"Exon\" \"Intron\" \"3_UTR\" \"Downstream\" \"Intergenic\")
    regions_label=()
    for i in \"\${!bed_files[@]}\"; do
        basefile=\$(basename \"\${bed_files[i]}\")
        count=\$(grep -c \"\$basefile\" \"\$output_dir/\${output_prefix}_bed_stacked_sorted_regions.bed\")
        regions_label+=(\"\${labels[i]} (\$count)\")
    done

    echo \"Running plotProfile...\"
    plotProfile -m \"\$output_dir/\${output_prefix}_bed_stacked_matrix.gz\" \
        -out \"\$output_dir/\${output_prefix}_bed_stacked_density_profile.pdf\" \
        --colors \"\${colors[@]}\" \
        --legendLocation upper-left \
        --startLabel \"Start (5')\" \
        --endLabel \"End (3')\" \
        --plotWidth 12 --plotHeight 15 --dpi 300 \
        --yAxisLabel \"Signal intensity (BPM) over scaled regions (200 bp)\" \
        --plotType lines \
        --numPlotsPerRow 1 \
        --regionsLabel \"\${regions_label[@]}\"

    echo \"Profile: \${output_prefix}_bed_stacked_density_profile.pdf\"
"


# Step 18: Generate Output Report using R Markdown
run_step "${outputDir}/18_Output_report" "
    mkdir -p ${outputDir}/18_Output_report

    # Run the Rmd script to generate the HTML report
    Rscript -e 'rmarkdown::render(\"EVscope_report.Rmd\",
    output_file=\"${outputDir}/18_Output_report/EVscope_report.html\")'

    echo \"HTML output report for sample '${Header}' generated successfully in '${outputDir}/18_Output_report/Final_Report.html'\"
"

echo "Pipeline completed successfully for sample: $Header"

