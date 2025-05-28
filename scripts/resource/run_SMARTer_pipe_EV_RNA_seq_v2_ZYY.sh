#!/usr/bin/env bash
#Usage: bash run_SMARTer_pipe_v3_ZYY.sh B4_TDP43_KD_caRNA_rep3 /home/yz2474/yiyong_2023/DB/example_data/fastq/1_raw/B4_TDP43_KD_caRNA_rep3/B4_TDP43_KD_caRNA_rep3_R1.fq.gz /home/yz2474/yiyong_2023/DB/example_data/fastq/1_raw/B4_TDP43_KD_caRNA_rep3/B4_TDP43_KD_caRNA_rep3_R2.fq.gz

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

# 运行某一步骤，如果step.done不存在则执行并生成step.done
run_step() {
    local step_dir=$1
    echo "runing in " ${step_dir}
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

# Step 1: 原始数据QC
run_step "${outputDir}/1_raw_QC" "
    fastqc -o ${outputDir}/1_raw_QC -t 20 $R1_path $R2_path
    multiqc ${outputDir}/1_raw_QC -o ${outputDir}/1_raw_QC
"

# Step 2: UMI提取与去接头
run_step "${outputDir}/2_UMI_trim" "
    umi_tools extract --bc-pattern='NNNNNNNN' \
      --stdin=$R2_path --stdout=${outputDir}/2_UMI_trim/R2_umi.fq.gz \
      --read2-in=$R1_path --read2-out=${outputDir}/2_UMI_trim/R1_umi.fq.gz \
      --log=${outputDir}/2_UMI_trim/umi_extract.log --umi-separator='_'
    cutadapt -u 14 -o ${outputDir}/2_UMI_trim/R2_umi_trim.fq.gz ${outputDir}/2_UMI_trim/R2_umi.fq.gz
    trim_galore --phred33 --fastqc --illumina --gzip --stringency 3 --cores 20 --paired --length 20 \
      ${outputDir}/2_UMI_trim/R1_umi.fq.gz ${outputDir}/2_UMI_trim/R2_umi_trim.fq.gz \
      -o ${outputDir}/2_UMI_trim
    fastqc -o ${outputDir}/2_UMI_trim -t 20 \
      ${outputDir}/2_UMI_trim/R1_umi_val_1.fq.gz ${outputDir}/2_UMI_trim/R2_umi_trim_val_2.fq.gz
    multiqc ${outputDir}/2_UMI_trim -o ${outputDir}/2_UMI_trim
"

# Step 3: Kraken下采样与物种注释
run_step "${outputDir}/3_kraken" "
    seqtk sample -s100 $R1_path 100000 | gzip > ${outputDir}/3_kraken/${Header}_R1_downsampled.fq.gz
    seqtk sample -s100 $R2_path 100000 | gzip > ${outputDir}/3_kraken/${Header}_R2_downsampled.fq.gz
    source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
    conda activate kraken2_env
    kraken2 --db /home/yz2474/yiyong_2023/soft/kraken2/standard_krakendb --threads 20 --report ${outputDir}/3_kraken/${Header}_report.txt --report-minimizer-data --paired --gzip-compressed ${outputDir}/3_kraken/${Header}_R1_downsampled.fq.gz ${outputDir}/3_kraken/${Header}_R2_downsampled.fq.gz
    python /home/yz2474/yiyong_2023/soft/KrakenTools/kreport2krona.py -r ${outputDir}/3_kraken/${Header}_report.txt -o ${outputDir}/3_kraken/${Header}_krona_input.txt
    ktImportText ${outputDir}/3_kraken/${Header}_krona_input.txt -o ${outputDir}/3_kraken/${Header}_krona.html
    conda deactivate
"

# Step 4: STAR比对
run_step "${outputDir}/4_STAR" "
    ulimit -n 65535
    STAR --genomeDir $genomeDir \
         --readFilesIn ${outputDir}/2_UMI_trim/R1_umi_val_1.fq.gz ${outputDir}/2_UMI_trim/R2_umi_trim_val_2.fq.gz \
         --outFileNamePrefix ${outputDir}/4_STAR/${Header}_ \
         --runThreadN 20 --twopassMode Basic --runMode alignReads --quantMode GeneCounts \
         --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
         --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \
         --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN 20
    samtools index -@ 10 ${outputDir}/4_STAR/${Header}_Aligned.sortedByCoord.out.bam
    python /home/yz2474/scripts/donglab/run_infer_experiment_ZYY.py ${outputDir}/4_STAR/${Header}_Aligned.sortedByCoord.out.bam
"

# Step 5: UMI去重
run_step "${outputDir}/5_umi_dedup_bam" "
    umi_tools dedup -I ${outputDir}/4_STAR/${Header}_Aligned.sortedByCoord.out.bam \
      -S ${outputDir}/5_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam \
      --log=${outputDir}/5_umi_dedup_bam/${Header}_umi_dedup.log \
      --extract-umi-method=read_id --paired
    samtools index -@ 20 ${outputDir}/5_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam
    samtools stats ${outputDir}/5_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam \
      > ${outputDir}/5_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup_samtools_stats.txt
    multiqc ${outputDir}/5_umi_dedup_bam -o ${outputDir}/5_umi_dedup_bam
"

# Step 6: 生成bigWig
run_step "${outputDir}/6_bw_density" "
    bamCoverage -b ${outputDir}/5_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam \
      -o ${outputDir}/6_bw_density/${Header}_umi_dedup_unstrand.bw \
      --normalizeUsing BPM --numberOfProcessors 20
    bash ~/scripts/donglab/run_deeptools_plotprfile_bed_stackd_202412_ZYY.sh ${outputDir}/6_bw_density/${Header}_umi_dedup_unstrand.bw
"

# Step 7: FeatureCounts计数
run_step "${outputDir}/7_featurecounts" "
    featureCounts -a $Gene_GTF \
                  -o ${outputDir}/7_featurecounts/${Header}_featureCounts.txt \
                  -p -T 20 -s 2 -g gene_id -t exon -B -C --countReadPairs \
                  ${outputDir}/5_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam
"

# Step 8: CIRCexplorer2分析
run_step "${outputDir}/8_CIRCexplorer2" "
    CIRCexplorer2 parse -t STAR -b ${outputDir}/8_CIRCexplorer2/${Header}_back_spliced_junction.bed \
      ${outputDir}/4_STAR/${Header}_Chimeric.out.junction
    CIRCexplorer2 annotate -r $refflat -g $genome_fa \
      -b ${outputDir}/8_CIRCexplorer2/${Header}_back_spliced_junction.bed \
      -o ${outputDir}/8_CIRCexplorer2/${Header}_circularRNA_known.txt
"
echo "Pipeline completed successfully for sample: $Header"
