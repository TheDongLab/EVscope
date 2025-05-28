#!/usr/bin/env bash
# Usage: bash EV-RNA-Profiler_v1.sh <SampleName> <R1.fq.gz> <R2.fq.gz>
set -eo pipefail
export TMPDIR=/home/yz2474/tmp
SampleName="$1"
R1_path="$2"
R2_path="$3"
STAR_genomeDir="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_STAR-2.7.10b_sjdbOverhang100_Gene"
combined_gene_gtf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/HG38_3631522_combined_RNAs_with_gold_standard_piRNAs.gtf"
combined_gene_refflat="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/HG38_3631522_combined_RNAs_with_gold_standard_piRNAs.refflat"
GeneID_meta_table="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/3631522_HG38_geneID_Symbol_RNAtype_with_gold_standard_piRNAs.tsv"
genecodeV45_gene_gtf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.gtf"
genecodeV45_gene_refflat="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.refflat"
genome_fa="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Sequence/hg38p14/hg38.p14.whole.genome.fa"
bwa_index_ref="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_BWA_V0.7.18/hg38.p14.whole.genome.fa"
rRNA_interval="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/rRNA/GeneCodeV45_combined_Nucl_Mt_rRNA_589.interval_list"
ref_myco="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/Bacterial/240_Mycoplasma_ge/Mycoplasma_240_ge.fa"
ref_ecoli="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/Bacterial/E.coli/E.coli_GCF_000005845.2_ASM584v2_ge.fa"
HG38_3UTR_saf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_3UTR_noOverlap.saf"
HG38_5UTR_saf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_5UTR_noOverlap.saf"
HG38_downstream_2kb_saf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_downstream_2kb_noOverlap.saf"
HG38_exon_saf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_exon_noOverlap.saf"
HG38_intergenic_saf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intergenic_noOverlap.saf"
HG38_intron_saf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intron_noOverlap.saf"
HG38_promoter_1500_500bp_saf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_promoter_1500_500bp_noOverlap.saf"
HG38_blacklist_saf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_ENCODE_blacklist_V2.saf"
HG38_3UTR_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_3UTR_noOverlap.bed"
HG38_5UTR_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_5UTR_noOverlap.bed"
HG38_downstream_2kb_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_downstream_2kb_noOverlap.bed"
HG38_exon_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_exon_noOverlap.bed"
HG38_intergenic_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intergenic_noOverlap.bed"
HG38_intron_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intron_noOverlap.bed"
HG38_promoter_1500_500bp_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_promoter_1500_500bp_noOverlap.bed"
HG38_blacklist_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_ENCODE_blacklist_V2.bed"

run_step() {
    local step_dir="$1"
    local step_cmds="$2"
    mkdir -p "$step_dir"
    if [ ! -f "$step_dir/step.done" ]; then
        echo "======> Starting step in $step_dir"
        start_time=$(date +%s)
        if eval "$step_cmds"; then
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            touch "$step_dir/step.done"
            echo "======> Finished step: $step_dir in $elapsed seconds"
        else
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            echo "Error: Step failed in $step_dir after $elapsed seconds. Check logs for details."
            exit 1
        fi
    else
        echo "======> Already finished step: $step_dir"
    fi
}

outputDir="${SampleName}"
mkdir -p "${outputDir}"

########################################
# Step 1: Raw FastQC QC
########################################
run_step "${outputDir}/01_raw_fastqc" "
fastqc -o ${outputDir}/01_raw_fastqc -t 20 $R1_path $R2_path
"
########################################
# Step 2: Extract 14-bp UMI motif from Read2
########################################
run_step "${outputDir}/02_UMI_motif" "
python /home/yz2474/scripts/donglab/get_SMARTer_Read2_14BP_motif/get_UMI_seq_motif_for_individual_fq.py \
  -head ${SampleName} \
  -fq $R2_path \
  -n 14 -r 100000 \
  -o ${outputDir}/02_UMI_motif \
  -p 20
"

########################################
# Step 3: UMI labeling for Read1/2 and adapter trimming
########################################
run_step "${outputDir}/03_trim_galore" "
umi_tools extract --bc-pattern='NNNNNNNN' \
  --stdin=$R2_path \
  --stdout=${outputDir}/03_trim_galore/${SampleName}_R2_UMI.fq.gz \
  --read2-in=$R1_path \
  --read2-out=${outputDir}/03_trim_galore/${SampleName}_R1_UMI.fq.gz \
  --log=${outputDir}/03_trim_galore/UMI_extract.log \
  --umi-separator='_'

cutadapt -u 14 -o ${outputDir}/03_trim_galore/${SampleName}_R2_UMI_trim.fq.gz \
  ${outputDir}/03_trim_galore/${SampleName}_R2_UMI.fq.gz
trim_galore --phred33 --illumina --gzip --stringency 3 --cores 20 --paired --length 10 \
  ${outputDir}/03_trim_galore/${SampleName}_R1_UMI.fq.gz \
  ${outputDir}/03_trim_galore/${SampleName}_R2_UMI_trim.fq.gz \
  -o ${outputDir}/03_trim_galore

mv ${outputDir}/03_trim_galore/${SampleName}_R1_UMI_val_1.fq.gz ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz
mv ${outputDir}/03_trim_galore/${SampleName}_R2_UMI_trim_val_2.fq.gz ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz

rm ${outputDir}/03_trim_galore/${SampleName}_R1_UMI.fq.gz
rm ${outputDir}/03_trim_galore/${SampleName}_R2_UMI.fq.gz
rm ${outputDir}/03_trim_galore/${SampleName}_R2_UMI_trim.fq.gz
"

########################################
# Step 4: QC after trim_galore
########################################
run_step "${outputDir}/04_trim_galore_fastqc" "
fastqc -o ${outputDir}/04_trim_galore_fastqc -t 20 \
  ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz \
  ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz
"

########################################
# Step 5: BBSplit (Bacterial read detection)
########################################
run_step "${outputDir}/05_BBSplit" "
bash /home/yz2474/yiyong_2023/soft/bbmap/bbsplit.sh \
  build=1 \
  in1=${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz \
  in2=${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz \
  ref=$ref_ecoli,$ref_myco \
  basename=${outputDir}/05_BBSplit/${SampleName}_%_R#.fq.gz \
  ambiguous=best \
  path=${outputDir}/05_BBSplit
rm -r ${outputDir}/05_BBSplit/ref
"

########################################
# Step 6: first_alignment by using STAR + UMI deduplication
########################################
run_step "${outputDir}/06_first_STAR_UMI_tools" "
#STAR alignment
ulimit -n 65535
STAR --genomeDir $STAR_genomeDir \
  --readFilesIn ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz \
  --outFileNamePrefix ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_ \
  --runThreadN 20 --twopassMode Basic --runMode alignReads --quantMode GeneCounts \
  --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
  --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \
  --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN 20
samtools index -@ 20 ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord.out.bam
#UMI deduplication
umi_tools dedup \
  -I ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord.out.bam \
  -S ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord_umi_dedup.out.bam \
  --log=${outputDir}/06_first_STAR_UMI_tools/${SampleName}_umi_dedup.log \
  --extract-umi-method=read_id --paired
samtools index -@ 20 ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord_umi_dedup.out.bam

#Extract UMI-dedup reads for downstream analysis
samtools view -@ 20 ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord_umi_dedup.out.bam | \
  cut -f1 | sort | uniq > ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_readnames.txt
seqtk subseq ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz \
  ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_readnames.txt | gzip > ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_R1_umi_dedup.clean.fq.gz
seqtk subseq ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz \
  ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_readnames.txt | gzip > ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_R2_umi_dedup.clean.fq.gz

rm ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord.out.bam*
"

########################################
# Step 7: 2nd STAR alignment
########################################
run_step "${outputDir}/07_second_STAR" "
#STAR alignment
ulimit -n 65535
STAR --genomeDir $STAR_genomeDir \
  --readFilesIn ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_R1_umi_dedup.clean.fq.gz ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_R2_umi_dedup.clean.fq.gz \
  --outFileNamePrefix ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_ \
  --runThreadN 20 --twopassMode Basic --runMode alignReads --quantMode GeneCounts \
  --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
  --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \
  --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN 20
samtools index -@ 20 ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam
samtools flagstat -@ 20 ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam > ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.flagstat
python /home/yz2474/scripts/donglab/run_infer_experiment_ZYY.py ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam
"


########################################
# Step 8: CIRCexplorer2
########################################
run_step "${outputDir}/08_CIRCexplorer2" "
CIRCexplorer2 parse -t STAR \
  -b ${outputDir}/08_CIRCexplorer2/${SampleName}_STAR_umi_dedup_back_spliced_junction.bed \
  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Chimeric.out.junction

CIRCexplorer2 annotate \
  -r $genecodeV45_gene_refflat \
  -g $genome_fa \
  -b ${outputDir}/08_CIRCexplorer2/${SampleName}_STAR_umi_dedup_back_spliced_junction.bed \
  -o ${outputDir}/08_CIRCexplorer2/${SampleName}_STAR_umi_dedup_circularRNA_known.txt

python /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/convert_CIRCexplorer2_out2CPM.py  --CIRCexplorer2_result ${outputDir}/08_CIRCexplorer2/${SampleName}_STAR_umi_dedup_circularRNA_known.txt --input_bam ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam --output ${outputDir}/08_CIRCexplorer2/${SampleName}_CIRCexplorer2_dedup_junction_readcounts_CPM.tsv

"


########################################
# Step 9: Alignment Using BWA + CIRI2
########################################
run_step "${outputDir}/09_BWA_CIRI2" "
bwa mem -t 20 -T 19 $bwa_index_ref ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_R1_umi_dedup.clean.fq.gz ${outputDir}/06_first_STAR_UMI_tools/${SampleName}_R2_umi_dedup.clean.fq.gz > ${outputDir}/09_BWA_CIRI2/${SampleName}_umi_dedup.bwa.sam
curr_dir=\$(pwd)
cd ${outputDir}/09_BWA_CIRI2
perl /home/yz2474/yiyong_2023/soft/CIRI_v2.0.6/CIRI2.pl -T 19 -I ${SampleName}_umi_dedup.bwa.sam -O ${SampleName}_CIRI2_out.tsv -F $genome_fa -A $genecodeV45_gene_gtf -G ${SampleName}_CIRI2.log
cd \$curr_dir
python /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/convert_CIRI2_out2CPM.py  --CIRI2_result ${outputDir}/09_BWA_CIRI2/${SampleName}_CIRI2_out.tsv --input_sam ${outputDir}/09_BWA_CIRI2/${SampleName}_umi_dedup.bwa.sam --output ${outputDir}/09_BWA_CIRI2/${SampleName}_CIRI2_dedup_junction_readcounts_CPM.tsv
"

########################################
# Step 10: merge_CIRI2_CIRCexplorer2
########################################
run_step "${outputDir}/10_merge_CIRI2_CIRCexplorer2" "
python /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/Combined_CIRCexplorer2_CIRI2.py --CIRCexplorer2 ${outputDir}/08_CIRCexplorer2/${SampleName}_STAR_umi_dedup_circularRNA_known.txt --CIRI2 ${outputDir}/09_BWA_CIRI2/${SampleName}_CIRI2_out.tsv --output ${outputDir}/10_merge_CIRI2_CIRCexplorer2/${SampleName}_combined_circRNA_known.txt
"

########################################
# Step 11: Picard RNA Metrics
########################################
run_step "${outputDir}/11_picard_metrics" "
source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
conda activate picard_env
picard -Xmx250g CollectRnaSeqMetrics \
  I=${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam \
  O=${outputDir}/11_picard_metrics/${SampleName}_STAR_UMI_dedup.picard_metrics.tsv \
  REF_FLAT=$genecodeV45_gene_refflat \
  STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
  RIBOSOMAL_INTERVALS=$rRNA_interval
picard CollectInsertSizeMetrics \
  I=${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam \
  O=${outputDir}/11_picard_metrics/${SampleName}_STAR_UMI_dedup.insert_size_metrics.tsv \
  H=${outputDir}/11_picard_metrics/${SampleName}_STAR_UMI_dedup.insert_size_histogram.pdf
conda deactivate
"

########################################
# Step 12: featureCounts for read counting
########################################
run_step "${outputDir}/12_featureCounts" "
#Fragment counts
featureCounts \
  -a $combined_gene_gtf \
  -o ${outputDir}/12_featureCounts/${SampleName}_fragment_level_featureCounts.tsv \
  -p -T 20 -s 2 -g gene_id -t exon -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

#Read counts (each read counted separately for meta-features)
featureCounts \
  -F SAF \
  -a $HG38_3UTR_saf \
  -o ${outputDir}/12_featureCounts/${SampleName}_HG38_3UTR_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_5UTR_saf \
  -o ${outputDir}/12_featureCounts/${SampleName}_HG38_5UTR_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_downstream_2kb_saf \
  -o ${outputDir}/12_featureCounts/${SampleName}_HG38_downstream_2kb_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction  --countReadPairs \
  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_blacklist_saf \
  -o ${outputDir}/12_featureCounts/${SampleName}_HG38_ENCODE_blacklist_V2_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_exon_saf \
  -o ${outputDir}/12_featureCounts/${SampleName}_HG38_exon_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_intergenic_saf \
  -o ${outputDir}/12_featureCounts/${SampleName}_HG38_intergenic_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_intron_saf \
  -o ${outputDir}/12_featureCounts/${SampleName}_HG38_intron_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_promoter_1500_500bp_saf \
  -o ${outputDir}/12_featureCounts/${SampleName}_HG38_promoter_1500_500bp_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam
"

########################################
# Step 13: Expression matrix generation with read counts and normalized TPM/CPM and plot RNA distribution
########################################
run_step "${outputDir}/13_plot_RNA_distribution" "
python /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/FeatureCountsToTPMMatrix.py  --featureCounts_out ${outputDir}/12_featureCounts/${SampleName}_fragment_level_featureCounts.tsv --GeneID_meta_table $GeneID_meta_table --output  ${outputDir}/13_plot_RNA_distribution/${SampleName}_Gene_readcounts_normalized_expression_matrix.tsv

python /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/combined_RNA_circRNA.py \
--gene_expr ${outputDir}/13_plot_RNA_distribution/${SampleName}_Gene_readcounts_normalized_expression_matrix.tsv \
--circRNA_expr ${outputDir}/09_BWA_CIRI2/${SampleName}_CIRI2_dedup_junction_readcounts_CPM.tsv \
--out_matrix ${outputDir}/13_plot_RNA_distribution/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv

#figure 1
python /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/StackedRNATypeCompositionPlot_1subplot.py  --sample_name ${SampleName} --Expr_matrix ${outputDir}/13_plot_RNA_distribution/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv --out_plot ${outputDir}/13_plot_RNA_distribution/${SampleName}_RNA_type_composition_1subplot.pdf

#figure 2
python /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/StackedRNATypeCompositionPlot_2subplots.py --sample_name ${SampleName} --Expr_matrix ${outputDir}/13_plot_RNA_distribution/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv --out_plot ${outputDir}/13_plot_RNA_distribution/${SampleName}_RNA_type_composition_2subplots.pdf

#figure 3
python /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/StackedRNATypeCompositionPlot_20subplots.py --sample_name ${SampleName} --Expr_matrix ${outputDir}/13_plot_RNA_distribution/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv --out_plot ${outputDir}/13_plot_RNA_distribution/${SampleName}_RNA_type_composition_20subplots.pdf
"

########################################
# Step 14: bedtools for read counting in GeneCodeV45 meta-features
########################################
run_step "${outputDir}/14_bedtools_meta_gene" "
samtools view -@ 20 -b -F 0x100  ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam > ${outputDir}/14_bedtools_meta_gene/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.primary.bam
samtools index -@ 20 ${outputDir}/14_bedtools_meta_gene/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.primary.bam
bedtools coverage -a $HG38_3UTR_bed -b ${outputDir}/14_bedtools_meta_gene/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.primary.bam > ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_3UTR_bedtools_cov.tsv
bedtools coverage -a $HG38_5UTR_bed -b ${outputDir}/14_bedtools_meta_gene/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.primary.bam > ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_5UTR_bedtools_cov.tsv
bedtools coverage -a $HG38_downstream_2kb_bed -b ${outputDir}/14_bedtools_meta_gene/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.primary.bam  > ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_downstream_2kb_bedtools_cov.tsv
bedtools coverage -a $HG38_exon_bed -b ${outputDir}/14_bedtools_meta_gene/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.primary.bam > ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_exon_bedtools_cov.tsv
bedtools coverage -a $HG38_intergenic_bed -b ${outputDir}/14_bedtools_meta_gene/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.primary.bam > ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_intergenic_bedtools_cov.tsv
bedtools coverage -a $HG38_intron_bed -b ${outputDir}/14_bedtools_meta_gene/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.primary.bam > ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_intron_bedtools_cov.tsv
bedtools coverage -a $HG38_promoter_1500_500bp_bed -b ${outputDir}/14_bedtools_meta_gene/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.primary.bam > ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_promoter_1500_500bp_bedtools_cov.tsv
bedtools coverage -a $HG38_blacklist_bed -b ${outputDir}/14_bedtools_meta_gene/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.primary.bam > ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_GENECODE_blacklist_bedtools_cov.tsv
"

########################################
# Step 15: Kraken for Taxonomic Classification (Downsampling)
########################################
run_step "${outputDir}/14_Kraken" "
seqtk sample -s100 ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz 100000 | gzip > ${outputDir}/15_Kraken/${SampleName}_R1_downsampled_100K.fq.gz
seqtk sample -s100 ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz 100000 | gzip > ${outputDir}/14_Kraken/${SampleName}_R2_downsampled_100K.fq.gz
source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
conda activate kraken2_env
kraken2 --db /home/yz2474/yiyong_2023/soft/kraken2/standard_krakendb --threads 20 \
  --report ${outputDir}/14_Kraken/${SampleName}_report.tsv --report-minimizer-data \
  --paired --gzip-compressed \
  ${outputDir}/14_Kraken/${SampleName}_R1_downsampled_100K.fq.gz \
  ${outputDir}/14_Kraken/${SampleName}_R2_downsampled_100K.fq.gz
python /home/yz2474/yiyong_2023/soft/KrakenTools/kreport2krona.py \
  -r ${outputDir}/14_Kraken/${SampleName}_report.tsv \
  -o ${outputDir}/14_Kraken/${SampleName}_krona_input.tsv
ktImportText ${outputDir}/14_Kraken/${SampleName}_krona_input.tsv \
  -o ${outputDir}/14_Kraken/${SampleName}_krona.html
conda deactivate
"

########################################
# Step 16: rRNA Detection
########################################
run_step "${outputDir}/16_RiboDetector" "
ribodetector_cpu -t 20 \
  -l 100 \
  -i ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz \
     ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz \
  -e rrna \
  --chunk_size 256 \
  -r ${outputDir}/16_RiboDetector/${SampleName}_rRNA_R1_trimed.fq.gz \
     ${outputDir}/16_RiboDetector/${SampleName}_rRNA_R2_trimed.fq.gz \
  -o /dev/null /dev/null
"

########################################
# Step 17: QC_matrix table generation
########################################
run_step "${outputDir}/17_QC_matrix" "
python /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/QC_matrix.py \
    --raw_R1_fastqc_zip ${outputDir}/01_raw_fastqc/${SampleName}*R1*_fastqc.zip \
    --raw_R2_fastqc_zip ${outputDir}/01_raw_fastqc/${SampleName}*R2*_fastqc.zip \
    --trimmed_R1_fastqc_zip ${outputDir}/04_trim_galore_fastqc/${SampleName}_R1_trimed_fastqc.zip \
    --trimmed_R2_fastqc_zip ${outputDir}/04_trim_galore_fastqc/${SampleName}_R2_trimed_fastqc.zip \
    --ecoli_R1_fastq ${outputDir}/05_BBSplit/${SampleName}_E.coli_GCF_000005845.2_ASM584v2_ge_R1.fq.gz \
    --ecoli_R2_fastq ${outputDir}/05_BBSplit/${SampleName}_E.coli_GCF_000005845.2_ASM584v2_ge_R2.fq.gz \
    --myco_R1_fastq ${outputDir}/05_BBSplit/${SampleName}_Mycoplasma_240_ge_R1.fq.gz \
    --myco_R2_fastq ${outputDir}/05_BBSplit/${SampleName}_Mycoplasma_240_ge_R2.fq.gz \
    --umi_dedup_bam ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam \
    --bam2strand_file ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out_bam2strand_file.xls \
    --STAR_log ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Log.final.out \
    --picard_insert_file ${outputDir}/11_picard_metrics/${SampleName}_STAR_UMI_dedup.insert_size_metrics.tsv \
    --picard_rnaseq_file ${outputDir}/11_picard_metrics/${SampleName}_STAR_UMI_dedup.picard_metrics.tsv \
    --featureCounts_3UTR ${outputDir}/12_featureCounts/${SampleName}_HG38_3UTR_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_5UTR ${outputDir}/12_featureCounts/${SampleName}_HG38_5UTR_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_downstream_2kb ${outputDir}/12_featureCounts/${SampleName}_HG38_downstream_2kb_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_exon ${outputDir}/12_featureCounts/${SampleName}_HG38_exon_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_ENCODE_blacklist ${outputDir}/12_featureCounts/${SampleName}_HG38_ENCODE_blacklist_V2_read_level_featureCounts.tsv \
    --featureCounts_intergenic ${outputDir}/12_featureCounts/${SampleName}_HG38_intergenic_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_intron ${outputDir}/12_featureCounts/${SampleName}_HG38_intron_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_promoter_1500_500bp ${outputDir}/12_featureCounts/${SampleName}_HG38_promoter_1500_500bp_noOverlap_read_level_featureCounts.tsv \
    --expression_matrix ${outputDir}/13_plot_RNA_distribution/${SampleName}_merged_matrix.tsv \
    --bedtools_cov_3UTR ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_3UTR_bedtools_cov.tsv \
    --bedtools_cov_5UTR ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_5UTR_bedtools_cov.tsv \
    --bedtools_cov_downstream_2kb ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_downstream_2kb_bedtools_cov.tsv \
    --bedtools_cov_exon ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_exon_bedtools_cov.tsv \
    --bedtools_cov_ENCODE_blacklist ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_GENECODE_blacklist_bedtools_cov.tsv \
    --bedtools_cov_intergenic ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_intergenic_bedtools_cov.tsv \
    --bedtools_cov_intron ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_intron_bedtools_cov.tsv \
    --bedtools_cov_promoter_1500_500bp ${outputDir}/14_bedtools_meta_gene/${SampleName}_HG38_promoter_1500_500bp_bedtools_cov.tsv \
    --kraken_report ${outputDir}/14_Kraken/${SampleName}_report.tsv \
    --ribo_R1_fastq ${outputDir}/16_RiboDetector/${SampleName}_rRNA_R1_trimed.fq.gz \
    --ribo_R2_fastq ${outputDir}/16_RiboDetector/${SampleName}_rRNA_R2_trimed.fq.gz \
    --output ${outputDir}/17_QC_matrix/${SampleName}_QC_matrix.tsv
"

########################################
# Step 18: bigWig Density Generation
########################################
#run_step "${outputDir}/18_bigWig" "
#bamCoverage -b ${outputDir}/07_second_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam \
  -o ${outputDir}/18_bigWig/${SampleName}_UMI_dedup_unstranded.bw \
  --normalizeUsing BPM --numberOfProcessors 20
#bash /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/run_deeptools_plotprfile_bed_stackd_meta_gene.sh \
  ${outputDir}/18_bigWig/${SampleName}_UMI_dedup_unstranded.bw
#bash /home/yz2474/scripts/donglab/EV-RNA-Profiler/scripts/run_deeptools_plotprfile_bed_stackd_RNA_types.sh \
  ${outputDir}/18_bigWig/${SampleName}_UMI_dedup_unstranded.bw
"
#echo "Pipeline completed successfully for sample: ${SampleName}"
