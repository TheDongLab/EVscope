#!/usr/bin/env bash
# Usage: bash EVscope_v1.sh <SampleName> <R1.fq.gz> <R2.fq.gz>
rm -r example_SMARTerv3_output/17_HTML_report
set -eo pipefail
export TMPDIR=/home/yz2474/yiyong_2023/TEMP
curr_dir=$(pwd)
SampleName="$1"
R1_path="$2"
R2_path="$3"
RSEM_bowtie2_index="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/Index/RSEM_index/RSEM_bowtie2_index/RSEM_REF_HG38_3659642_RNNAs"
STAR_genomeDir="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_STAR-2.7.10b_sjdbOverhang99_GeneCode45"
combined_gene_gtf="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/combine/HG38_3659642_combined_RNAs_with_gold_standard_piRNAs.gtf"
combined_gene_refflat="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/combine/HG38_3659642_combined_RNAs_with_gold_standard_piRNAs.refflat"
GeneID_meta_table="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/combine/3659642_HG38_geneID_Symbol_RNAtype_with_gold_standard_piRNAs.tsv"
genecodeV45_gene_gtf="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/combine/gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.gtf"
genecodeV45_gene_refflat="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/combine/gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.refflat"
genecodeV45_bed="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/3_gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.geneID.bed6"
genome_fa="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Sequence/hg38p14/hg38.p14.whole.genome.fa"
bwa_index_ref="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_BWA_V0.7.18/hg38.p14.whole.genome.fa"
rRNA_interval="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/combine/rRNA/GeneCodeV45_combined_Nucl_Mt_rRNA_589.interval_list"
ref_myco="/home/yz2474/scripts/donglab/EVscope/genome_anno/Bacterial/240_Mycoplasma_ge/Mycoplasma_240_ge.fa"
ref_ecoli="/home/yz2474/scripts/donglab/EVscope/genome_anno/Bacterial/E.coli/E.coli_GCF_000005845.2_ASM584v2_ge.fa"
HG38_3UTR_saf="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_3UTR_noOverlap.saf"
HG38_5UTR_saf="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_5UTR_noOverlap.saf"
HG38_downstream_2kb_saf="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_downstream_2kb_noOverlap.saf"
HG38_exon_saf="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_exon_noOverlap.saf"
HG38_intergenic_saf="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intergenic_noOverlap.saf"
HG38_intron_saf="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intron_noOverlap.saf"
HG38_promoter_1500_500bp_saf="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_promoter_1500_500bp_noOverlap.saf"
HG38_blacklist_saf="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_ENCODE_blacklist_V2.saf"
HG38_3UTR_bed="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_3UTR_noOverlap.bed"
HG38_5UTR_bed="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_5UTR_noOverlap.bed"
HG38_downstream_2kb_bed="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_downstream_2kb_noOverlap.bed"
HG38_exon_bed="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_exon_noOverlap.bed"
HG38_intergenic_bed="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intergenic_noOverlap.bed"
HG38_intron_bed="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intron_noOverlap.bed"
HG38_promoter_1500_500bp_bed="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_promoter_1500_500bp_noOverlap.bed"
HG38_blacklist_bed="/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_ENCODE_blacklist_V2.bed"

run_step() {
    local step_dir="$1"
    local step_cmds="$2"
    local ignore_error="${3:-false}"  # 新参数，默认为false

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
            echo "Error: Step failed in   $step_dir after $elapsed seconds. Check logs for details."

            if [ "$ignore_error" != "true" ]; then
                exit 1
            else
                echo "Ignoring the error as per configuration, continuing pipeline."
            fi
        fi
    else
        echo "======> Already finished step: $step_dir"
    fi
}


outputDir="${SampleName}_output"
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
python /home/yz2474/scripts/donglab/EVscope/scripts/01_plot_fastq2UMI_motif.py \
  -head ${SampleName} \
  -fq $R2_path \
  -n 14 -r 1000000 \
  -o ${outputDir}/02_UMI_motif
"

########################################
# Step 3: UMI labeling for Read1/2 and adapter trimming
########################################
run_step "${outputDir}/03_cutadapt" "
umi_tools extract --bc-pattern='NNNNNNNNNNNNNN' \
  --stdin=$R2_path \
  --stdout=${outputDir}/03_cutadapt/${SampleName}_R2_umi_tools.fq.gz \
  --read2-in=$R1_path \
  --read2-out=${outputDir}/03_cutadapt/${SampleName}_R1_umi_tools.fq.gz \
  --log=${outputDir}/03_cutadapt/UMI_extract.log \
  --umi-separator='_'

cutadapt \
  -a AGATCGGAAGAGC \
  -A AGATCGGAAGAGC \
  --overlap 3 \
  --minimum-length 10 \
  -j 20 \
  -o ${outputDir}/03_cutadapt/${SampleName}_R1_adapter_trimmed.fq.gz \
  -p ${outputDir}/03_cutadapt/${SampleName}_R2_adapter_trimmed.fq.gz \
  ${outputDir}/03_cutadapt/${SampleName}_R1_umi_tools.fq.gz \
  ${outputDir}/03_cutadapt/${SampleName}_R2_umi_tools.fq.gz

python ~/scripts/donglab/EVscope/scripts/UMIAdapterTrimR1.py \
    --input_R1_fq ${outputDir}/03_cutadapt/${SampleName}_R1_adapter_trimmed.fq.gz \
    --input_R2_fq ${outputDir}/03_cutadapt/${SampleName}_R2_adapter_trimmed.fq.gz \
    --output_R1_fq ${outputDir}/03_cutadapt/${SampleName}_R1_adapter_UMI_trimmed.fq.gz \
    --output_R2_fq ${outputDir}/03_cutadapt/${SampleName}_R2_adapter_UMI_trimmed.fq.gz \
    --output_tsv ${outputDir}/03_cutadapt/${SampleName}_R1_readthrough_UMI_trimming.log \
    --min-overlap 3 \
    --min-length 10 \
    --chunk-size 100000\
    --error-rate 0.1

cutadapt \
  -q 20 \
  --minimum-length 10 \
  -j 20 \
  -o ${outputDir}/03_cutadapt/${SampleName}_R1_clean.fq.gz \
  -p ${outputDir}/03_cutadapt/${SampleName}_R2_clean.fq.gz \
  ${outputDir}/03_cutadapt/${SampleName}_R1_adapter_UMI_trimmed.fq.gz \
  ${outputDir}/03_cutadapt/${SampleName}_R2_adapter_UMI_trimmed.fq.gz
"

########################################
# Step 4: QC after trim_galore
########################################
run_step "${outputDir}/04_cutadapt_fastqc" "
fastqc -o ${outputDir}/04_cutadapt_fastqc -t 20 \
  ${outputDir}/03_cutadapt/${SampleName}_R1_clean.fq.gz \
  ${outputDir}/03_cutadapt/${SampleName}_R2_clean.fq.gz
"

########################################
# Step 5: BBSplit (Bacterial read detection)
########################################
run_step "${outputDir}/05_BBSplit" "
bash /home/yz2474/yiyong_2023/soft/bbmap/bbsplit.sh \
  build=1 \
  threads=20 \
  in1=${outputDir}/03_cutadapt/${SampleName}_R1_clean.fq.gz \
  in2=${outputDir}/03_cutadapt/${SampleName}_R2_clean.fq.gz \
  ref=$ref_ecoli,$ref_myco \
  basename=${outputDir}/05_BBSplit/${SampleName}_%_R#.fq.gz \
  ambiguous=best \
  path=${outputDir}/05_BBSplit
rm -r ${outputDir}/05_BBSplit/ref
"

########################################
# Step 6-1: STAR + UMI deduplication
########################################
run_step "${outputDir}/06_STAR/1st_STAR_UMI_tools" "
#STAR alignment
ulimit -n 65535
STAR --genomeDir $STAR_genomeDir \
  --readFilesIn ${outputDir}/03_cutadapt/${SampleName}_R1_clean.fq.gz ${outputDir}/03_cutadapt/${SampleName}_R2_clean.fq.gz \
  --outFileNamePrefix ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_ \
  --runThreadN 20 --twopassMode Basic --runMode alignReads \
  --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
  --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \
  --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN 20
samtools index -@ 20 ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord.out.bam
#UMI deduplication
umi_tools dedup \
  -I ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord.out.bam \
  -S ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord_umi_dedup.out.bam \
  --log=${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_umi_dedup.log \
  --extract-umi-method=read_id --paired
samtools index -@ 20 ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord_umi_dedup.out.bam

#Extract UMI-dedup reads for downstream analysis
samtools view -@ 20 ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord_umi_dedup.out.bam | \
  cut -f1 | sort | uniq > ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_readnames.txt
seqtk subseq ${outputDir}/03_cutadapt/${SampleName}_R1_clean.fq.gz \
  ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_readnames.txt | gzip > ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_R1_umi_dedup.clean.fq.gz
seqtk subseq ${outputDir}/03_cutadapt/${SampleName}_R2_clean.fq.gz \
  ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_readnames.txt | gzip > ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_R2_umi_dedup.clean.fq.gz

rm ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_Aligned.sortedByCoord.out.bam
"

########################################
# Step 6-2: 2nd STAR alignment
########################################
run_step "${outputDir}/06_STAR/2nd_STAR" "
#STAR alignment
ulimit -n 65535
STAR --genomeDir $STAR_genomeDir \
  --readFilesIn ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_R1_umi_dedup.clean.fq.gz ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_R2_umi_dedup.clean.fq.gz \
  --outFileNamePrefix ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_ \
  --runThreadN 20 --twopassMode Basic --runMode alignReads --quantMode GeneCounts \
  --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
  --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \
  --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN 20
samtools index -@ 20 ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam
samtools flagstat -@ 20 ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam > ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.flagstat
"

########################################
# step 7 RNA strandness detection
########################################
run_step "${outputDir}/07_strand_detection" "
cd ${outputDir}/07_strand_detection
python /home/yz2474/scripts/donglab/EVscope/scripts/02_bam2strand1.py --input_bam ../06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam --bed $genecodeV45_bed --test_read_num 100000000
cd \$curr_dir
"

########################################
# Step 8-1: circRNA_detection by CIRCexplorer2
########################################
run_step "${outputDir}/08_circRNA_detection/CIRCexplorer2" "
CIRCexplorer2 parse -t STAR \
  -b ${outputDir}/08_circRNA_detection/CIRCexplorer2/${SampleName}_STAR_umi_dedup_back_spliced_junction.bed \
  ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Chimeric.out.junction

CIRCexplorer2 annotate \
  -r $genecodeV45_gene_refflat \
  -g $genome_fa \
  -b ${outputDir}/08_circRNA_detection/CIRCexplorer2/${SampleName}_STAR_umi_dedup_back_spliced_junction.bed \
  -o ${outputDir}/08_circRNA_detection/CIRCexplorer2/${SampleName}_STAR_umi_dedup_circularRNA_known.txt

python /home/yz2474/scripts/donglab/EVscope/scripts/03_convert_CIRCexplorer2CPM.py  --CIRCexplorer2_result ${outputDir}/08_circRNA_detection/CIRCexplorer2/${SampleName}_STAR_umi_dedup_circularRNA_known.txt --input_bam ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam --GeneID_meta_table $GeneID_meta_table --output ${outputDir}/08_circRNA_detection/CIRCexplorer2/${SampleName}_CIRCexplorer2_dedup_junction_readcounts_CPM.tsv 
"

########################################
# Step 8-2: circRNA_detection by Using BWA + CIRI2
########################################
run_step "${outputDir}/08_circRNA_detection/BWA_CIRI2" "
bwa mem -t 20 -T 19 $bwa_index_ref ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_R1_umi_dedup.clean.fq.gz ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_R2_umi_dedup.clean.fq.gz > ${outputDir}/08_circRNA_detection/BWA_CIRI2/${SampleName}_umi_dedup.bwa.sam
cd ${outputDir}/08_circRNA_detection/BWA_CIRI2
perl /home/yz2474/yiyong_2023/soft/CIRI_v2.0.6/CIRI2.pl -T 19 -I ${SampleName}_umi_dedup.bwa.sam -O ${SampleName}_CIRI2_out.tsv -F $genome_fa -A $genecodeV45_gene_gtf -G ${SampleName}_CIRI2.log
cd \$curr_dir
python /home/yz2474/scripts/donglab/EVscope/scripts/04_convert_CIRI2CPM.py  --CIRI2_result ${outputDir}/08_circRNA_detection/BWA_CIRI2/${SampleName}_CIRI2_out.tsv --input_sam ${outputDir}/08_circRNA_detection/BWA_CIRI2/${SampleName}_umi_dedup.bwa.sam --output ${outputDir}/08_circRNA_detection/BWA_CIRI2/${SampleName}_CIRI2_dedup_junction_readcounts_CPM.tsv --GeneID_meta_table $GeneID_meta_table
rm ${outputDir}/08_circRNA_detection/BWA_CIRI2/${SampleName}_umi_dedup.bwa.sam
"

########################################
# Step 8-3: merge_CIRI2_CIRCexplorer2
########################################
run_step "${outputDir}/08_circRNA_detection/merge_CIRI2_CIRCexplorer2" "
python /home/yz2474/scripts/donglab/EVscope/scripts/05_combined_CIRCexplorer2_CIRI2.py --CIRCexplorer2 ${outputDir}/08_circRNA_detection/CIRCexplorer2/${SampleName}_CIRCexplorer2_dedup_junction_readcounts_CPM.tsv --CIRI2 ${outputDir}/08_circRNA_detection/BWA_CIRI2/${SampleName}_CIRI2_dedup_junction_readcounts_CPM.tsv --output_matrix ${outputDir}/08_circRNA_detection/merge_CIRI2_CIRCexplorer2/${SampleName}_combined_CIRCexplorer2_CIRI2.tsv --out_venn ${outputDir}/08_circRNA_detection/merge_CIRI2_CIRCexplorer2/${SampleName}_Venn_diagram_of_circRNAs_identified_between_CIRCexplorer2_CIRI2.png
"

########################################
# Step 9: Picard RNA Metrics
########################################
run_step "${outputDir}/09_picard_metrics" "
conda run -n picard_env picard -Xmx250g CollectRnaSeqMetrics \
  I=${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam \
  O=${outputDir}/09_picard_metrics/${SampleName}_STAR_UMI_dedup.picard_metrics.tsv \
  REF_FLAT=$genecodeV45_gene_refflat \
  STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
  RIBOSOMAL_INTERVALS=$rRNA_interval
conda run -n picard_env picard CollectInsertSizeMetrics \
  I=${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam \
  O=${outputDir}/09_picard_metrics/${SampleName}_STAR_UMI_dedup.insert_size_metrics.tsv \
  H=${outputDir}/09_picard_metrics/${SampleName}_STAR_UMI_dedup.insert_size_histogram.pdf

conda run -n picard_env convert -density 300 \
  -background white -alpha remove \
  ${outputDir}/09_picard_metrics/${SampleName}_STAR_UMI_dedup.insert_size_histogram.pdf \
  ${outputDir}/09_picard_metrics/${SampleName}_STAR_UMI_dedup.insert_size_histogram.png

"

########################################
# Step 10-1: featureCounts for read counts
########################################
run_step "${outputDir}/10_RNA_expr_quantification/featureCounts" "
#Fragment counts
featureCounts \
  -a $combined_gene_gtf \
  -o ${outputDir}/10_RNA_expr_quantification/featureCounts/${SampleName}_fragment_level_featureCounts.tsv \
  -p -T 20 -s 2 -g gene_id -t exon -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam
"

########################################
# Step 10-2: RSEM for read counts
########################################
run_step "${outputDir}/10_RNA_expr_quantification/RSEM" " 
/home/yz2474/yiyong_2023/soft/RSEM_v1.3.3/rsem-calculate-expression \
  --paired-end \
  --bowtie2 \
  --strandedness reverse \
  --bowtie2-k 2 \
  -p 20 \
  --no-bam-output \
  --seed 12345 \
  ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_R1_umi_dedup.clean.fq.gz \
  ${outputDir}/06_STAR/1st_STAR_UMI_tools/${SampleName}_R2_umi_dedup.clean.fq.gz \
  $RSEM_bowtie2_index \
  ${outputDir}/10_RNA_expr_quantification/RSEM/${SampleName}_RSEM \
  --temporary-folder ${outputDir}/10_RNA_expr_quantification/RSEM/tmp \
  2> ${outputDir}/10_RNA_expr_quantification/RSEM/${SampleName}_RSEM.log"

########################################
# Step 11-1: FeatureCounts: Expression matrix generation with read counts and normalized TPM/CPM and plot RNA distribution
########################################
run_step "${outputDir}/11_RNA_distribution/featureCounts/" "
python /home/yz2474/scripts/donglab/EVscope/scripts/06_featureCounts2TPM.py  --featureCounts_out ${outputDir}/10_RNA_expr_quantification/featureCounts/${SampleName}_fragment_level_featureCounts.tsv --GeneID_meta_table $GeneID_meta_table --output  ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_Gene_readcounts_normalized_expression_matrix_featureCounts.tsv
#Combine linear RNA and circRNA expression matrix

python /home/yz2474/scripts/donglab/EVscope/scripts/07_combine_total_RNA_expr_matrix.py \
--gene_expr ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_Gene_readcounts_normalized_expression_matrix_featureCounts.tsv \
--circRNA_expr ${outputDir}/08_circRNA_detection/merge_CIRI2_CIRCexplorer2/${SampleName}_combined_CIRCexplorer2_CIRI2.tsv \
--out_matrix ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv

#figure 1
python /home/yz2474/scripts/donglab/EVscope/scripts/08_plot_RNA_distribution_1subplot.py  --sample_name ${SampleName} --Expr_matrix ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv --out_plot ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_RNA_type_composition_1subplot.pdf
#figure 2
python /home/yz2474/scripts/donglab/EVscope/scripts/09_plot_RNA_distribution_2subplots.py --sample_name ${SampleName} --Expr_matrix ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv --out_plot ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_RNA_type_composition_2subplots.pdf
#figure 3
python /home/yz2474/scripts/donglab/EVscope/scripts/10_plot_RNA_distribution_20subplots.py --sample_name ${SampleName} --Expr_matrix ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv --out_plot ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_RNA_type_composition_20subplots.pdf
"

########################################
# Step 11-2: RSEM: Expression matrix generation with read counts and normalized TPM/CPM and plot RNA distribution
########################################
run_step "${outputDir}/11_RNA_distribution/RSEM" "
python ~/scripts/donglab/EVscope/scripts/06_RSEM2expr_matrix.py --RSEM_out ${outputDir}/10_RNA_expr_quantification/RSEM/${SampleName}_RSEM.genes.results --GeneID_meta_table $GeneID_meta_table --output  ${outputDir}/11_RNA_distribution/RSEM/${SampleName}_Gene_readcounts_normalized_expression_matrix_RSEM.tsv

python /home/yz2474/scripts/donglab/EVscope/scripts/07_combine_total_RNA_expr_matrix.py \
--gene_expr ${outputDir}/11_RNA_distribution/RSEM/${SampleName}_Gene_readcounts_normalized_expression_matrix_RSEM.tsv \
--circRNA_expr ${outputDir}/08_circRNA_detection/merge_CIRI2_CIRCexplorer2/${SampleName}_combined_CIRCexplorer2_CIRI2.tsv \
--out_matrix ${outputDir}/11_RNA_distribution/RSEM/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv

#figure 1
python /home/yz2474/scripts/donglab/EVscope/scripts/08_plot_RNA_distribution_1subplot.py  --sample_name ${SampleName} --Expr_matrix ${outputDir}/11_RNA_distribution/RSEM/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv --out_plot ${outputDir}/11_RNA_distribution/RSEM/${SampleName}_RNA_type_composition_1subplot.pdf
#figure 2
python /home/yz2474/scripts/donglab/EVscope/scripts/09_plot_RNA_distribution_2subplots.py --sample_name ${SampleName} --Expr_matrix ${outputDir}/11_RNA_distribution/RSEM/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv --out_plot ${outputDir}/11_RNA_distribution/RSEM/${SampleName}_RNA_type_composition_2subplots.pdf
#figure 3
python /home/yz2474/scripts/donglab/EVscope/scripts/10_plot_RNA_distribution_20subplots.py --sample_name ${SampleName} --Expr_matrix ${outputDir}/11_RNA_distribution/RSEM/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv --out_plot ${outputDir}/11_RNA_distribution/RSEM/${SampleName}_RNA_type_composition_20subplots.pdf
"

########################################
# Step 12: meta-features
########################################
run_step "${outputDir}/12_reads_mapping_stats" "
#Read counts (each read counted separately for meta-features)
featureCounts \
  -F SAF \
  -a $HG38_3UTR_saf \
  -o ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_3UTR_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_5UTR_saf \
  -o ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_5UTR_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_downstream_2kb_saf \
  -o ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_downstream_2kb_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction  --countReadPairs \
  ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_blacklist_saf \
  -o ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_ENCODE_blacklist_V2_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_exon_saf \
  -o ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_exon_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_intergenic_saf \
  -o ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_intergenic_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_intron_saf \
  -o ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_intron_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

featureCounts \
  -F SAF \
  -a $HG38_promoter_1500_500bp_saf \
  -o ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_promoter_1500_500bp_noOverlap_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam

cd ${outputDir}/12_reads_mapping_stats
python /home/yz2474/scripts/donglab/EVscope/scripts/plot_reads_mapping_stats.py \
  --sampleName ${SampleName} \
  --input_5UTR_readcounts ${SampleName}_HG38_5UTR_noOverlap_read_level_featureCounts.tsv \
  --input_exon_readcounts ${SampleName}_HG38_exon_noOverlap_read_level_featureCounts.tsv \
  --input_3UTR_readcounts ${SampleName}_HG38_3UTR_noOverlap_read_level_featureCounts.tsv \
  --input_intron_readcounts ${SampleName}_HG38_intron_noOverlap_read_level_featureCounts.tsv \
  --input_promoters_readcounts ${SampleName}_HG38_promoter_1500_500bp_noOverlap_read_level_featureCounts.tsv \
  --input_downstream_2Kb_readcounts ${SampleName}_HG38_downstream_2kb_noOverlap_read_level_featureCounts.tsv \
  --input_intergenic_readcounts ${SampleName}_HG38_intergenic_noOverlap_read_level_featureCounts.tsv \
  --input_ENCODE_blacklist_readcounts ${SampleName}_HG38_ENCODE_blacklist_V2_read_level_featureCounts.tsv
cd \$curr_dir
"

########################################
# Step 13: Kraken for Taxonomic Classification (Downsampling)
########################################
run_step "${outputDir}/13_Kraken" "
seqtk sample -s100 ${outputDir}/03_cutadapt/${SampleName}_R1_clean.fq.gz 100000 | gzip > ${outputDir}/13_Kraken/${SampleName}_R1_downsampled.fq.gz
seqtk sample -s100 ${outputDir}/03_cutadapt/${SampleName}_R2_clean.fq.gz 100000 | gzip > ${outputDir}/13_Kraken/${SampleName}_R2_downsampled.fq.gz
conda run -n kraken2_env kraken2 --db /home/yz2474/yiyong_2023/soft/kraken2/standard_krakendb --threads 20 \
  --report ${outputDir}/13_Kraken/${SampleName}_report.tsv --report-minimizer-data \
  --paired --gzip-compressed \
  ${outputDir}/13_Kraken/${SampleName}_R1_downsampled.fq.gz \
  ${outputDir}/13_Kraken/${SampleName}_R2_downsampled.fq.gz
python /home/yz2474/yiyong_2023/soft/KrakenTools/kreport2krona.py \
  -r ${outputDir}/13_Kraken/${SampleName}_report.tsv \
  -o ${outputDir}/13_Kraken/${SampleName}_krona_input.tsv
conda run -n kraken2_env ktImportText ${outputDir}/13_Kraken/${SampleName}_krona_input.tsv \
  -o ${outputDir}/13_Kraken/${SampleName}_krona.html
" true

########################################
# Step 14-1: bulk-RNA-seq deconvolution for featurecounts
########################################
run_step "${outputDir}/14_bulk_RNA_deconvolution/featureCounts" "
awk '{print \$1\",\" \$5}' ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_Gene_readcounts_normalized_expression_matrix_featureCounts.tsv > ${outputDir}/14_bulk_RNA_deconvolution/featureCounts/${SampleName}_Ensembl_ID_TPM_matrix_featureCounts.csv
cd ${outputDir}/14_bulk_RNA_deconvolution/featureCounts
python /home/yz2474/scripts/donglab/EVscope/scripts/14_run_RNA_deconvolution_ARIC.py --bulk_expr_file ${SampleName}_Ensembl_ID_TPM_matrix_featureCounts.csv  --ref_expr_file /home/yz2474/scripts/donglab/EVscope/genome_anno/reference_RNA_deconvolution/GTex_v10/GTEx_v10_human_SMTS_AVG_TPM.csv

python /home/yz2474/scripts/donglab/EVscope/scripts/14_run_RNA_deconvolution_ARIC.py --bulk_expr_file ${SampleName}_Ensembl_ID_TPM_matrix_featureCounts.csv  --ref_expr_file /home/yz2474/scripts/donglab/EVscope/genome_anno/reference_RNA_deconvolution/GTex_v10/GTEx_v10_huaman_SMTSD_AVG_TPM.csv

python /home/yz2474/scripts/donglab/EVscope/scripts/14_run_RNA_deconvolution_ARIC.py --bulk_expr_file ${SampleName}_Ensembl_ID_TPM_matrix_featureCounts.csv  --ref_expr_file /home/yz2474/scripts/donglab/EVscope/genome_anno/reference_RNA_deconvolution/Human_Brain_Single_Cell_Atlas_V1/human_brain_single_cell_atlas_V1_31_superclusters_AVP_CPM.csv
cd \$curr_dir
" true

########################################
# Step 14-2: bulk-RNA-seq deconvolution for RSEM
########################################
run_step "${outputDir}/14_bulk_RNA_deconvolution/RSEM" "
awk '{print \$1\",\" \$5}' ${outputDir}/11_RNA_distribution/RSEM/${SampleName}_Gene_readcounts_normalized_expression_matrix_RSEM.tsv > ${outputDir}/14_bulk_RNA_deconvolution/RSEM/${SampleName}_Ensembl_ID_TPM_matrix_RSEM.csv

cd ${outputDir}/14_bulk_RNA_deconvolution/RSEM
python /home/yz2474/scripts/donglab/EVscope/scripts/14_run_RNA_deconvolution_ARIC.py --bulk_expr_file ${SampleName}_Ensembl_ID_TPM_matrix_RSEM.csv  --ref_expr_file /home/yz2474/scripts/donglab/EVscope/genome_anno/reference_RNA_deconvolution/GTex_v10/GTEx_v10_human_SMTS_AVG_TPM.csv

python /home/yz2474/scripts/donglab/EVscope/scripts/14_run_RNA_deconvolution_ARIC.py --bulk_expr_file ${SampleName}_Ensembl_ID_TPM_matrix_RSEM.csv  --ref_expr_file /home/yz2474/scripts/donglab/EVscope/genome_anno/reference_RNA_deconvolution/GTex_v10/GTEx_v10_huaman_SMTSD_AVG_TPM.csv

python /home/yz2474/scripts/donglab/EVscope/scripts/14_run_RNA_deconvolution_ARIC.py --bulk_expr_file ${SampleName}_Ensembl_ID_TPM_matrix_RSEM.csv  --ref_expr_file /home/yz2474/scripts/donglab/EVscope/genome_anno/reference_RNA_deconvolution/Human_Brain_Single_Cell_Atlas_V1/human_brain_single_cell_atlas_V1_31_superclusters_AVP_CPM.csv
cd \$curr_dir
"

########################################
# Step 15: rRNA Detection
########################################
cd $curr_dir
run_step "${outputDir}/15_RiboDetector" "
ribodetector_cpu -t 20 \
  -l 100 \
  -i ${outputDir}/13_Kraken/${SampleName}_R1_downsampled.fq.gz  \
     ${outputDir}/13_Kraken/${SampleName}_R2_downsampled.fq.gz  \
  -e rrna \
  --chunk_size 800 \
  -r ${outputDir}/15_RiboDetector/${SampleName}_rRNA_R1.fq.gz \
     ${outputDir}/15_RiboDetector/${SampleName}_rRNA_R2.fq.gz \
  -o /dev/null /dev/null
"

########################################
# Step 16: QC_matrix table generation
########################################
run_step "${outputDir}/16_QC_matrix_generation" "
python /home/yz2474/scripts/donglab/EVscope/scripts/11_generate_QC_matrix.py \
    --raw_R1_fastqc_zip ${outputDir}/01_raw_fastqc/*${SampleName}*R1*_fastqc.zip \
    --raw_R2_fastqc_zip ${outputDir}/01_raw_fastqc/*${SampleName}*R2*_fastqc.zip \
    --trimmed_R1_fastq ${outputDir}/03_cutadapt/${SampleName}_R1_clean.fq.gz \
    --trimmed_R2_fastq ${outputDir}/03_cutadapt/${SampleName}_R2_clean.fq.gz \
    --ecoli_R1_fastq ${outputDir}/05_BBSplit/${SampleName}_E.coli_GCF_000005845.2_ASM584v2_ge_R1.fq.gz \
    --ecoli_R2_fastq ${outputDir}/05_BBSplit/${SampleName}_E.coli_GCF_000005845.2_ASM584v2_ge_R2.fq.gz \
    --myco_R1_fastq ${outputDir}/05_BBSplit/${SampleName}_Mycoplasma_240_ge_R1.fq.gz \
    --myco_R2_fastq ${outputDir}/05_BBSplit/${SampleName}_Mycoplasma_240_ge_R2.fq.gz \
    --umi_dedup_bam ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam \
    --bam2strand_file ${outputDir}/07_strand_detection/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out_bam2strandness.tsv \
    --STAR_log ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Log.final.out \
    --picard_insert_file ${outputDir}/09_picard_metrics/${SampleName}_STAR_UMI_dedup.insert_size_metrics.tsv \
    --picard_rnaseq_file ${outputDir}/09_picard_metrics/${SampleName}_STAR_UMI_dedup.picard_metrics.tsv \
    --featureCounts_3UTR ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_3UTR_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_5UTR ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_5UTR_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_downstream_2kb ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_downstream_2kb_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_exon ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_exon_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_ENCODE_blacklist ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_ENCODE_blacklist_V2_read_level_featureCounts.tsv \
    --featureCounts_intergenic ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_intergenic_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_intron ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_intron_noOverlap_read_level_featureCounts.tsv \
    --featureCounts_promoter_1500_500bp ${outputDir}/12_reads_mapping_stats/${SampleName}_HG38_promoter_1500_500bp_noOverlap_read_level_featureCounts.tsv \
    --expression_matrix ${outputDir}/11_RNA_distribution/featureCounts/${SampleName}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM.tsv \
    --kraken_report ${outputDir}/13_Kraken/${SampleName}_report.tsv \
    --downsampled_trimmed_R1_fastq  ${outputDir}/13_Kraken/${SampleName}_R1_downsampled.fq.gz \
    --downsampled_trimmed_R2_fastq  ${outputDir}/13_Kraken/${SampleName}_R2_downsampled.fq.gz \
    --ribo_R1_fastq ${outputDir}/15_RiboDetector/${SampleName}_rRNA_R1.fq.gz \
    --ribo_R2_fastq ${outputDir}/15_RiboDetector/${SampleName}_rRNA_R2.fq.gz \
    --output ${outputDir}/16_QC_matrix_generation/${SampleName}_QC_matrix.tsv
"

########################################
 #Step 17: Generate Output Report using R Markdown
########################################
run_step "${outputDir}/17_HTML_report" "
Rscript -e \"rmarkdown::render('/home/yz2474/scripts/donglab/EVscope/scripts/html_reprot.Rmd',
              params=list(output_dir='${outputDir}'),
              output_file=file.path('${outputDir}', '17_HTML_report', '${SampleName}_final_report.html'))\"
"


########################################
 #Step 18-1: bigWig/bedGraph Density Generation
########################################
#run_step "${outputDir}/18_bigWig/bamCoverage" "
#bamCoverage -b ${outputDir}/06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam \
#  -o ${outputDir}/18_bigWig/bamCoverage/${SampleName}_UMI_dedup_unstranded.bw \
#  --normalizeUsing BPM --numberOfProcessors 20
#bash /home/yz2474/scripts/donglab/EVscope/scripts/12_density_plot_over_meta_gene.sh ${outputDir}/18_bigWig/bamCoverage/${SampleName}_UMI_dedup_unstranded.bw
#bash /home/yz2474/scripts/donglab/EVscope/scripts/13_density_plot_over_RNA_types.sh ${outputDir}/18_bigWig/bamCoverage/${SampleName}_UMI_dedup_unstranded.bw
#"

########################################
 #Step 18-2: bigWig/bedGraph Density Generation
########################################
#run_step "${outputDir}/18_bigWig/EMapper" "
#cd ${outputDir}/18_bigWig/EMapper
#python /home/yz2474/scripts/donglab/EVscope/scripts/EMapper.py --input_bam ../../06_STAR/2nd_STAR/${SampleName}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam --sample_name ${SampleName} --bam_aligner STAR --split_by_strand yes --mode multi
#"

echo "Pipeline completed successfully for sample: ${SampleName}"

