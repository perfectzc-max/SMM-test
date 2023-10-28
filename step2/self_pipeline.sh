#!/bin/bash

# 设置环境变量
export GATK_PATH=/path/to/gatk
export TRIM_GALORE_PATH=/path/to/trim_galore
export FASTQC_PATH=/path/to/fastqc
export REFERENCE=reference.fasta
export DBSNP=dbsnp.vcf
export TARGET_BED=chrY_9055175-9057608.bed
export KNOWN_INDELS=Homo_sapiens_assembly38.known_indels.vcf.gz
export MILLS_GOLD_INDELS=Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# 设置Java虚拟机选项
JAVA_OPTS="-Xmx7g"

# 提取样本ID
INPUT_R1="id_R1.fastq.gz"
INPUT_R2="id_R2.fastq.gz"
SAMPLE_ID=$(echo "$INPUT_R1" | sed 's/_R1.*//')

# 运行FastQC以进行质量评估
$FASTQC_PATH -o fastqc_reports "$INPUT_R1" "$INPUT_R2"

# 运行Trim Galore以去除低质量的碱基和接头
$TRIM_GALORE_PATH --paired --output_dir trimmed_data "$INPUT_R1" "$INPUT_R2"

# 数据预处理：Fastq到BAM格式的映射，包括带有样本ID的bwa mem命令
bwa mem -K 100000000 -R "@RG\tID:1\tPU:1\tSM:$SAMPLE_ID\tLB:$SAMPLE_ID\tPL:illumina" -t 16 -M $REFERENCE trimmed_data/"${SAMPLE_ID}_R1_val_1.fq" trimmed_data/"${SAMPLE_ID}_R2_val_2.fq" | samtools view -bS - > mapped.bam

# 创建BAM文件索引
samtools index mapped.bam

# 标记重复序列
gatk --java-options "$JAVA_OPTS" MarkDuplicates --INPUT mapped.bam --METRICS_FILE "$SAMPLE_ID.bam.metrics" --TMP_DIR . --ASSUME_SORT_ORDER coordinate --CREATE_INDEX true --OUTPUT "$SAMPLE_ID.md.bam"

# 基本质量分数校正（BQSR）
gatk --java-options "$JAVA_OPTS" BaseRecalibrator -R $REFERENCE -I "$SAMPLE_ID.md.bam" -O "$SAMPLE_ID.recal_data.table" --known-sites $DBSNP --known-sites $KNOWN_INDELS --known-sites $MILLS_GOLD_INDELS -L $TARGET_BED --verbosity INFO

# 应用基本质量分数校正（BQSR），包括优化的Java虚拟机选项
gatk --java-options "$JAVA_OPTS" ApplyBQSR -R $REFERENCE --input "$SAMPLE_ID.md.bam" --output "$SAMPLE_ID.recalibrated.bam" -L $TARGET_BED --bqsr-recal-file "$SAMPLE_ID.recal_data.table"

# 变异调用（HaplotypeCaller）
$GATK_PATH HaplotypeCaller -R $REFERENCE -I "$SAMPLE_ID.recalibrated.bam" -O "$SAMPLE_ID.vcf"

# 现在添加 Genotyping 和 Indexing 步骤
# 创建gVCF文件的索引
$GATK_PATH IndexFeatureFile -I "$SAMPLE_ID.g.vcf"

# 对gVCF文件进行基因型调用，生成VCF文件
$GATK_PATH GenotypeGVCFs -R $REFERENCE -L $TARGET_BED --D $DBSNP -V "$SAMPLE_ID.g.vcf" -O "$SAMPLE_ID.vcf"

# 运行 bcftools stats 生成统计信息
bcftools stats merged_variants.vcf > merged_variants.bcf.tools.stats.out

# 使用 vcftools 进行统计
vcftools --gzvcf merged_variants.vcf --TsTv-by-count --out HaplotypeCaller_SRR15669403
vcftools --gzvcf merged_variants.vcf --TsTv-by-qual --out HaplotypeCaller_SRR15669403
vcftools --gzvcf merged_variants.vcf --FILTER-summary --out HaplotypeCaller_SRR15669403

# 使用 Qualimap 进行质量评估
qualimap --java-mem-size=128G bamqc -bam "$SAMPLE_ID.recalibrated.bam" --paint-chromosome-limits --genome-gc-distr HUMAN -nt 16 --skip-duplicated --skip-dup-mode 0 -outdir "$SAMPLE_ID.recal" -outformat HTML

# 过滤和注释变异
# 这一步通常需要使用额外的工具来过滤和注释变异，如GATK的VariantFiltration，Annovar，或其他工具，具体根据您的需求选择。
# 运行 VEP 进行变异注释
vep -i merged_variants.vcf -o "$SAMPLE_ID_VEP.ann.vcf" --assembly GRCh38 --species homo_sapiens --offline --cache --cache_version 99 --dir_cache /.vep --everything --filter_common --fork 4 --format vcf --per_gene --stats_file "$SAMPLE_ID_VEP.summary.html" --total_length --vcf

# 结果报告
# 生成包含变异信息的结果报告，通常以VCF格式。此处根据需求和工具自定义。
# 合并变异记录
$GATK_PATH GatherVcfs -I "$SAMPLE_ID.vcf" -O merged_variants.vcf

