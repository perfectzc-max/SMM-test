#!/bin/bash
mark=$1
# Extract sample ID, pattern "id_R2.fastq.gz"
INPUT_R1=$2
INPUT_R2=$3
SAMPLE_ID=$(echo "$INPUT_R1" | sed 's/_R1.*//')

# Set environment variables
export GATK_PATH=/path/to/gatk
export TRIM_GALORE_PATH=/path/to/trim_galore
export FASTQC_PATH=/path/to/fastqc
export REFERENCE=reference.fasta
export DBSNP=dbsnp.vcf
export TARGET_BED=chrY_9055175-9057608.bed
export KNOWN_INDELS=Homo_sapiens_assembly38.known_indels.vcf.gz
export MILLS_GOLD_INDELS=Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

echo 

# Set Java virtual machine options
JAVA_OPTS="-Xmx7g"

# Create directories for results
mkdir -p fastqc_reports
mkdir -p trimmed_data

# Run FastQC for quality assessment
$FASTQC_PATH -o fastqc_reports "$INPUT_R1" "$INPUT_R2"

# Run Trim Galore to remove low-quality bases and adapters
$TRIM_GALORE_PATH --paired --output_dir trimmed_data "$INPUT_R1" "$INPUT_R2"

# Data preprocessing: Mapping from Fastq to BAM
bwa mem -K 100000000 -R "@RG\tID:1\tPU:1\tSM:$SAMPLE_ID\tLB:$SAMPLE_ID\tPL:illumina" -t 16 -M $REFERENCE trimmed_data/"${SAMPLE_ID}_R1_val_1.fq" trimmed_data/"${SAMPLE_ID}_R2_val_2.fq" | samtools view -bS - > mapped.bam

# Create an index for the mapped BAM file
samtools index mapped.bam

# Conditional execution of MarkDuplicates
if [ "$mark" == "T" ]; then
  $GATK_PATH MarkDuplicates --INPUT mapped.bam --METRICS_FILE "$SAMPLE_ID.bam.metrics" --TMP_DIR . --ASSUME_SORT_ORDER coordinate --CREATE_INDEX true --OUTPUT "$SAMPLE_ID.md.bam"
else
  # If mark is "F", skip the MarkDuplicates step
  mv mapped.bam "$SAMPLE_ID.md.bam"
fi

# Base quality score recalibration (BQSR)
$GATK_PATH BaseRecalibrator -R $REFERENCE -I "$SAMPLE_ID.md.bam" -O "$SAMPLE_ID.recal_data.table" --known-sites $DBSNP --known-sites $KNOWN_INDELS --known-sites $MILLS_GOLD_INDELS -L $TARGET_BED --verbosity INFO

# Apply BQSR, including optimized Java virtual machine options
$GATK_PATH ApplyBQSR -R $REFERENCE --input "$SAMPLE_ID.md.bam" --output "$SAMPLE_ID.recalibrated.bam" -L $TARGET_BED --bqsr-recal-file "$SAMPLE_ID.recal_data.table"

# Variant calling (HaplotypeCaller)
$GATK_PATH HaplotypeCaller -R $REFERENCE -I "$SAMPLE_ID.recalibrated.bam" -O "$SAMPLE_ID.g.vcf"

# Index the gVCF file
$GATK_PATH IndexFeatureFile -I "$SAMPLE_ID.g.vcf"

# Genotype gVCF files to create a VCF file
$GATK_PATH GenotypeGVCFs -R $REFERENCE -L $TARGET_BED --D $DBSNP -V "$SAMPLE_ID.g.vcf" -O "$SAMPLE_ID.vcf"

# Merge variant records
$GATK_PATH GatherVcfs -I "$SAMPLE_ID.vcf" -O merged_variants.vcf

# Run bcftools stats to generate statistics
bcftools stats merged_variants.vcf > merged_variants.bcftools.stats.out

# Use vcftools for statistics
vcftools --gzvcf merged_variants.vcf --TsTv-by-count --out HaplotypeCaller_SRR15669403
vcftools --gzvcf merged_variants.vcf --TsTv-by-qual --out HaplotypeCaller_SRR15669403
vcftools --gzvcf merged_variants.vcf --FILTER-summary --out HaplotypeCaller_SRR15669403

# Perform quality assessment using Qualimap
qualimap --java-mem-size=128G bamqc -bam "$SAMPLE_ID.recalibrated.bam" --paint-chromosome-limits --genome-gc-distr HUMAN -nt 16 --skip-duplicated --skip-dup-mode 0 -outdir "$SAMPLE_ID.recal" -outformat HTML

# Filter and annotate variants
# This step typically requires additional tools such as GATK's VariantFiltration, Annovar, or others, depending on your specific needs.
# Run VEP for variant annotation
vep -i merged_variants.vcf -o "$SAMPLE_ID_VEP.ann.vcf" --assembly GRCh38 --species homo_sapiens --offline --cache --cache_version 99 --dir_cache /.vep --everything --filter_common --fork 4 --format vcf --per_gene --stats_file "$SAMPLE_ID_VEP.summary.html" --total_length --vcf