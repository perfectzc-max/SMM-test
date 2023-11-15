#!/bin/bash
#SBATCH -J smm #job name
#SBATCH -p normal #partition
#SBATCH -n 3 #ntasks
#SBATCH -c 15 #cpus per task
#SBATCH -o smm-%J.out
#SBATCH -e smm-%J.err
#SBATCH -w node03

BULK=$1
# Extract sample ID, pattern "id_R2.fastq.gz"
INPUT_R1=$2
INPUT_R2=$3
SAMPLE_ID=$4
OUT_DIR=$5

echo -e "\nout:$OUT_DIR \nsample: $SAMPLE_ID \n BULK SAMPLE:$BULK"

# Set environment variables
module load samtools/1.17
module load gatk/4.3.0
module load anaconda3/2022.10
module load java/19.0.2
source activate /cluster/home/qiangyu/.conda/envs/nf-core/envs/trim
export GATK_PATH=/cluster/apps/gatk/4.3.0/gatk
export PICARD_PATH=/cluster/apps/picard/3.0.0/picard.jar
export TRIM_GALORE_PATH=/cluster/home/qiangyu/.conda/envs/nf-core/envs/trim/bin/trim_galore
export FASTQC_PATH=/cluster/home/qiangyu/.conda/envs/nf-core/envs/trim/bin/fastqc
export CUTADAPT_PATH=/cluster/home/qiangyu/.conda/envs/nf-core/envs/trim/bin/cutadapt
export BWA_PATH=/cluster/apps/bwa/0.7.17/bwa
export SAMTOOLS_PATH=/cluster/apps/samtools/1.17/bin/samtools
export BCFTOOLS_PATH=/cluster/apps/bcftools/1.17/bin/bcftools
export VCFTOOLS_PATH=/cluster/apps/vcftools/0.1.16/src/cpp/vcftools
export VEP_PATH=/cluster/apps/vep/ensembl-vep-release-109/vep
#reference dir
export REFERENCE=/cluster/groups/Jan-Lab/qiangyu/myref/Homo_sapiens_assembly38.fasta
export DBSNP=/cluster/groups/Jan-Lab/qiangyu/myref/dbsnp_146.hg38.vcf.gz
export KNOWN_INDELS=/cluster/groups/Jan-Lab/qiangyu/myref/Homo_sapiens_assembly38.known_indels.vcf.gz
export MILLS_GOLD_INDELS=/cluster/groups/Jan-Lab/qiangyu/myref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

echo "DNA mapping for human start on"; date; echo "====="
# Set Java virtual machine options
JAVA_OPTS="-Xmx7g"

# Create directories for results
mkdir -p $OUT_DIR
cd $OUT_DIR
mkdir ./processing/1-fastqc
mkdir ./processing/2-trim
mkdir ./processing/3-fastqc
mkdir ./processing/4-bwa

# 1.Run FastQC for quality assessment
echo "1th fastqc start on"; date; echo "====="
$FASTQC_PATH -o ./processing/1-fastqc "$INPUT_R1" "$INPUT_R2"
echo "====="; echo "1th fastqc end on"; date

# 2.Run Trim Galore to remove low-quality bases and adapters
echo "2th trim galore start on"; date; echo "====="
$TRIM_GALORE_PATH --illumina --paired --gzip --output_dir ./processing/2-trim "$INPUT_R1" "$INPUT_R2" --path_to_cutadapt $CUTADAPT_PATH
echo "====="; echo "2th trim galore end on"; date

# 3.Run FastQC for quality assessment
echo "3th fastqc start on"; date; echo "====="
$FASTQC_PATH -o ./processing/1-fastqc ./processing/2-trim/"${SAMPLE_ID}_1_val_1.fq.gz" ./processing/2-trim/"${SAMPLE_ID}_2_val_2.fq.gz"
echo "====="; echo "3th fastqc end on"; date

# 4.Data preprocessing: Mapping from Fastq to BAM
echo "4th bwa start on"; date; echo "====="
$BWA_PATH mem -K 100000000 -R '@RG\tID:"$SAMPLE_ID"\tPL:"$SAMPLE_ID"\tPU:0\tLB:"$SAMPLE_ID"\tSM:"$SAMPLE_ID"' -t 16 -M $REFERENCE ./processing/2-trim/"${SAMPLE_ID}_1_val_1.fq.gz" ./processing/2-trim/"${SAMPLE_ID}_2_val_2.fq.gz" | $SAMTOOLS_PATH view -bS - > ./processing/4-bwa/mapped.bam

# Create an index for the mapped BAM file
$SAMTOOLS_PATH sort -@ 10 -m 4G -o ./processing/4-bwa/sorted.bam ./processing/4-bwa/mapped.bam
$SAMTOOLS_PATH index ./processing/4-bwa/sorted.bam
echo "====="; echo "4th bwa end on"; date

# 5.Picard
cd processing/4-bwa
mkdir tmp
chmod 777 tmp

java -Xmx20g -jar $PICARD_PATH CleanSam \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=`pwd`/tmp \
I=sorted.bam \
O=$SAMPLE_ID-pe.bam

java -Xmx20g -jar $PICARD_PATH AddOrReplaceReadGroups \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=`pwd`/tmp \
I=$SAMPLE_ID-pe.bam \
O=$SAMPLE_ID-pe.rd.bam \
RGID=$SAMPLE_ID RGLB=$SAMPLE_ID RGPL=illumina RGPU=$SAMPLE_ID RGSM=$SAMPLE_ID

java -Xmx20g -jar $PICARD_PATH FixMateInformation \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=`pwd`/tmp \
I=$SAMPLE_ID-pe.rd.bam \
O=$SAMPLE_ID-pe.rd.fixed.bam \
SO=coordinate

java -Xmx20g -jar $PICARD_PATH SortSam \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=`pwd`/tmp \
I=$SAMPLE_ID-pe.rd.fixed.bam \
O=$SAMPLE_ID-pe.rd.fixed.sorted.bam \
SO=coordinate \
CREATE_INDEX=true

java -Xmx20g -jar $PICARD_PATH ReorderSam \
TMP_DIR=`pwd`/tmp \
I=$SAMPLE_ID-pe.rd.fixed.sorted.bam \
O=$SAMPLE_ID-pe.rd.fixed.sorted.reorder.bam \
CREATE_INDEX=true \
REFERENCE=$REFERENCE

# Conditional execution of MarkDuplicates
if [ "$BULK" == "T" ]; then
  $GATK_PATH MarkDuplicates --INPUT sorted.bam --METRICS_FILE "$SAMPLE_ID.bam.metrics" --TMP_DIR . --ASSUME_SORT_ORDER coordinate --CREATE_INDEX true --OUTPUT "$SAMPLE_ID.md.bam"
else
  # If mark is "F", skip the MarkDuplicates step
  ln -s sorted.bam "$SAMPLE_ID.md.bam"
fi

# Base quality score recalibration (BQSR)
$GATK_PATH BaseRecalibrator -R $REFERENCE -I "$SAMPLE_ID.md.bam" -O "$SAMPLE_ID.recal_data.table" --known-sites $DBSNP --known-sites $KNOWN_INDELS --known-sites $MILLS_GOLD_INDELS --verbosity INFO

# Apply BQSR, including optimized Java virtual machine options
$GATK_PATH ApplyBQSR -R $REFERENCE --input "$SAMPLE_ID.md.bam" --output "$SAMPLE_ID.recalibrated.bam"  --bqsr-recal-file "$SAMPLE_ID.recal_data.table"

if [ "$BULK" == "T" ]; then
# Variant calling (HaplotypeCaller)
 $GATK_PATH HaplotypeCaller -R $REFERENCE -I "$SAMPLE_ID.recalibrated.bam" -O "$SAMPLE_ID.g.vcf"

# Index the gVCF file
 $GATK_PATH IndexFeatureFile -I "$SAMPLE_ID.g.vcf"

# Genotype gVCF files to create a VCF file
 $GATK_PATH GenotypeGVCFs -R $REFERENCE --D $DBSNP -V "$SAMPLE_ID.g.vcf" -O "$SAMPLE_ID.vcf"

# Merge variant records
 $GATK_PATH GatherVcfs -I "$SAMPLE_ID.vcf" -O merged_variants.vcf

# Run bcftools stats to generate statistics
 $BCFTOOLS_PATH stats merged_variants.vcf > merged_variants.bcftools.stats.out

# Use vcftools for statistics
 $VCFTOOLS_PATH --gzvcf merged_variants.vcf --TsTv-by-count --out HaplotypeCaller_SRR15669403
 $VCFTOOLS_PATH --gzvcf merged_variants.vcf --TsTv-by-qual --out HaplotypeCaller_SRR15669403
 $VCFTOOLS_PATH --gzvcf merged_variants.vcf --FILTER-summary --out HaplotypeCaller_SRR15669403

# Perform quality assessment using Qualimap
# $QUALIMAP_PATH --java-mem-size=128G bamqc -bam "$SAMPLE_ID.recalibrated.bam" --paint-chromosome-limits --genome-gc-distr HUMAN -nt 16 --skip-duplicated --skip-dup-mode 0 -outdir "$SAMPLE_ID.recal" -outformat HTML

# Filter and annotate variants
# This step typically requires additional tools such as GATK's VariantFiltration, Annovar, or others, depending on your specific needs.
# Run VEP for variant annotation
 $VEP_PATH -i merged_variants.vcf -o "$SAMPLE_ID_VEP.ann.vcf" --assembly GRCh38 --species homo_sapiens --offline --cache --cache_version 99 --dir_cache /.vep --everything --filter_common --fork 4 --format vcf --per_gene --stats_file "$SAMPLE_ID_VEP.summary.html" --total_length --vcf
 echo "Bulk sample finished."
else
 echo "Sample pipeline finished."
fi

