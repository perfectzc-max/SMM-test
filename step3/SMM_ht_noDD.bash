#!/bin/bash
#SBATCH -p normal
#SBATCH -J SMM_ht
#SBATCH -n 1
#SBATCH --mem=30GB
#SBATCH -o SMM_ht_%A.log
#SBATCH -t 1-00:00:00

source ./SMM_env.bash

echo "HaplotypeCaller start on"; date; echo "====="

#Provide  sample name+ 'gatk'  as Data
Data=$1
Chr=$2

mkdir ./ht


#Remove duplicates, sort and index BAM file
#python main.py rmdup $Data.bam $Chr
cd ./ht
#mv ../$Data.$Chr.rmdup.bam .

ln -s ../$Data.bam $Data.$Chr.rmdup.bam 
ln -s ../$Data.bam.bai $Data.$Chr.rmdup.bam.bai 

#samtools sort -o $Data.$Chr.rmdup.bam $Data.$Chr.rmdup.bam
#samtools index $Data.$Chr.rmdup.bam    


#Find germline variants
java -Xmx20g -jar $gatk \
-R $refgenome \
-T HaplotypeCaller \
-I $Data.$Chr.rmdup.bam \
--dbsnp $dbsnp \
-stand_call_conf 30 \
-stand_emit_conf 10 \
-L $Chr \
-o $Data.$Chr.vcf  


#Sort, compress and index VCF file
bcftools sort -Oz -o $Data.$Chr.vcf.gz $Data.$Chr.vcf
tabix $Data.$Chr.vcf.gz 


#java -Xmx20g -jar $gatk \
#-R $refgenome \
#-T HaplotypeCaller \
#-I $Data.bam \
#--dbsnp $dbsnp \
#-stand_call_conf 30 \
#-stand_emit_conf 10 \
#-L $Chr \
#-o $Data.ht.$Chr.vcf
