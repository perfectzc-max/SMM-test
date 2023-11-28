# SMM-tset
try to build a pipeline for smm-seq
# Data processing and variant calling
Raw sequence reads were adapter- and quality-trimmed, aligned to human reference genome, realigned, and recalibrated on the basis of known indels as we described previously (7) except that deduplication step was omitted.
For variant calling, we developed a set of filters that were applied to each position in SMM-seq data. Only reads in proper pairs, with mapping quality not less than 60 and without secondary alignments, were taken in consideration. Positions in SMM-seq data were considered as qualified for variant calling if it is covered by UMI family containing not less than seven reads from each strand and this position is covered at least 20× in regular sequencing data. The qualified position was considered as a potential variant if all the reads within a given UMI family reported the same base at this position and this base was different from the corresponding reference genome. Next, to filter out germline variants, we checked if a found potential variant is in a list of SNPs of this DNA sample as well as in dbSNP. A list of sample specific germline SNPs was prepared by analysis of conventional sequencing data with Genome Analysis Toolkit (GATK) HaplotypeCaller. Last, a variant was rejected if one or more reads of a different UMI family in SMM data or in conventional data contained the same variant. SNV frequency was calculated as a ratio of the number of identified variants to the total number of qualified positions.
# step1. download data
Data may be accessed using the following link: https://dataview.ncbi.nlm.nih.gov/object/PRJNA758911.  
We use the SMM-seq raw data of ENU 50 sample and their control to build this pipeline.
### 1.Download data and tansfer data into fastq file.(download.sh)
### 2.Trim UMIs from reads in FASTQ files and make new ones with UMIs in the read names. Paired reads receive the same UMI (UMI1+UMI2). Automatically process files with read 1 and read 2.
(variants calling)
```
sbatch umi_trimmer.bash FileName_[1.fastq.gz] Prefix_len UMI_len Postfix_len
```

# step2. Alignment and quality trimming
### nf-core (unavailable)
### our gatk pipeline
This pipline includ fastQC, trimgalore, mapping and sort.
Use hapllotype caller call germline mutation, also preformed BQSR, bamQC.
In this step we omit the MarkDuplicate step in GATK pipeline for sample raw data, but still run it in bulkdata. So we perform Sarek pipeline for bulk data and another custom pipeline for sample data.
```
sbatch smm_gatk.batch <path to R1.fastq> <path to R2.fastq> <sample ID, suffix of fast1 file> <path to outdir>
```

### process summary：
Raw sequence reads were adapter and quality trimmed using Trim Galore (version 0.3.7), and aligned to human reference genome build 37 using BWA MEM (version 0.7.10)16. The initial mapped reads were indel realigned based on known indels from the 1000 Genomes Project (Phase I)17, and base quality score recalibrated based on known indels from the 1000 Genomes Project (Phase I) and SNVs from dbSNP (build 144) using GATK (version 3.4.46)8.


# step3. UMI summary and variants calling
## UMI summary
```
bash umi_qc.sh
```


## variants calling
### 1.Only reads in proper pairs, with mapping quality not less than 60 and without secondary alignments, were taken in consideration. 

### 2.Positions in SMM-seq data were considered as qualified for variant calling if it is covered by UMI family containing not less than seven reads from each strand and this position is covered at least 20× in regular sequencing data. 
UMI family containing not less than seven reads from each strand
covered at least 20× in regular sequencing data
### 3.The qualified position was considered as a potential variant if all the reads within a given UMI family reported the same base at this position and this base was different from the corresponding reference genome. 
all reads within a umi family reported same base
### 4.Next, to filter out germline variants, we checked if a found potential variant is in a list of SNPs of this DNA sample as well as in dbSNP. 
