### 1.umi count for S-lib
umi_extractor.bash bamdir S-file.bam chr1
python SMM_qc.py umi S-file.1.bam.umi

### 2.pileup sum for G_lib
sbatch SMM_pileup.bash bamdir G-file.bam
python SMM_qc.py plp G-file.mpu
