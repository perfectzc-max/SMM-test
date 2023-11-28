#!/bin/bash
#SBATCH -p normal
#SBATCH -J SMM_plp
#SBATCH --mem=4GB
#SBATCH -o SMM_plp_%A.log
#SBATCH -t 2-00:00:00


source SMM_env.bash

echo "DNA mapping for human start on"; date; echo "====="

mkdir pileup

samplename=$1


#refgenome=/gs/gsfs0/users/amaslov/ref/homo-sapiens/b37/human_g1k_v37_decoy.fasta
#refgenome2=/gs/gsfs0/users/amaslov/ref/homo-sapiens/b37/samtools/human_g1k_v37_decoy.fasta
#refvar=/gs/gsfs0/users/amaslov/ref/homo-sapiens/gatk-b37


echo "samplename=$samplename;np=$np"
echo "The reference will be used: $refgenome"
echo "Mpileup started"; date
echo "samtools mpileup -f $refgenome $samplename > ./pileup/$samplename.mpu"
samtools mpileup -f $refgenome $samplename > ./pileup/$samplename.mpu

echo "Mpileup completed"; date
