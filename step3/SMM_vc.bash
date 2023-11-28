#!/bin/bash
#SBATCH -p ht
#SBATCH -J SMM_vc
#SBATCH -n 1
#SBATCH --mem=12GB
#SBATCH -o SMM_vc_%A.log
#SBATCH -t 1-00:00:00

#MODULEPATH=/gs/gsfs0/users/vijg-lab/2022-SlurmSoftware/modulefiles:$MODULEPATH

##module load Python/2.7.18-Alex
#module load Python/2.7.18-control

source activate py2

python main.py vc $1 $2 $3 $4 $5 $6

