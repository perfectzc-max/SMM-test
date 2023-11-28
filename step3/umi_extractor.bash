#!/bin/bash
#SBATCH -J umi
#SBATCH -p normal #partition
#SBATCH -o umi-%J.out
#SBATCH -e umi-%J.err

source activate py2
which python

bamdir=$1
Data=$2
Chr=$3
mkdir -p ./umi
python main.py umi $bamdir/$Data $Chr > ./umi/$Data.$Chr.umi
