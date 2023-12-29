# https://github.com/OpenGene/fastp#unique-molecular-identifier-umi-processing
#!/bin/bash
#SBATCH -t 08:00:00
#SBATCH -p normal
#SBATCH -J umi
#SBATCH -o umi_%A.log
#SBATCH -n 1
#SBATCH --mem=4GB

File=$1
Pre=$2
UMI=$3
Post=$4

source activate py2
python UMItrimmer.py $File $Pre $UMI $Post
