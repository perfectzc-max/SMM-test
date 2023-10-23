#!/bin/bash

### 
### This is a help document.
###
### Usage:
###   run_nfcore_haplotypecaller_vep.sh <refs> <inputfile> <outdir>
###
### Options:
###   <refs>	GRCh37 or GRCh38 or GRCm38.
###   <inputfile>	The absolute path to samplesheet.
###   <outdir>	Must be an absolute path.
###   -h	Help document.
###
### Path to this script: /cluster/apps/pipelines/nfcore/sarek/2.7
### Get more information of this pipeline: https://nf-co.re/sarek/2.7
###
###
###

help() {
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]
	then
		help
		exit
	else
		if [[ "$2" == "" ]]
			then
				refs=$1 && echo -e "\nRefs set to $refs. Inputfile and outdir not specified. \n"
				exit
		else
			if [[ "$3" == "" ]]
				then
					refs=$1 && inputfile=$2 && outdir=$(pwd)
					echo -e "\nRefs set to $refs, Inputfile set to $inputfile. \nOutdir not specified, set to current directory: $outdir. \n"
			else
				refs=$1 && inputfile=$2 && outdir=$3
				echo -e "\nRefs set to $refs, Inputfile to $inputfile. \nOutdir set to $outdir. \n"
			fi
		fi
fi



cd $outdir
export NXF_SINGULARITY_CACHEDIR=/cluster/apps/nf_core/singularity_imgs/sarek_v2.7

/cluster/apps/nf_core/NEXTFLOW/nextflow run /cluster/apps/nf_core/nf-core-sarek-2.7/workflow/ \
--input $inputfile \
--tools HaplotypeCaller,VEP \
--skip_tools markduplicates \
--genome $refs \
--igenomes_base /cluster/apps/Refs/references \
-profile singularity \
-c /cluster/apps/nf_core/slurm_config/slurm_submit.config \
-bg > nfcore_haplotypecaller_vep.log 2>&1 && echo -e "\nJob is submitted. Please run 'squeue' to check the status. \n"
