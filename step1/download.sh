srun -w node06 --pty /bin/bash
module load sratoolkit/3.0.2

#download sra file
nohup prefetch --max-size 100G --option-file list.txt &

#transfer into fastq
for i in *sra; do fastq-dump ${i} --split-3 --gzip -O ./; done
