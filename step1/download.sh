srun -w node06 --pty /bin/bash
module load sratoolkit/3.0.2

#download sra file
nohup prefetch --max-size 100G --option-file list.txt &
#or (faster)
awk '{print "nohup prefetch --max-size 100G " $0 " &"}' list.txt|bash

#transfer into fastq
for i in *sra; do fastq-dump ${i} --split-3 --gzip -O ./; done
#OR 
awk '{print "nohup fastq-dump --split-3 --gzip /cluster/home/qiangyu/smmdata/" $1 "/" $1 ".sra -O ./ &"}' list.txt |bash

#submit sbatch job
awk '{print "nohup fastq-dump --split-3 --gzip /cluster/home/qiangyu/smmdata/" $1 "/" $1 ".sra -O ./ &"}' list.txt> code.txt
cat sbach.head code.sh >job.slurm
sbatch job.slurm
