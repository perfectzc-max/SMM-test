srun -w node06 --pty /bin/bash
module load sratoolkit/3.0.2
nohup prefetch --max-size 100G --option-file list.txt &
