module load FastQC/0.11.9 
module load bwa/0.7.17
module load bcftools/1.17
module load tabix/0.2.5
module load samtools/1.17
module load gatk/4.3.0
module load anaconda3/2022.10
module load java/19.0.2
# source activate /cluster/home/qiangyu/.conda/envs/nf-core/envs/trim
source activate /cluster/home/qiangyu/.conda/envs/py2

picard=/cluster/apps/picard/3.0.0/picard.jar
gatk=/cluster/apps/gatk/4.3.0/gatk
python2=/cluster/home/qiangyu/.conda/envs/py2/bin/python
## module load python/2.7.8/gcc.4.4.7

which python2

refgenome=/cluster/groups/Jan-Lab/qiangyu/myref/Homo_sapiens_assembly38.fasta
refvar=/cluster/groups/Jan-Lab/qiangyu/myref
dbsnp=/cluster/groups/Jan-Lab/qiangyu/myref/dbsnp_146.hg38.vcf
dbsnp_gz=/cluster/groups/Jan-Lab/qiangyu/myref/dbsnp_146.hg38.vcf.gz

chr_file=/cluster/groups/Jan-Lab/qiangyu/myref/chr_human.txt

