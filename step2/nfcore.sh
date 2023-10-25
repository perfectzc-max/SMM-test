#creat samplesheet like samplesheet.tsv 
ls /cluster/home/qiangyu/smmdata/*_1*.gz|awk -F '[_/]' '{print $6  "\t" "NA" "\t" '0' "\t"  $6 "\t" '1' "\t" $0 "\t" $0 }'>samplelist.tsv
#vi change the filename to your real fastq filename
vi samplelist.tsv

#bulk
tail samplelist.tsv -n 3 > sample.tsv
#sample
head -n 1 samplelist.tsv > bulk.tsv

#help：/cluster/apps/nf_core/nf-core-sarek-2.7/workflow/main.nf
#run nf-core for bulk data
module load gatk/4.3.0
module load nfcore/sarek/2.7
mkdir bulk
sh run_nfcore_haplotypecaller_vep_bam.sh GRCh38 /cluster/groups/Jan-Lab/qiangyu/smmnfcore/bulk.tsv /cluster/groups/Jan-Lab/qiangyu/smmnfcore/bulk

#run without markdupulicate step for sample data
mkdir sample
sh run_nfcore_haplotypecaller_vep_nodup.sh GRCh38 /cluster/groups/Jan-Lab/qiangyu/smmnfcore/sample.tsv /cluster/groups/Jan-Lab/qiangyu/smmnfcore/sample

#if not working, we need build our own pipeline
