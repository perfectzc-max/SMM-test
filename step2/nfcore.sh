#creat samplesheet like samplesheet.tsv 
ls /cluster/home/qiangyu/smmdata/*R1*.gz|awk -F '[_/]' '{print $7  "\t" "NA" "\t" '0' "\t"  $7 "\t" '1' "\t" $0 "\t" $0 }'>samplelist.tsv

#run nf-core
module load module load nfcore/sarek/2.7
run_nfcore_haplotypecaller_vep.sh GRCh38 /cluster/home/qiangyu/smmdata/samplelist.tsv /cluster/home/qiangyu/bulkdata/nfcore/
