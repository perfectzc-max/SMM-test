import os
import sys
import subprocess
import Bio
from Bio import pairwise2
import pysam
import time

def read_fq(fq):
    try:
        lines = [fq.readline(), fq.readline(), fq.readline(), fq.readline()]
        if all(line == '' for line in lines):  
            return []
        return lines
    except:
        return []


def write_fq(fq, read_name, read, l_umi):
	fq.write('%s %s\n' % (read_name,read[0].split(' ')[1]))
	fq.write(read[1][l_umi:])
	fq.write(read[2])
	fq.write(read[3][l_umi:])
	
print sys.argv[1]
def main():
    l_pre=int(sys.argv[2])
    l_umi=int(sys.argv[3])
    l_post=int(sys.argv[4])       
    in_data1 = gzip.open(sys.argv[1] + '1.fastq.gz', 'r')
    in_data2 = gzip.open(sys.argv[1] + '2.fastq.gz', 'r')
    out_data1 = gzip.open(sys.argv[1] + 'umi_1.fastq.gz', 'w')
    out_data2 = gzip.open(sys.argv[1] + 'umi_2.fastq.gz', 'w')
	reads_total=0

	while True:
		read1=read_fq(in_data1)
		read2=read_fq(in_data2)
		if read1!=[] and read2!=[]:
			reads_total+=1
			if read1[0].split(' ')[0]==read2[0].split(' ')[0]:
                            read_name='%s:%s+%s' %(read1[0].split(' ')[0], read1[1][l_pre:l_pre+l_umi], read2[1][l_pre:l_pre+l_umi])
			    write_fq(out_data1, read_name, read1, l_pre+l_umi+l_post)
			    write_fq(out_data2, read_name, read2, l_pre+l_umi+l_post)
				
			else:
			    sys.exit("Name conflict:\t%s,%s" %(read1[0], read2[0]))
				
		else:
			in_data1.close()
			in_data2.close()
			out_data1.close()
			out_data2.close()
			print 'Processed:\t', sys.argv[2] 
			print 'Trimmed: %s symbols' %(sys.argv[1])
			print 'Total reads: ', reads_total
			break				
	


#--------------------------------
if __name__ == "__main__":
    main()
