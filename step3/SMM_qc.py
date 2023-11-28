import os
import sys
import subprocess
import Bio
from Bio import pairwise2
import pysam
import time

def umi_analysis(FileName):
    umi_file=FileName
    count_all=0
    count_good=0
    #count_intervals=0
    count_reads=0
    length_good=0
    duplex=0
    min_family_size=7
    f_input=open(umi_file)
    f_input.readline()
    current_interval=""
    spectra={}
    for line in f_input:
        fields=line.split()
        count_all+=1

        try:
            fields[2]=int(fields[2])
            fields[3]=int(fields[3])
            count_reads+=fields[2]+fields[3]
            
            if min(fields[2], fields[3])>0:
                duplex+=1

                if int(fields[2])>=min_family_size and int(fields[3])>=min_family_size: 
                    count_good+=1
                length_good+=int(fields[2])+int(fields[3])
                family_size=int(fields[2])+int(fields[3])
                count_reads+=family_size
                if family_size in spectra:
                    spectra[family_size]+=1
                else:
                    spectra[family_size]=1
        except: pass
    
    classes=[0,0,0,0]
    for item in spectra:
        if item<10: classes[0]+=spectra[item]
        elif item<20: classes[1]+=spectra[item]
        elif item<30: classes[2]+=spectra[item]
        else: classes[3]+=spectra[item]
        #print item, spectra[item]
    print 'Sample, Reads, Families, DuplexFamilies (%s), GoodFam >=%s, PrcGoodFam, <10, <20, <30, >30'  %('%', min_family_size)   
    print "%s,%s,%s,%s (%2.2f),%s,%2.2f,%s" %(umi_file, count_reads, count_all, duplex, duplex*100.0/count_all, count_good, count_good*100.0/count_all, ",".join(map(lambda item: str(item*100.0/sum(classes)), classes)))

    
def plp_analysis(FileName):     #Analysis og mpu file after pileup
    classes=[0,0,0]
    for line in open(FileName):
        if int(line.split()[3])<10:
            classes[0]+=1
        elif int(line.split()[3])<20:
            classes[1]+=1
        else:
            classes[2]+=1
        total=sum(classes)
        if total>1000000:
            break
        
    print '<10, >=10, >=20'
    print '%2.2f,%2.2f,%2.2f' %(classes[0]*100.0/total, (classes[1]+classes[2])*100.0/total, classes[2]*100.0/total)    

def read_mpu_file(MpuFile): #Read MPU file and return as a dictionary
    mpu={}
    for line in open(MpuFile):
        Chr=line.split()[0]
        Pos=int(line.split()[1])
        Cov=int(line.split()[3])
        if Chr in mpu:
            mpu[Chr][Pos]=Cov
        else:
            mpu[Chr]={Pos:Cov}
    return mpu
    


def main():
    if sys.argv[1]=='umi':
        umi_analysis(sys.argv[2])
    elif sys.argv[1]=='plp':
        plp_analysis(sys.argv[2])
        
    else:
        sys.exit('Unknown command:\t%s')  %sys.argv[1]
    
        

#--------------------------------
if __name__ == "__main__":
    main()

