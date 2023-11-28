import os
import sys
import subprocess
import Bio
from Bio import pairwise2
import pysam
import time


def start_analysis(SampleName, ChrList):
	for chr in ChrList:
		string='qsub LSDS.bash %s %s'	%(SampleName, chr)
		print string
		subprocess.call(string, shell=True)
		
		
def convert_VCF_string(VCF_string):	#Make a list from string describing VCF
	#print VCF_string
	fields=VCF_string.split('\t')
	fields[0]=fields[0]
	fields[1]=int(fields[1])
	try:
		fields[-1]=dict(map(lambda x: x.split('=') ,fields[-1].split(';')))		
	finally:		
		return fields	
		

def convert_VCF_list(VCF_list):		#Make a string from list describing VCF
	#print VCF_list
	VCF_list[0]=str(VCF_list[0])
	VCF_list[1]=str(VCF_list[1])
	VCF_list[7]=';'.join(map(lambda x: '='.join(x), VCF_list[7].items()))
	return '\t'.join(VCF_list)+'\n'		
	
	

def read_vcf_file(FileName):
	stat={}
	variants=[]
	par=''
	fi=open(FileName)
	if FileName.split('.')[-1]!='vcf':
		sys.exit('%s is not VCF file' %FileName)
		
	#print 'Reading %s' %FileName
	for rec in fi:
		
		if rec[:2]!='##':
			variants.append(convert_VCF_string(rec.strip()))	
			
		else:
			if 'Stat' in rec:
				stat=dict(map(lambda x: [x.split('=')[0], int(x.split('=')[1])], rec.split(':')[1].split(';')))
			elif 'Param' in rec:
				par=rec
				#par[chr]=dict(map(lambda x: x.split('='), rec.split(':')[1].split(';')))

	return par, stat, variants
			
	
def save_vcf_file(FileName, par, summary_stat, variants):
	print 'Writing %s'	%FileName
	fo=open(FileName, 'w')
	
	fo.write(par)
	fo.write('##Stat:'+';'.join(map(lambda x: '='.join([x[0], str(x[1])]), summary_stat.items()))+'\n')
	
	for var in variants:
		fo.write(convert_VCF_list(var))
		
	fo.close()
	

def combine(SampleName, ChrList):
	stat={}
	variants=[]
	summary_stat={}
	par=[]
	#par={}
	for chr in ChrList:
		ChrFile='%s.%s.vcf' %(SampleName.split('.bam')[0], chr)
		print ChrFile, '\t', 
		try:
			chr_par, chr_stat, chr_variants=read_vcf_file(ChrFile)
			for variable in chr_stat:
				if chr_stat[variable]==0:
					print 'File is corrupted'
					continue					
			par.append(chr_par)
			stat[chr]=chr_stat
			variants+=chr_variants
			print 'OK'

		except:
			print 'File does not exist'

	if len(set(par))!=1: 
		print 'Error: Parameters are not consistent for sample %s.'	%SampleName
		print set(par)
		return
	
	for chr in stat:
		for variable in stat[chr]:
			if variable in summary_stat:
				summary_stat[variable]+=stat[chr][variable]
			else:
				summary_stat[variable]=stat[chr][variable]
	
	OutputFileName='%s.vcf' %SampleName.split('.bam')[0]	
	
	save_vcf_file(OutputFileName, par[0], summary_stat, variants)
	#sys.exit()
	
	
	'''
	print summary_stat	
	grm_v=len(filter(lambda x: x[7]['SAO']=='1', variants))
	som_v=len(filter(lambda x: x[7]['SAO']=='2', variants))
	snv_v=len(filter(lambda x: x[7]['SAO']=='2' and x[7]['VC']=='SNV', variants))
	ins_v=len(filter(lambda x: x[7]['SAO']=='2' and x[7]['VC']=='INS', variants))
	del_v=len(filter(lambda x: x[7]['SAO']=='2' and x[7]['VC']=='DEL', variants))
	
	print(summary_stat['TBS'], summary_stat['QBS'], summary_stat['FBS']*1.0/summary_stat['QBS'], summary_stat['SBS'], grm_v, som_v, snv_v, ins_v, del_v, 
										snv_v*1000000.0/summary_stat['SBS'],
										ins_v*1000000.0/summary_stat['SBS'],
										del_v*1000000.0/summary_stat['SBS'])
 	
		
	'''	


		
		
def is_same_variant(var1, var2):
	if var1[0]==var2[0] and int(var1[1])==int(var2[1]) and var1[3]==var2[3] and var1[4]==var2[4]:
		return True
	else:
		return False




def get_bulk_bases(chr, pos, aln_data):		#get list of calls from bulk 
	bases=[]
	UMIs=[]
	for pileupcolumn in aln_data.pileup(chr, pos-1, pos, truncate=True):
		for pileupread in pileupcolumn.pileups:
			if pileupcolumn.pos==pos-1 and not pileupread.is_del and not pileupread.is_refskip:
				if pileupread.alignment.query_name.split(':')[-1] not in UMIs:
					UMIs.append(pileupread.alignment.query_name.split(':')[-1])
					bases.append(pileupread.alignment.query_sequence[pileupread.query_position])
	return bases	


'''		

def get_bulk_bases(chr, pos, aln_data):		#get list of calls from bulk 
	bases=[]
	for pileupcolumn in aln_data.pileup(chr, pos-1, pos, truncate=True):
		for pileupread in pileupcolumn.pileups:
			if pileupcolumn.pos==pos-1 and not pileupread.is_del and not pileupread.is_refskip:
				bases.append(pileupread.alignment.query_sequence[pileupread.query_position])
				#self.pileup_array.append([pileupread.alignment, pileupread.alignment.query_sequence[pileupread.query_position], pileupread.query_position])
	return bases			
				
		
def germline(SampleName, BulkNames):
	
    VcfFileName='%svc7.vcf' %SampleName.split('bam')[0]
    par, summary_stat, variants=read_vcf_file(VcfFileName)		#load identified variants
	
    print 'Marking germline SNVs'
    for BulkName in open(BulkNames):
        if BulkName.strip()!=SampleName:
            aln_data=pysam.AlignmentFile(BulkName.strip(), 'rb')
            print '\tSource:\t%s' %BulkName.strip()
            print '\tSomatic SNVs found before filtering:\t%s' %len(filter(lambda x: x[7]['SAO']=='2' and x[7]['VC']=='SNV', variants)) 
            for i in range(len(variants)):
                if variants[i][7]['SAO']=='2' and variants[i][7]['VC']=='SNV':
                    bases=get_bulk_bases(variants[i][0], variants[i][1], aln_data)
                    if variants[i][4] in bases:
                        alt_count=bases.count(variants[i][4])
                        #if float(alt_count)/len(bases)>0.1:
                        if alt_count>2:
                            variants[i][7]['SAO']='1'
            print '\tSomatic SNVs found after filtering:\t%s' %len(filter(lambda x: x[7]['SAO']=='2' and x[7]['VC']=='SNV', variants))
	
    OutputFileName='%sgrmrm.vcf' %VcfFileName.split('vcf')[0]
    save_vcf_file(OutputFileName, par, summary_stat, variants)				


def germline(SampleName, BulkName):
	
	VcfFileName='%svc7.vcf' %SampleName.split('bam')[0]
	par, summary_stat, variants=read_vcf_file(VcfFileName)		#load identified variants
	
	print 'Marking germline SNVs'
	#for BulkName in open(BulkNames):
		#if BulkName.strip()!=SampleName:
	aln_data=pysam.AlignmentFile(BulkName.strip(), 'rb')
	print '\tSource:\t%s' %BulkName.strip()
	print '\tSomatic SNVs found before filtering:\t%s' %len(filter(lambda x: x[7]['SAO']=='2' and x[7]['VC']=='SNV', variants)) 
	for i in range(len(variants)):
		if variants[i][7]['SAO']=='2' and variants[i][7]['VC']=='SNV':
			bases=get_bulk_bases(variants[i][0], variants[i][1], aln_data)
			if len(bases)>=10:
				if variants[i][4] in bases:
					alt_count=bases.count(variants[i][4])
					#if float(alt_count)/len(bases)>0.1:
					if alt_count>1:
						variants[i][7]['SAO']='1'
			else:
				variants[i][7]['SAO']='1'

	print '\tSomatic SNVs found after filtering:\t%s' %len(filter(lambda x: x[7]['SAO']=='2' and x[7]['VC']=='SNV', variants))
	
	OutputFileName='%sgrmrm.vcf' %VcfFileName.split('vcf')[0]
	save_vcf_file(OutputFileName, par, summary_stat, variants)				
'''
					
def germline(pairs, ChrList):
	print 'Germline'
	Counter=0
	for line in open(pairs):
		pair=line.strip().split()
		for Chr in ChrList:
			#Command='qsub -N SMM_grm SMM_ht.bash %s %s' %(pair[1].split('.bam')[0], Chr)
                        Command='sbatch SMM_ht.bash %s %s' %(pair[1].split('.bam')[0], Chr)
			print Command
                        subprocess.call(Command, shell=True)
			#Counter+=1
			#print Counter
			#time.sleep(0.3)
			#if Counter==1500:
			#	time.sleep(3600)
			#	Counter=0
			#print Command

def germline_noDD(pairs, ChrList):
	print 'Germline'
	for line in open(pairs):
		pair=line.strip().split()
		for Chr in ChrList:
			#Command='qsub -N SMM_grm SMM_ht_noDD.bash %s %s' %(pair[1].split('.bam')[0], Chr)
                        Command='sbatch SMM_ht_noDD.bash %s %s' %(pair[1].split('.bam')[0], Chr)
			print Command
                        subprocess.call(Command, shell=True)
		    
	
def mut_type(Mut):
	if Mut=='CG' or Mut=='GC':
		return 0
	elif Mut=='CA' or Mut=='GT':
		return 1
	elif Mut=='AC' or Mut=='TG':
		return 2
	elif Mut=='AT' or Mut=='TA':
		return 3
	elif Mut=='AG' or Mut=='TC':
		return 4
	elif Mut=='CT' or Mut=='GA':
		return 5
	else:
		#print "Insert", Mut
		return 6



def get_stat(SampleName):
	par, summary_stat, variants=read_vcf_file(SampleName)		#load identified variants
	Spectra=['CG', 'CA', 'AC', 'AT', 'AG', 'CT', 'UKN']
	counter_s=[0]*7
	counter_g=[0]*7
	
	
	
	snv_som=filter(lambda x: x[7]['SAO']=='2' and x[7]['VC']=='SNV', variants)
	
	for var in snv_som:
		Mut=(var[3]+var[4]).upper()
		counter_s[mut_type(Mut)]+=1
	
	snv_grm=filter(lambda x: x[7]['SAO']=='1' and x[7]['VC']=='SNV', variants)
		
	for var in snv_grm:
		Mut=(var[3]+var[4]).upper()
		counter_g[mut_type(Mut)]+=1	
	
	#print 'Parameters'
	#print '%s' %(par[2:].split(':')[1].strip()) 



	print 'Sample,SF Bases,QSF Bases,sSNVs,gSNVs, sSNVs frq, gSNVs frq,CG,CA,AC,AT,AG,CT,UKN,Parameters' 
	print '%s,%s,%s,%s,%s,%s,%s,%s,%s' %(SampleName,summary_stat['TB'], summary_stat['QB'], len(snv_som), len(snv_grm), len(snv_som)*1000000.0/summary_stat['QB'], len(snv_grm)*1000000.0/summary_stat['QB'], ','.join(str(x) for x in counter_s),par[2:].split(':')[1].strip())
	
	

'''

def get_stat(SampleName):
	par, summary_stat, variants=read_vcf_file(SampleName)		#load identified variants
	Spectra=['CG', 'CA', 'AC', 'AT', 'AG', 'CT', 'UKN']
	counter_s=[0]*7
	counter_g=[0]*7
	
	
	
	snv_som=filter(lambda x: x[7]['SAO']=='2' and x[7]['VC']=='SNV', variants)
	
	for var in snv_som:
		Mut=(var[3]+var[4]).upper()
		counter_s[mut_type(Mut)]+=1
	
	snv_grm=filter(lambda x: x[7]['SAO']=='1' and x[7]['VC']=='SNV', variants)
		
	for var in snv_grm:
		Mut=(var[3]+var[4]).upper()
		counter_g[mut_type(Mut)]+=1	
	
	print 'Sample,Parameters'
	print '%s,%s' %(SampleName, par[2:].split(':')[1].strip()) 
	# print par[2:].strip()
	print 'SF Bases,QSF Bases,sSNVs,gSNVs, sSNVs frq, gSNVs frq' 
	print '%s,%s,%s,%s,%s,%s' %(summary_stat['SBS'], summary_stat['QSF'], len(snv_som), len(snv_grm), len(snv_som)*1000000.0/summary_stat['QSF'], len(snv_grm)*1000000.0/summary_stat['QSF'])
	
	print	
	print ',%s' %','.join(Spectra)
	print 'Somatic SNVs,%s' %','.join(str(x) for x in counter_s)
	print 'Germline SNVs,%s' %','.join(str(x) for x in counter_g)
	print



'''

def variants(Pairs, ChrList):
	#Get Reference and dbsnp
	for line in open('SMM_env.bash'):
		fields=line.strip().split('=')
		if fields[0]=='refgenome':
			refgenome=fields[1]	
		elif fields[0]=='dbsnp_gz':
			dbsnp=fields[1]
		else:
			pass	
	#start variant caller
	for pair in open(Pairs):
		s_file, g_file=pair.strip().split()
		for chr in ChrList:
			if chr=='-': 
				g_vcf='./ht/%s.vcf.gz' %(g_file.split('.bam')[0])
                        	g_rmdup='./ht/%s.%s.rmdup.bam' %(g_file.split('.bam')[0], 'all')
			else: 
				g_vcf='./ht/%s.%s.vcf.gz' %(g_file.split('.bam')[0], chr)	
				g_rmdup='./ht/%s.%s.rmdup.bam' %(g_file.split('.bam')[0], chr)	
			


			Command='sbatch SMM_vc.bash %s %s %s %s %s %s' %(s_file, chr, g_rmdup, refgenome, dbsnp, g_vcf)
                        #Command='qsub SMM_vc.bash %s %s %s %s %s %s' %(s_file, chr, g_rmdup, refgenome, dbsnp, g_vcf)
			#Command='qsub -cwd -b y -j y -l h_vmem=12G -N SMM_vc python main.py vc %s %s %s %s %s %s' %(s_file, chr, g_rmdup, refgenome, dbsnp, g_vcf)
                        print Command[:50]
                        subprocess.call(Command, shell=True)
			
	
def get_chr_list(ChrFile):
	ChrList=[]
	for rec in open(ChrFile):
		ChrList.append(rec.split()[0].strip())
	return ChrList

		
def main():
	#ChrList=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22', 'X', 'Y']
	#ChrList=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
	command=sys.argv[1]
	#SampleName=sys.argv[2]

	for line in open('SMM_env.bash'):
                fields=line.strip().split('=')
                if fields[0]=='refgenome':
                        refgenome=fields[1]
                elif fields[0]=='dbsnp_gz':
                        dbsnp=fields[1]
		elif fields[0]=='chr_file':
			ChrList=get_chr_list(fields[1])
                else:
                        pass



	if command=='vc':
		pairs=sys.argv[2]		
		variants(pairs, ChrList)
		#print 'Command is not defined yet'
	
	elif command=='germline':
		pairs=sys.argv[2]
		germline(pairs, ChrList)	

	elif command=='germline_noDD':
		pairs=sys.argv[2]
		germline_noDD(pairs, ChrList)	

	elif command=='merge':
		pairs=sys.argv[2]	
		for line in open(pairs):
			SampleName=line.split()[0]
			print SampleName
			combine(SampleName, ChrList)
	
	elif command=='stat':
		pairs=sys.argv[2]
		for line in open(pairs):
			SampleName=line.split()[0]
			VcfFileName='%svcf' %SampleName.split('bam')[0]
			get_stat(VcfFileName)
		
		
	else:
		print 'Unrecognized command: ', command
		print 'Commands: vc, merge, germline, germline_noDD, stat'


	

#--------------------------------
if __name__ == "__main__":
    main()
