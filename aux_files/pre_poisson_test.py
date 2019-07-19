import os
import sys
import re
import time
import os.path
import subprocess
import glob
from subprocess import call
from scipy.stats import poisson

path = r'./'
import numpy

#do this for all SNV files of interest
for filename in glob.glob(os.path.join(path, 'SNV_Profile_contig_inc/*SNV.txt')):
	#process files
	sFile = open(filename,'r+')
	filename_delete = str(filename).split('/')
	filename = filename_delete[2]
	filename_delete = str(filename).split('_SNV')
	filename = filename_delete[0]
	
	#the resulting filename is in the form of sample_bin_numb
	filename_split = filename.split('_')
	(mega,bins) = (filename_split[0],str(filename_split[1]+'_'+filename_split[2]))
	
	#now go into the specific gene calls file to find the connection between contigs and gene calls
	allgenes_file = str('SUMMARY_'+mega+'/bin_by_bin/'+bins+'/'+bins+'-gene_calls.txt')
	print(allgenes_file)
	allFile = open(allgenes_file,'r+')
	
	#process contigs file
	contigs_path = str('SUMMARY_' + mega +'/bin_by_bin/' + bins + '/' + bins + '-contigs.fa')
	contigs_file = open(contigs_path,'r+')
	contigs_lines = contigs_file.readlines()
	
	#dictionary for contig length
	contigs = {}
	for i in range(len(contigs_lines)):
		#only count the title (id) line
		if (contigs_lines[i][0] == '>'):
			#the id of the contig is there
			con_id = contigs_lines[i][1:].rstrip('\n')
			con_now = con_id.split('_')
			con_id = str(con_now[0]+'_'+con_now[1])
			
			#initialize the length
			contigs[con_id]= 0
			j = i+1
			while (contigs_lines[j][0] != '>'):
				#the length of the line
				length = len(contigs_lines[j])
				con_length = contigs[con_id]
				
				#sum it up to the length of the contig
				contigs[con_id] = length + con_length
				j=j+1
				if(j>len(contigs_lines)-1):
					break
	
	#now we are going to find the length of the genome by just summing up the lengths of all the contigs in that genome
	genome_length=0
	for contig in contigs:
		genome_length = genome_length + contigs[contig]
	contigs_file.close()
	
	#now process the real file
	sFile_lines = sFile.readlines()
	allFile_lines = allFile.readlines()
	
	#number of SNVs in the genome
	genome_SNV = 0
	
	#gene-wise dictionaries
	numb_SNV_genes = {} #number of SNVs in the genes
	genes_outside_length = {} #the length of contig outside the gene itself
	genes_outside_SNV = {} #number of SNVs in the contig but outside the gene
	genes_in_contigs = {} #the contig id where the gene is
	coverage = {} #the coverage of the gene
	
	#contigs-wise dictionaries
	numb_SNV = {} #number of SNV in the contig
	
	for line in sFile_lines[1:]:
		lines_split = line.split('\t')
		(gene_length,sample,pos,gene,contig_id,cover)=(int(lines_split[2]),lines_split[3],lines_split[5],lines_split[6],lines_split[28],lines_split[11])
		
		#find the sample id
		sample_id_split = sample.split('_')
		sample_id = sample_id_split[0]
		
		#make sure that we are only looking at self-to-self mapping
		if (sample_id == mega):
			
			#we know that the # of SNV increases in the genome
			genome_SNV = genome_SNV+1
			
			#new numb_SNV if contig id is not there, and in that case also initialize new numb_SNV_genes
			if (contig_id not in numb_SNV):
				numb_SNV[contig_id]= 1
				numb_SNV_genes[gene] = 1
				
			else:
				old_numb = numb_SNV[contig_id]
				#add 1 SNV to the contig
				numb_SNV[contig_id] = old_numb + 1
				
				#add 1 SNV to the gene
				if (gene not in numb_SNV_genes):
					numb_SNV_genes[gene] = 1
				else:
					old = numb_SNV_genes[gene]
					numb_SNV_genes[gene] = old+1
	
	for line in allFile_lines[1:]:
		lines_split = line.split('\t')
		(gene,contig,gene_length) = (lines_split[0],lines_split[1],int(lines_split[3])-int(lines_split[2])+1)
		
		#compute the length of contig outside the gene
		genes_outside_length[gene] = contigs[contig]-gene_length
		
		#figure out in which contig the gene is
		genes_in_contigs[gene] = contig
		
		if gene in numb_SNV_genes:
			#compute numb of SNVs outside the gene inside the contig
			genes_outside_SNV[gene] = numb_SNV[contig] - numb_SNV_genes[gene]
			
		elif contig in numb_SNV:
			#compute the same thing if there is no SNV inside the gene
			genes_outside_SNV[gene] = numb_SNV[contig]
			numb_SNV_genes[gene] = 0
			
		else:
			#otherwise, there is no SNV in the contig itself
			genes_outside_SNV[gene] = 0
			numb_SNV_genes[gene] = 0

	
	#also gene wise dictionary for p_val (or chance)
	p_val = {}
	exp_val = {}
	#calculate this for all genes with some region in the contig outside the gene itself
	for gene in genes_outside_length:
		if (gene != '-1'):
			contig_length = genes_outside_length[gene]
			expected_val = float(genome_SNV) * float(contig_length) / float(genome_length)
			exp_val[gene] = expected_val
            
			#possible_val = float(genome_length)-float(contig_length)
			#number = genes_outside_SNV[gene] used this in analyses before
			contig = genes_in_contigs[gene]
			if (contig not in numb_SNV):
				numb_SNV[contig] = 0
			number = numb_SNV[contig]
			print(number)
			#chance, not probability
			p_val[gene] = str(float(poisson.cdf(number, expected_val)))
		else:
			p_val[gene] = 'NaN'
			exp_val[gene] = 'NaN'

	#write everything up
	outFile = open(str(filename + '_Poisson_07_10.txt'), 'w')
	#outFile.write('gene	p-val	SNV_inside	SNV_outside	length_outside	expected' + '\n')
	outFile.write('gene	p-val	total_SNV	expected' + '\n')
	for gene in p_val:
		outFile.write(gene + '\t' + p_val[gene] + '\t' + str(numb_SNV[genes_in_contigs[gene]])+'\t'+str(exp_val[gene])+'\n')
	outFile.close()
	sFile.close()