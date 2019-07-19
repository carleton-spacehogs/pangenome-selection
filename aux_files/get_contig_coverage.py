#Codes to get the self-to-self contig (or split--anvi'o terms for contig splitted into pieces) coverage for all the Sulfurovum ORFs in the pangenome.

import os
import sys
import re
import time
import os.path
import subprocess
import glob
import numpy
from collections import OrderedDict
from subprocess import call

path = r'./'
run = str(sys.argv[1])
#first, self-to-self contig coverage table.
self_to_self_coverage = {}

#all samples containing Sulfurovum MAGs
all_samples = ['FS844','FS854','FS856','FS866','FS881','FS896','FS898','FS900','FS903','FS904','FS906','FS908','FS917']
#get this for each sample

for sample in all_samples:
	contigs_path = sample+'_contigs.db'
	
	#run the anvi'o provess when necessary
	if(run=="run"):
		os.system('anvi-export-splits-and-coverages -p %s_scv_profile.db -c %s_contigs.db -O %s --report-contigs --use-Q2Q3-coverages'%(sample,sample,sample))
	
	#anvi'o uses uppercase
	out_sample = sample.upper()
	print(out_sample)
	#get the coverage for this sample first
	self_to_self_coverage[sample] = {}
	with open('%s-COVs.txt'%(out_sample),'r+') as f:
		lines = f.readlines()
		
		#first line is important
		title_split = lines[0].split('\t')
		index=0
		for title in title_split:
			if(title.split('_')[0] == out_sample):
				#find the right index to get only self-to-self coverage
				index = title_split.index(title)
				break
		#the main job		
		for line in lines[1:]:
			line_split = line.split('\t')
			self_to_self_coverage[sample][line_split[0]] = str(line_split[index])

#now look up all ORFs in the pangenome
with open('Sulfurovum_gene_clusters_summary_cleaned.txt','r+') as cluster_file:
	cluster_file_lines = cluster_file.readlines()
	sulfurovum_id = {}
	gene_to_contig_id = {}
	
	#sulfurovum id created and then at the same time also map gene to contig (or contig to gene)
	with open('sulfurovum-genomes.txt','r+') as genome_file:
		genome_lines = genome_file.readlines()
		for line in genome_lines[1:]:
			separate = line.split('\t')
			(genome,bin_id,collection_id) = (separate[0],separate[1],separate[2])
			sample_id = collection_id.split('_bins')[0]
			sulfurovum_id[genome] = sample_id+'-'+bin_id
			gene_to_contig_id[genome] = {}
			
			#the gene caller file helps us relate contig to gene and vice-versa
			allgenes_file = str('SUMMARY_'+sample_id+'/bin_by_bin/'+bin_id+'/'+bin_id+'-gene_calls.txt')
			with open(allgenes_file,'r+') as gene_call:
				for call_line in gene_call.readlines()[1:]:
					call_split = call_line.split('\t')
					(gene,contig) = (call_split[0],call_split[1])
					gene_to_contig_id[genome][gene]=contig
	
	coverage = {} #the all-important dictionary	
	
	#operate this for all ORFs
	for cluster_line in cluster_file_lines[1:]:
		cluster_split = cluster_line.split('\t')
		(unique_id,genome,caller_id) = (cluster_split[0],cluster_split[3],cluster_split[4])
		contig = gene_to_contig_id[genome][caller_id] #determine which contig it belongs to
		id_split = sulfurovum_id[genome].split('-')
		sample_id = id_split[0] #get the sample_id
		coverage[unique_id] = str(self_to_self_coverage[sample_id][contig]).rstrip('\n') #done
	
with open('All_Sulfurovum_self_to_self_coverage.txt','w') as all_cov:
	all_cov.write('unique_id	coverage'+'\n')
	for unique in coverage:
		all_cov.write(str(unique)+'\t'+str(coverage[unique])+'\n')