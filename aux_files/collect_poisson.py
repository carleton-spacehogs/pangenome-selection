# this script collects all the 

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

genus = "Sulfurovum" # default genus
poisson_file_table={}
poisson_table = {}
SNV_table ={}
exp_table = {}

if (genus=="Sulfurovum"):
	cluster_file = open('Sulfurovum_gene_clusters_summary_cleaned.txt','r+')
	genome_file = open('sulfurovum-genomes.txt','r+')
	
else:
	cluster_file = open('Arcobacter_gene_clusters_summary_cleaned.txt','r+')
	genome_file = open('arcobacter-genomes.txt','r+')
	
genome_lines = genome_file.readlines()

# the following lines create the table to map MAGs to specific files and samples
for line in genome_lines[1:]:
	separate = line.split('\t')
	(genome,bin_id,collection_id) = (separate[0],separate[1],separate[2])
	collection_ids = collection_id.split('_bins') 
	collection_id = collection_ids[0]
	poisson_file = str(collection_id + '_' + bin_id)
	poisson_file_table[genome] = poisson_file

for filename in glob.glob(os.path.join(path, '*_Poisson_07_10.txt')):
	poisson_file = open(filename,'r+') # read the Poisson score file
	name_split = str(filename).split('_Poisson')
	name = name_split[0]
	names = name.split('/')
	name = names[1]
	poisson_table[name] = {}
	SNV_table[name] = {}
	exp_table[name] = {}
	poisson_read = poisson_file.readlines()
	for line in poisson_read[1:]: # put all the information into dictionaries
		linet = line.split('\t')
		(gene,poiss,snv,exp) = (linet[0],linet[1],linet[2],linet[3])
		poisson_table[name][gene] = str(poiss)
		SNV_table[name][gene] = str(snv)
		exp_table[name][gene] = str(exp)
	poisson_file.close()
	
cluster_reads = cluster_file.readlines()

# create global dictionaries which summarize the data on a gene-by-gene basis
gene_poiss = {}
num_genomes = {}
SNV = {}
expected = {}

for line in cluster_reads[1:]: #run for all genes
	linet = line.split('\t')
	(unique,cluster_id,genome,gene_call,num_genome,num_occurrence,cat,func) = (linet[0],linet[1],linet[3],linet[4],linet[5],linet[6],linet[9],linet[11])
	
	poisson_file = poisson_file_table[genome]
	if (gene_call in poisson_table[poisson_file]): # update global dictionaries
		gene_poiss[unique]=poisson_table[poisson_file][gene_call]
		num_genomes[unique]=num_genome
		SNV[unique]=SNV_table[poisson_file][gene_call]
		expected[unique]=exp_table[poisson_file][gene_call]
		
outFile = open('%s_poisson_sweep.txt'%(genus),'w')
outFile.write('unique_id	pv_val	numb	SNV	expected'+'\n')
for cluster in gene_poiss:
	outFile.write(cluster + '\t' + gene_poiss[cluster] + '\t'+num_genomes[cluster]+ '\t'+ SNV[cluster] + '\t'+ expected[cluster] +'\n')
outFile.close()