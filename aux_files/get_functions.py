#Get the count table for COG categories across Number of MAGs

import os
import sys
import re
import time
import os.path
import subprocess
import glob
from subprocess import call

Genome = 'Sulfurovum_gene_clusters_summary_cleaned.txt' #the cluster file
max_genome = 22 #number
taxon = Genome.split('_')[0]
inFile = open(Genome,'r+')
lines = inFile.readlines()

max_function = {} #dictionary for count

for line in lines[1:]:
	row=line.split('\t')
	unique_id,cluster_id,group,genome,gene_callers_id,num_genomes,cog,function=row[0],row[1],row[2],row[3],row[4],row[5],row[9],row[10]
	functions = function.split('|') #some gene clusters have more than just one category
	for function_ind in functions:
		if function_ind not in max_function:
			max_function[function_ind] = [0]*int(max_genome) #count initializing
		max_function[function_ind][int(num_genomes)-1] += 1 #count + 1

inFile.close()
				
max_func = open('%s_function_list.txt'%(taxon),'w')
grouping = []
for i in range(int(max_genome)):
	grouping.append(str(i+1))
max_func.write('Category' + '\t' + '\t'.join(grouping) + '\n')

for func_max in max_function:
	for i in range(int(max_genome)):
		max_function[func_max][i] = str(max_function[func_max][i])
	max_func.write(func_max + '\t' +'\t'.join(max_function[func_max]) + '\n')
max_func.close()
