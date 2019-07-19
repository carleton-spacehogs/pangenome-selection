#this script collects all the pN/pS ratio files for all the 

import os
import sys
import re
import time
import os.path
import subprocess
import glob
from subprocess import call

#get the genome name info first
genomeFile = open('sulfurovum-genomes.txt','r+')
genomeLines = genomeFile.readlines()
genomeDB = {}

#create a mini dictionary
for line in genomeLines[1:]:
	genomeSplit = line.split('\t')
	genome,bin_id,collection_id = genomeSplit[0],genomeSplit[1],genomeSplit[2]
	bin_name = collection_id.split('_')[0] + '_' + bin_id
	genomeDB[genome] = bin_name
	
genomeFile.close()


clusterFile = open('Sulfurovum_gene_clusters_summary.txt','r+')
clusterLines = clusterFile.readlines()
outFile = open('Sulfurovum_pn_ps_summary.txt','w')
outFile.write('cluster_ID	pN_pS	genome_name	gene_caller	num_genomes	category	cog_function' + '\n')

# start summarizing
for line in clusterLines[1:]:
	pn_ps=''
	clusterSplit = line.split('\t')
	cluster,genome_name,gene_caller,num_genomes,category,cog_function = clusterSplit[1],clusterSplit[3],clusterSplit[4],clusterSplit[5],clusterSplit[9],clusterSplit[11]
	genomeInfo = genomeDB[genome_name]
	# all files are in a specific folder
	with open('pn_ps_03_25/%s/pN_pS_ratio.txt'%(genomeInfo),'r+') as readingFile:
		readingLines = readingFile.readlines()
		firstSplit = readingLines[0].split('\t')
		for i in range(1,len(firstSplit)+1):
			if firstSplit[i].split('_')[0] == genomeInfo.split('_')[0]: #self-to-self only
				num = i
				break
		for line in readingLines[1:]:
			lineSplit = line.split('\t')
			if lineSplit[0] == gene_caller:
				pn_ps = lineSplit[num]
				break
	summary = [str(cluster),str(pn_ps),genome_name,str(gene_caller),str(num_genomes),category,cog_function]
	outFile.write('\t'.join(summary)+ '\n')

clusterFile.close()
outFile.close()
				