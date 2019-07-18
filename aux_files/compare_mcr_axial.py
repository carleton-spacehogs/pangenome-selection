# summarize the region count file

import os
import sys
import re
import time
import subprocess
import glob
import numpy
from subprocess import call
from scipy.stats import binom

genus = 'Sulfurovum' #the genus, which can be substituted
inFile = str(genus+'_vent_count_per_cluster.txt')
cluster_file = str(genus+'_gene_clusters_summary_cleaned.txt')
axial_lim = int(sys.argv[1]) #min number of genome in axial; default to be 12
mcr_lim = int(sys.argv[2]) #min number of genome in MCR; default to be 8

genus=cluster_file.split('_')[0]

axial_enriched={}
axial_enriched_list = []
mcr_enriched_list = []
cluster_dict={}
category_dict={}
cog_dict={}

cluster = open(cluster_file,'r+')

inFile = open(inFile,'r+') 
for line in inFile.readlines()[1:]:
	parsed = line.split('\t')
	(cluster_id,axial_count,mcr_count) = (parsed[0],float(parsed[1]),float(parsed[2]))
	if (axial_count >= axial_lim-1 or mcr_count >= mcr_lim-1):
		axial_enriched[cluster_id]=binom.cdf(axial_count, axial_count+mcr_count, float(13)/float(22))

score={}
		
for cluster_id in axial_enriched:
	for seq in cluster_dict:
		if cluster_dict[seq][0] == cluster_id:
			initial = cluster_dict[seq]
			initial.append(axial_enriched[cluster_id])
			axial_enriched_list.append(initial)
			if initial[0] not in score:
				score[initial[0]] = axial_enriched[cluster_id]
			

for line in axial_enriched_list:
	axial_out.write('\t'.join(str(e) for e in line) + '\n')
	
axial_out.close()

out2 = open('%s_Axial_enriched_compressed.txt'%(genus),'w')
out2.write('cluster id	Enrichment	Category	COG' + '\n')
for line in score:
	cat = category_dict[line]
	cat_split = cat.split('|')
	cog = cog_dict[line]
	for category in cat_split:
		out2.write(line+'\t'+str(score[line])+'\t'+category+'\t'+cog+ '\n')
