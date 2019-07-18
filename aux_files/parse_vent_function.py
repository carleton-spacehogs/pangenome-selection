# count the number of gene clusters found in each vent region

import os
import sys
import re
import time
import subprocess
import glob
import numpy
from subprocess import call

cluster_file = 'sulfurovum_vent.txt' #vent ID file
vent_file = 'Sulfurovum_gene_clusters_summary_cleaned.txt' #original cluster file

genus=vent_file.split('_')[0]

axial_list=[]
mcr_list=[]
vent_count_dict={}
genome_count_dict={}
axial_genome_list_dic={}
mcr_genome_list_dic={}

vents=open(vent_file,'r+')

clusters=open(cluster_file,'r+')
for line in clusters.readlines():
	line=line.rstrip('\n')
	parsed_line = line.split('\t')
	print(parsed_line[2])
	if(parsed_line[2]=='Axial'):
		axial_list.append(parsed_line[0]) #append Axial to the list
	else:
		mcr_list.append(parsed_line[0])
clusters.close()

print(axial_list)

for line in vents.readlines()[1:]: #now look through the cluster file and do the counting
	parsed_line = line.split('\t')
	cluster=parsed_line[1]
	if(parsed_line[1] not in vent_count_dict): #new count
		if parsed_line[3] in axial_list:
			vent_count_dict[parsed_line[1]]=[1,0]
			axial_genome_list_dic[cluster]=[parsed_line[3]]
			mcr_genome_list_dic[cluster]=[]
		else:
			vent_count_dict[parsed_line[1]]=[0,1]
			mcr_genome_list_dic[cluster]=[parsed_line[3]]
			axial_genome_list_dic[cluster]=[]
		genome_count_dict[parsed_line[1]]=1
	else: #existent count
		if parsed_line[3] in axial_list:
			if parsed_line[3] not in axial_genome_list_dic[cluster]:
				vent_count_dict[parsed_line[1]]=[vent_count_dict[parsed_line[1]][0]+1,vent_count_dict[parsed_line[1]][1]]
				axial_genome_list_dic[cluster].append(parsed_line[3])
		else:
			if parsed_line[3] not in mcr_genome_list_dic[cluster]:
				vent_count_dict[parsed_line[1]]=[vent_count_dict[parsed_line[1]][0],vent_count_dict[parsed_line[1]][1]+1]
				mcr_genome_list_dic[cluster].append(parsed_line[3])
		genome_count_dict[parsed_line[1]]=genome_count_dict[parsed_line[1]]+1
		
outfile=open('%s_vent_count_per_cluster.txt'%(genus),'w')
outfile.write('cluster	axial count	mcr count	total occurrence' + '\n')
for cluster in vent_count_dict:
	outfile.write('\t'.join([cluster,str(vent_count_dict[cluster][0]),str(vent_count_dict[cluster][1]),str(genome_count_dict[cluster])])+'\n')
outfile.close()
		