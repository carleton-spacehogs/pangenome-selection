#Codes to get a similarity matrix between the genomes in the pangenome
from __future__ import division

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

inFile = open('Sulfurovum_gene_clusters_summary_cleaned.txt','r+')
#fileLines = inFile.readlines
countCluster = {}
#counting gene cluster dictionary for each cluster

for line in inFile.readlines()[1:]:
	lineSplit = line.split('\t')
	(cluster_id,genome_id) = lineSplit[1],lineSplit[3]
	genome_number = int(genome_id.split('_')[1]) - 1 #indexing based on genome number
	if cluster_id not in countCluster:
		countCluster[cluster_id] = [0]*22
		countCluster[cluster_id][genome_number] += 1
	else:
		countCluster[cluster_id][genome_number] += 1

countGenome = [[0 for x in range(22)] for y in range(22)]
for cluster_id in countCluster:
	for i in range(22):
		if countCluster[cluster_id][i] > 0:
			for j in range(22):
				if countCluster[cluster_id][j] > 0:
					countGenome[i][j] += 1

proportionGenome = [[0 for x in range(22)] for y in range(22)]
for i in range(22):
	limit = countGenome[i][i]
	for j in range(22):
		proportionGenome[i][j] = countGenome[i][j]/limit
		
def write_matrix_to_textfile(a_matrix, file_to_write):
    def compile_row_string(a_row):
        return str(a_row).strip(']').strip('[').replace(', ','\t')

    with open(file_to_write, 'w') as f:
        for row in a_matrix:
            f.write(compile_row_string(row)+'\n')

    return True

write_matrix_to_textfile(proportionGenome, 'similarity_matrix.txt')
