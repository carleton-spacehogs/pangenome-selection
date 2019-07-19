# get all the necessary variability profiles

import os
import sys
import re
import time
import os.path
import subprocess
import glob
from subprocess import call

fs = sys.argv[1]
inAxial = open('axial-bins.txt','r+')
inMCR = open('mcr-bins.txt', 'r+')

lines = inAxial.readlines()
lines2 = inMCR.readlines()

for line in lines:
	lit = line.split('\t')
	sample,mag = lit[0],lit[1]
	mag = mag.rstrip('\n')
	if sample == fs:
		os.system('anvi-gen-variability-profile -p %s_scv_profile.db -c %s_contigs.db -C %s -b %s -o SCV_Profile/%s_%s_SCV.txt --engine CDN'%(sample,sample,sample,mag,sample,mag))
inAxial.close()

for line in lines2:
	lit = line.split('\t')
	sample,mag = lit[0],lit[1]
	mag = mag.rstrip('\n')
	if sample == fs:
		os.system('anvi-gen-variability-profile -p %s_scv_profile.db -c %s_contigs.db -C %s -b %s -o SCV_Profile/%s_%s_SCV.txt --engine CDN'%(sample,sample,sample,mag,sample,mag))
inMCR.close()