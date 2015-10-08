#!/usr/bin/python

import os
import re
import subprocess

# run QC on the library-demultiplexed ddRADseq data with DBRs:

#m.group(1) for l in lines for m in [regex.search(l)] if m

# check for barcodes, read 1
def barcode_r1(directory, out_name):
	files = os.listdir(directory)
	print files
	#r1_files = []
	for f in files:
		r1 = re.findall(r'R1', f)
		if r1:
			#r1_files.append(f)
			barcode_check_call = 'zcat ' + directory + f + " | sed -n '2~4p' | cut -c 1-5 | sort | uniq -c | sort -nr -k 1 > " + directory + f + '_' + out_name
			subprocess.call(barcode_check_call, shell=True)

# check for cutsite, read 1
#cutsite_check_call = 'zcat ' + library_name + " | sed -n '2~4p' | cut -c 6-10 | sort | uniq -c | sort -nr -k 1 > " + out_name
#subprocess.call(cutsite_call, shell=True)

# check for degeneracy, read 2
#degeneracy_check_call = 'zcat ' + library_name + " | sed -n 2~4p | cut -c 1-8 | sort | uniq -c | sort -nr -k 1 > " + out_name
#subprocess.call(degeneracy_call, shell=True)

barcode_r1('/home/antolinlab/Downloads/CWD_RADseq/', 'barcode_check')