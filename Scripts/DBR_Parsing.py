#!/usr/bin/python

#################################################################################
##																			   ##
##  Parse raw Illumina Read 1 FASTQ files with degenerate base regions (dbrs)  ##
##  to dictionaries for Maximum Likelihood Estimation of PCR duplicates		   ##
##																			   ##
#################################################################################

from collections import defaultdict
from subprocess import call, Popen, PIPE
import subprocess
import os, os.path
import sys
import json
import re
import itertools

# pre-process read 2 file -- currently commented out because it only needs to be done once (eventually will parse cmd args for file paths with option to re-construct this file or not)
# extract the headers
#subprocess.call("zcat /home/antolinlab/Downloads/Sample1_ACTTGA_L008_R2_001.fastq.gz | sed -n '1~4p' > /home/antolinlab/Downloads/headers_parse.txt", shell=True)
# extract the dbr sequences
#subprocess.call("zcat /home/antolinlab/Downloads/Sample1_ACTTGA_L008_R2_001.fastq.gz | sed -n '2~4p' | cut -c 1-8 > /home/antolinlab/Downloads/dbr_seqs_parse.txt", shell=True)
# interleave files so that dbr sequence is first, header is second (makes dictionary construction with DBRs as keys easier because the DBR line is the first encountered
#subprocess.call("paste -d'\n' /home/antolinlab/Downloads/dbr_seqs_parse.txt /home/antolinlab/Downloads/headers_parse.txt > header_seq_only.txt", shell=True)

        
# set default dictionary value type to list
dbr = defaultdict(list)

# set initial value of toggle to keep track of line number (even or odd)
isEven = True

# load file (read 2 contains DBR) and build dictionary
# keys will be DBRs -- the first 8 items in the line

with open('/home/antolinlab/Downloads/header_seq_only.txt', 'r') as f:
	for line in f:
		if isEven: # even lines = dbr sequences
			dbr.fromkeys([line])
			last_line=line
			#if the sequence file has the full seq (not just the dbr), do this:
			#dbr.fromkeys(line[:8])
		else: # odd lines = IDs
			ID = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', line)[1]
			dbr[last_line].append(ID)
		isEven = not isEven

x = itertools.islice(dbr.iteritems(), 0, 4)

for key, value in x:
	print key, value
		
#json.dump(dbr, '/home/antolinlab/Downloads/dbr_dict')

#myDict = json.load(myFilehandle)
