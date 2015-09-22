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

# pre-process read 2 file

# extract DBR sequences?
DBRs=False
# extract R1 sequences?
R1_seqs=True
# save the collated FASTQ/(DBR|R1_seq) file? (this file is needed to facilitate dictionary generation 
collate=True
# create the dictionary with the DBR as key?
generate_dict_DBR=False
# create the dictionary with the FASTQ header as key?
generate_dict_FASTQ=False
# create a dictionary of {FASTQ:R1} ?
generate_dict_R1_seq=True
# print output of dictionary to screen?
test_dict=True
# save dictionary?
save_DBR_dict=False
save_R1_dict=True

if DBRs: # extract the dbr sequences
    print 'Extracting FASTQ headers and DBR sequences from Read 2.'
    subprocess.call("zcat /home/antolinlab/Downloads/Sample1_ACTTGA_L008_R2_001.fastq.gz | sed -n '1~4p' > /home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/headers_parse.txt", shell=True)
    subprocess.call("zcat /home/antolinlab/Downloads/Sample1_ACTTGA_L008_R2_001.fastq.gz | sed -n '2~4p' | cut -c 1-8 > /home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/dbr_seqs_parse.txt", shell=True)

if R1_seqs:
    print 'Extracting FASTQ headers and full sequences from Read 1.'
    subprocess.call("zcat /home/antolinlab/Downloads/Sample1_ACTTGA_L008_R1_001.fastq.gz | sed -n '1~4p' > /home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_headers_parse.txt", shell=True)
    subprocess.call("zcat /home/antolinlab/Downloads/Sample1_ACTTGA_L008_R1_001.fastq.gz | sed -n '2~4p' > /home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_seqs_parse.txt", shell=True)

if generate_dict_DBR:
    
    keyFile = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/dbr_seqs_parse.txt'
    valueFile = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/headers_parse.txt'
    
    if collate: # interleave files so that dbr sequence is first, header is second (makes dictionary construction with DBRs as keys easier because the DBR line is the first encountered
        print 'Parsing Read 2 headers and DBRs for dictionary building.'
        collate_call = "paste -d '\\n' " + keyFile + ' ' + valueFile + " > /home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/header_seq_only.txt"
        #print collate_call
        subprocess.call(collate_call, shell=True)
    
    print 'Building dictionary {DBR:FASTQ}.'        
    # set default dictionary value type to list
    dbr = defaultdict(list)

    # set initial value of toggle to keep track of line number (even or odd)
    isEven = True
    
    # load file (read 2 contains DBR) and build dictionary
    # keys will be DBRs -- the first 8 items in the line

    with open('/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/header_seq_only.txt', 'r') as f:
        for line in f:
            if isEven:
                dbr.fromkeys([line])
                last_line=line
                #if the sequence file has the full seq (not just the dbr), do this:
                #dbr.fromkeys(line[:8])
            else: # odd lines = IDs
                ID = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', line)[1]
                dbr[last_line].append(ID)
            isEven = not isEven
            
if generate_dict_FASTQ:
    
    keyFile = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/headers_parse.txt'
    valueFile = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/dbr_seqs_parse.txt'
    
    if collate: 
        print 'Parsing Read 2 headers and DBRs for dictionary building.'
        collate_call = "paste -d '\\n' " + keyFile + ' ' + valueFile + " > /home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/header_seq_only.txt"
        #print collate_call
        subprocess.call(collate_call, shell=True)
    
    print 'Building dictionary {FASTQ:DBR}.'        
    # set default dictionary value type to list
    dbr = defaultdict(list)

    # set initial value of toggle to keep track of line number (even or odd)
    isEven = True
    
    with open('/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/header_seq_only.txt', 'r') as f:
        for line in f:
            if isEven:
                ID = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', line)[1]
                dbr.fromkeys([ID])
                last_line=ID
                #if the sequence file has the full seq (not just the dbr), do this:
                #dbr.fromkeys(line[:8])
            else: # odd lines = IDs
                dbr[last_line].append(line.rstrip('\n'))
            isEven = not isEven

if generate_dict_R1_seq:
    
    keyFile = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_headers_parse.txt'
    valueFile = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_seqs_parse.txt'
    
    if collate:
        print 'Parsing Read 1 headers and sequences for dictionary building.'
        collate_call = "paste -d '\\n' " + keyFile + ' ' + valueFile + " > /home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/header_R1_seq_only.txt"
        #print collate_call
        subprocess.call(collate_call, shell=True)
    
    print 'Building dictionary {FASTQ:R1_seq}.'        
    # set default dictionary value type to list
    R1 = defaultdict(list)

    # set initial value of toggle to keep track of line number (even or odd)
    isEven = True
    
    # load file (read 2 contains DBR) and build dictionary
    # keys will be DBRs -- the first 8 items in the line

    with open('/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/header_R1_seq_only.txt', 'r') as f:
        for line in f:
            if isEven:
                ID = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', line)[1]
                R1.fromkeys([ID])
                last_line=ID
            else: # odd lines = sequences
                R1[last_line].append(line.rstrip('\n'))
            isEven = not isEven

if test_dict: # check construction by printing first entries to screen
    if generate_dict_R1_seq:
        print 'Checking R1 dictionary format.'
        x = itertools.islice(R1.iteritems(), 0, 4)
        for key, value in x:
            print key, value
    else:
        print 'Checking DBR dictionary format.'
        x = itertools.islice(dbr.iteritems(), 0, 4)
        for key, value in x:
            print key, value
           
if save_DBR_dict:
    dbr_dict_out = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/dbr_dict.json'
    print 'Writing dictionary to ' + dbr_dict_out
    with open('/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/dbr_dict', 'w') as fp:	  	
        json.dump(dbr, fp)
        
if save_R1_dict:
    R1_dict_out = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_dict.json'
    print 'Writing dictionary to ' + R1_dict_out
    with open('/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_dict', 'w') as fp:          
        json.dump(R1, fp)

