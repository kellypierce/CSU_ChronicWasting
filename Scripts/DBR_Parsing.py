#!/usr/bin/python

#################################################################################
##																			   ##
##  Parse raw Illumina Read 1 FASTQ files with degenerate base regions (dbrs)  ##
##  to dictionaries for Maximum Likelihood Estimation of PCR duplicates		   ##
##																			   ##
#################################################################################

import json
import re
import itertools
import gzip

## New version 23 September 2015: DBR and R1 dictionaries made directly from FASTQ files
## No need for subprocess calls to SED for FASTQ file parsing
        
def R1_dict(fq_path, test_dict = False, save_path = None): 
    print 'Creating {ID: (full seq, +, qual)} dictionary.'
    R1 = {}
    fq_line = 1
    # gzip.open() will work with both compressed and uncompressed files
    with gzip.open(fq_path, 'r') as r1:
        for line in r1:
            if fq_line == 1:
                ID = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', line)[1]
                fq_line = 2
            elif fq_line == 2:
                seq = line
                fq_line = 3
            elif fq_line == 3:
                fq_line = 4
            elif fq_line == 4:
                qual = line
                R1[ID]=(seq, '+', qual)
                fq_line = 1
    if test_dict:
        print 'Checking Read 1 dictionary format.'
        x = itertools.islice(R1.iteritems(), 0, 4)
        for key, value in x:
            print key, value
    if save_path:
        print 'Writing dictionary to ' + save_path
        with open(save_path, 'w') as fp:          
            json.dump(R1, fp)

def DBR_dict(fq_path, dbr_start, dbr_stop, test_dict = False, save = None):
    if not checkFile(fq_path):
        raise IOError("where is the input file: %s" % fq_path)
    info('Creating {ID: dbr} dictionary from %s.' % fq_path)
    dbr = {}
    fq_line = 1
    if fq_path.endswith('gz'):
        openFxn = gzip.open
    else:
        openFxn = open
    with openFxn(fq_path, 'r') as db:
        for line in db:
            if fq_line == 1:
                ID = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', line)[1]
                fq_line = 2
            elif fq_line == 2:
                seq = list(line) # split the sequence line into a list
                tag = ''.join(seq[dbr_start:dbr_stop])
                dbr[ID] = tag
                fq_line = 3
            elif fq_line == 3:
                fq_line = 4
            elif fq_line == 4:
                fq_line = 1
    if test_dict:
        print 'Checking DBR dictionary format.'
        x = itertools.islice(dbr.iteritems(), 0, 4)
        for key, value in x:
            print key, value
        #print dbr['8:1101:15808:1492'] # this is the first entry in /home/antolinlab/Downloads/CWD_RADseq/pear_merged_Library12_L8.assembled.fastq
    if save:
    	fq_name = os.path.splitext(fq_path)[0]
    	save_path = save + '_DBRdict'
        print 'Writing dictionary to ' + save_path
        with open(save_path, 'w') as fp:          
            json.dump(dbr, fp)

# test functions:
#DBR_dict(fq_path = '/home/antolinlab/Downloads/Sample1_ACTTGA_L008_R2_001.fastq', 
#         dbr_start=0, dbr_stop=8, test_dict=True)
# if using Pear assembled data, need to read the string from the other end...
#DBR_dict(fq_path = '/home/antolinlab/Downloads/CWD_RADseq/pear_merged_Library12_L8.assembled.fastq', 
#         dbr_start=-9, dbr_stop=None, test_dict=True) #-9:None captures the last 8 values in the sequence string
#R1_dict('/home/antolinlab/Downloads/Sample1_ACTTGA_L008_R1_001_filtered.fastq', True)





