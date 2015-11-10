#!/usr/bin/python

###################################################
### Influence of quality filters on DBR distr.  ###
###################################################

# Note: hard-coded file paths to be replaced with cmd args

from subprocess import call, Popen, PIPE
import subprocess
import os, os.path
import sys
import numpy as np


# Vary minimum median quality score
def dbr_qual_sens(infile, outdir, outprefix, qual):
    
    if not os.path.exists(outdir):
            os.makedirs(outdir)
    
    while qual < 40:
    
        print ['Filtering at median quality = ' + str(qual)]
    
        # Quality filter with FASTX Toolkit
        out_file = outdir + '/' + outprefix + str(qual) + ".fastq.gz"
        filter_call = "zcat " + infile + " | fastq_quality_filter -q " + str(qual) + " -p 50 -z -o " + out_file
        subprocess.call(filter_call, shell=True)
    
        # Generate DBR distribution
        count_out = str(qual) + "median_quality_degeneracy_count.txt"
        count_call = "zcat " + out_file + " | sed -n '2~4p' | cut -c 2-9 | sort | uniq -c | sort -nr -k 1 > " + outdir + '/' + count_out
        subprocess.call(count_call, shell=True)
    
        print filter_call
        print count_call
    
        # Recursively call function so that the filtering happens on iteratively smaller input datasets (hopefully to increase speed and avoid recalculation median quality multiple times for low quality sequences
        return dbr_qual_sens(out_file, outdir, outprefix, qual+5)

#out_dir = "/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/QualDBR_Sensitivity"
#original_data = "~/Downloads/Sample1_ACTTGA_L008_R1_001.fastq.gz"
#out_prefix = "Sample1_ACTTGA_L008_R1_001_filtered_temp"

out_dir = "/home/antolinlab/Downloads/CWD_RADseq/Library12_DBR_distr/"
original_data = "/home/antolinlab/Downloads/CWD_RADseq/Library12_S65_L008_R2_001.fastq.gz"
out_prefix = "Library12_S65_L008_R2_001_dbr_distr"

dbr_qual_sens(original_data, out_dir, out_prefix, 35)


