#!/usr/bin/python

# Note: hard-coded file paths to be replaced with cmd args

from subprocess import call, Popen, PIPE
import subprocess
import os, os.path
import sys
import numpy as np
import pandas as pd

'''
# Quality filter with FASTX Toolkit
subprocess.call("zcat ~/Downloads/Sample1_ACTTGA_L008_R1_001.fastq.gz | fastq_quality_filter -q 25 -p 50 -z -o ~/Downloads/Sample1_ACTTGA_L008_R1_001_filtered.fastq.gz", shell=True)
# to-do: use zcat | wc -l to determine how many lines were removed after filtering (print to screen)

# Demultiplex with FASTX Toolkit
subprocess.call("zcat /home/antolinlab/Downloads/Sample1_ACTTGA_L008_R1_001_filtered.fastq.gz | /usr/local/bin/fastx_barcode_splitter.pl --bcfile /home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/pilot_barcode_file --prefix /home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed/pilot_demultiplex_ --bol", shell=True)
# to-do: automatically generate a new directory "Demultiplexed" to house the new data (and ONLY the new data)
'''

# Generate "-s" arguments to Stacks denovo_map.pl
mp_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed/'
s_list=[]
for i in os.listdir(mp_dir):
    # the proper file name
    proper_file = i+'.fq'
    # full path to file for renaming
    full_path=os.path.join(mp_dir, i)
    new_path=os.path.join(mp_dir, proper_file)
    # rename the file as .fq so stacks can find it
    os.rename(full_path, new_path)
    # generate a string for the -s stacks argument
    s_list.append('-s ')
    s_list.append(new_path)
    s_list.append(' ')
# the last file in the dir is "unmatched"... remove it and the final trailing space
s_list = s_list[0:len(s_list)-4]
# join list into a string to pass to denovo_map.pl
formatted_list = ''.join(s_list)

# Run denovo_map.pl
build_args = ['/home/antolinlab/Downloads/stacks-1.31/scripts/denovo_map.pl ', 
              '-e ', '/home/antolinlab/Downloads/stacks-1.31/scripts '
              '-o ', '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/ ',
              '-m ', '10 ',
              '-n ', '2 ',
              '-t ',
              '-b ', '1 ',
              '-D ', 'initial_assembly ',
              '-S ',
              formatted_list]
denovo_call = ''.join(build_args)

#subprocess.call(denovo_call, shell=True)

# Extract consensus sequence from tags.tsv Stacks output, save as FASTQ

ref_path = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/batch_1.catalog.tags.tsv'
col_names = ['SQLID', 'SampleID', 'LocusID', 'Chrom', 
             'BP', 'Strand', 'SeqType', 'Stack', 'SeqID', 
             'Seq', 'Deleveraged', 'Blacklisted', 'Lumberjack', 'LogLik']
ref = pd.read_csv(ref_path, comment='#', sep='\t', names=col_names)
#print ref
locus = ref.ix[:,'LocusID']
new_header = pd.concat([pd.DataFrame(['> pseudoref_']*len(locus)), pd.DataFrame(locus)], axis=1)
reformatted = [None]*len(locus)
for i in range(len(locus)-1):
    #print new_header.iloc[i:,0]   
    reformatted[i]=new_header.iloc[i:,0]+str(new_header.iloc[i:,1])
print reformatted
#seq = ref.ix[:,'Seq']
#print seq

# Call BWA-MEM for reference based mapping
