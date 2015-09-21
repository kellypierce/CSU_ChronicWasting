#!/usr/bin/python

###################################################
### De novo pipeline for processing READ 1 only ###
###################################################

# Note: hard-coded file paths to be replaced with cmd args

from subprocess import call, Popen, PIPE
import subprocess
import os, os.path
import sys
import numpy as np
import re

filter = False
demultiplex = False
trim = False
assemble = False
pseudoref = False
refmapping = True

mp_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed/'
prefix = 'pilot_demultiplex_' # eventually make cmd arg
suffix = '_trimmed.fq' # eventually make cmd arg
prefix_path = mp_dir + 'pilot_demultiplex_'
new_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed_Trimmed/'            
pseudoref_path = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/pseudoref.fa'
    
if filter:
    
    print 'Quality filtering.'
    
    # Quality filter with FASTX Toolkit
    subprocess.call("zcat ~/Downloads/Sample1_ACTTGA_L008_R1_001.fastq.gz | fastq_quality_filter -q 25 -p 50 -z -o ~/Downloads/Sample1_ACTTGA_L008_R1_001_filtered.fastq.gz", shell=True)
    # to-do: use zcat | wc -l to determine how many lines were removed after filtering (print to screen)

if demultiplex:
    
    print 'Demultiplexing sequence data.'
    
    # Demultiplex with FASTX Toolkit (note: does not remove barcodes)
    if not os.path.exists(mp_dir):
        os.makedirs(mp_dir)
    demultiplex_call = 'zcat /home/antolinlab/Downloads/Sample1_ACTTGA_L008_R1_001_filtered.fastq.gz | /usr/local/bin/fastx_barcode_splitter.pl --bcfile /home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/pilot_barcode_file --prefix ' + prefix_path + ' --bol'
    subprocess.call(demultiplex_call, shell=True)

if trim:
    
    print 'Trimming DBR and enzyme cut sites.'
    
    # new directory for trimmed files
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    for i in os.listdir(mp_dir):
    
        # the proper file name
        proper_file_trim = i+suffix
    
        # full path to file for trimming
        full_path=os.path.join(mp_dir, i)
        new_path=os.path.join(new_dir, proper_file_trim)
    
        # Remove barcodes and R1 enzyme cut site (first
        trim_call = "/usr/local/bin/fastx_trimmer -f 11 -i " + full_path + " -o " + new_path
        subprocess.call(trim_call, shell=True)

if assemble:
    
    print 'Assembling sequences de novo'
    
    s_list=[]
    
    # Generate "-s" arguments to Stacks denovo_map.pl; trim out barcode and cut site sequences
    for i in os.listdir(mp_dir):

        # the proper file name
        proper_file_trim = i+suffix
    
        # full path to trimmed file
        new_path=os.path.join(new_dir, proper_file_trim)
    
        # rename the file as .fq so stacks can find it -- not necessary if trimming first as the new name can be specified
        # os.rename(full_path, new_path)
    
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
                  '-S ', formatted_list]

    denovo_call = ''.join(build_args)
    subprocess.call(denovo_call, shell=True)

if pseudoref:
    
    print "Preparing pseudoreference genome."
    # Extract consensus sequence from tags.tsv Stacks output, save as FASTQ
    # This actually works! no pandas or numpy needed
    ref_path = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/batch_1.catalog.tags.tsv'
    ids = []
    seqs = []
    with open(ref_path, 'r') as y:
        next(y)
        for l in y:
            fields = l.split("\t")
            new_id = '>'+fields[2]+'_pseudoref'
            ids.append(new_id)
            seqs.append(fields[9])
    fastq_flat = [val for pair in zip(ids, seqs) for val in pair]
    fastq_interleaved = '\n'.join(fastq_flat)

    print "Writing pseudoreference genome to file " + pseudoref_path
    with open(pseudoref_path, 'a') as fq: # need to check if it already exists
        fq.write(fastq_interleaved)
    
    # Index the new pseudoreference genome
    print 'Indexing the reference genome.'
    index_call = '/home/antolinlab/Downloads/bwa.kit/bwa index ' + pseudoref_path
    print index_call
    subprocess.call(index_call, shell = True)

if refmapping:
    
    # Call BWA-MEM for reference based mapping
    
    assembly_path = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Assembled/'
    if not os.path.exists(assembly_path):
        os.makedirs(assembly_path)

    for i in os.listdir(new_dir): # for every demultiplexed & trimmed fastq...
        rex = re.compile(r'\d+') # the regex should be a cmd line argument probably, since the sample name is likely to vary....
        if rex.search(i):
            fname, fext = os.path.splitext(i)
            print 'Reference mapping ' + fname
            read_group_header = '"@RG\\tID:' + fname + '\\tPL:Illumina\\tLB:' + fname + '"' 
            bwa_mem_call = '/home/antolinlab/Downloads/bwa.kit/bwa mem -M -R ' + read_group_header + " " + pseudoref_path + ' ' + new_dir + i + ' > ' + assembly_path + fname + '.sam'
            print bwa_mem_call
            subprocess.call(bwa_mem_call, shell=True)
            

