#!/usr/bin/python

###################################################
### De novo pipeline for processing READ 1 only ###
###################################################

# Note: hard-coded file paths to be replaced with cmd args

import subprocess
import os
import re
import DBR_Parsing
from DBR_Parsing import R1_dict, DBR_dict

def PEAR_assemble(in_dir, forward, reverse, out_dir, out_name,  extra_params=None):
    print 'Merging overlapping reads 1 & 2 with PEAR.'
    # extra params is a list of optional args to pass to PEAR
    if extra_params:
        pear_call = 'pear' + ' -f ' + in_dir + forward + ' -r ' + in_dir + reverse + ' -o ' + out_dir + out_name + ' ' + extra_params
    else:
        pear_call = 'pear' + ' -f ' + in_dir + forward + ' -r ' + in_dir + reverse + ' -o ' + out_dir + out_name
    subprocess.call(pear_call, shell=True)

def FASTQ_R1_R2_merge(in_dir, fq_r1, fq_r2, fq_out):
    print 'Taking reverse complement of read 2 with FASTX Toolkit.\n'
    rc_call = 'fastx_reverse_complement -i ' + in_dir + fq_r2 + ' -o ' + in_dir + fq_out
    subprocess.call(rc_call, shell=True)
    
    print 'Merging read 1 and read 2 reverse complement into a single FASTQ file.'
    merge_call = 'paste ' + in_dir + fq_r1 + ' ' + in_dir + fq_out + ' > ' + in_dir + fq_out + '_merged.fq'
    subprocess.call(merge_call, shell=True)
    rm_tab_call = "sed 's/[ \t]//g' " + in_dir + fq_out + '_merged.fq' + ' >> ' + in_dir + fq_out + '_untabbed_merged.fq'
    subprocess.call(rm_tab_call, shell=True)
    
def FASTQ_quality_filter(fq_in, fq_out, q, p):    
    print 'Quality filtering with FASTX Toolkit.\n'
    if 'gz' in fq_in:
        fqc_call = "zcat " + fq_in + " | fastq_quality_filter -q " + str(q) + " -p " + str(p) + " -z -o " + fq_out 
    else:
        fqc_call = "fastq_quality_filter -q " + str(q) + " -p " + str(p) + " -z -i " + fq_in + " -o " + fq_out 
    print fqc_call
    subprocess.call(fqc_call, shell=True)
    # to-do: use zcat | wc -l to determine how many lines were removed after filtering (print to screen)

def Demultiplex(in_file, barcode_file, out_dir, out_prefix):    
    print 'Demultiplexing sequence data with FASTX Toolkit.\n'
    
    prefix_path = out_dir + out_prefix
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    demultiplex_call = 'zcat ' + in_file + ' | /usr/local/bin/fastx_barcode_splitter.pl --bcfile ' + barcode_file + ' --prefix ' + prefix_path + ' --bol'
    
    subprocess.call(demultiplex_call, shell=True)

def Trim(in_dir, out_dir, suffix, first_base, last_base=None):    
    print 'Trimming DBR and enzyme cut sites with FASTX Toolkit.\n'
    
    # new directory for trimmed files
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for i in os.listdir(in_dir):
        # the proper file name
        proper_file_trim = i+suffix
    
        # full path to file for trimming
        full_path=os.path.join(in_dir, i)
        new_path=os.path.join(out_dir, proper_file_trim)
    
        # Remove barcodes and R1 enzyme cut site (first
        trim_call = "/usr/local/bin/fastx_trimmer -f " + str(first_base) + ' -l ' + str(last_base) + " -i " + full_path + " -o " + new_path
        subprocess.call(trim_call, shell=True)

def denovo_Stacks(in_dir, denovo_path, stacks_executables, out_dir, m, n, b, D):    
    print 'Assembling sequences de novo using Stacks denovo_map.pl\n'
    
    # Generate "-s" arguments to Stacks denovo_map.pl; trim out barcode and cut site sequences
    s_list=[]
    rm_unmatched = False
    for i in os.listdir(in_dir):
        # don't proceed if the inputs don't have the proper extension
        assert '.fq' in i, "%s needs extension .fq for Stacks compatibility" % i

        if 'unmatched' in i: # don't process the unmatched sequences
            print 'Skipping sequences with unmatched barcodes.\n'
            rm_unmatched = True
        else:
            # full path to trimmed file
            new_path=os.path.join(in_dir, i)
    
            # generate a string for the -s stacks argument
            s_list.append('-s ')
            s_list.append(new_path)
            s_list.append(' ')
    
    # the last file in the dir is "unmatched"... remove it and the final trailing space
    #s_list = s_list[0:len(s_list)-4]
    
    # join list into a string to pass to denovo_map.pl
    formatted_list = ''.join(s_list)

    # Run denovo_map.pl
    build_args = [denovo_path, ' -e ', stacks_executables, ' -o ', out_dir, ' -m ', str(m), ' -n ', str(n), ' -t ', '-b ', str(b), ' -D ', D, ' -S ', formatted_list]
    denovo_call = ''.join(build_args)
    subprocess.call(denovo_call, shell=True)

def GeneratePseudoref(in_dir, out_dir, out_name, BWA_path):    
    print "Preparing pseudoreference genome by extracting Stacks denovo consensus sequence.\n"

    ref_path = in_dir + 'batch_1.catalog.tags.tsv'
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

    out_path=os.path.join(out_dir, out_name)
    print "Writing pseudoreference genome to file " + out_path + '\n'
    with open(out_path, 'w') as fq: # need to check if it already exists
        fq.write(fastq_interleaved)
    
    # Index the new pseudoreference genome
    print 'Indexing the reference genome using BWA.\n'
    index_call =  BWA_path + ' index ' + out_path # index = BWA function to use
    print index_call
    subprocess.call(index_call, shell = True)

def refmap_BWA(in_dir, out_dir, BWA_path, pseudoref_full_path):    
    
    print 'Mapping sequence data to pseudoreference genome using BWA.\n'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for i in os.listdir(in_dir): # for every demultiplexed & trimmed fastq...
        rex = re.compile(r'\d+') # the regex should be a cmd line argument probably, since the sample name is likely to vary....
        if rex.search(i):
            fname, fext = os.path.splitext(i)
            print 'Reference mapping ' + fname + '\n'
            read_group_header = '"@RG\\tID:' + fname + '\\tPL:Illumina\\tLB:' + fname + '"' 
            bwa_mem_call = BWA_path + ' mem -M -R ' + read_group_header + " " + pseudoref_full_path + ' ' + in_dir + i + ' > ' + out_dir + fname + '.sam'
            print bwa_mem_call
            subprocess.call(bwa_mem_call, shell=True)

### Test functions
'''
# ASSEMBLE WITH PEAR (OPTIONAL)
in_dir = '/home/antolinlab/Downloads/'
forward = 'Sample1_ACTTGA_L008_R1_001.fastq.gz'
reverse = 'Sample1_ACTTGA_L008_R2_001.fastq.gz'
out_dir = in_dir
out_name = 'pear_merged_Sample1_ACTTGA'
PEAR_assemble(in_dir, forward, reverse, out_dir, out_name)
'''
'''
# MERGE READ 1 AND READ 2
in_dir = '/home/antolinlab/Downloads/'
fq_r1 = 'Sample1_ACTTGA_L008_R1_001.fastq'
fq_r2 = 'Sample1_ACTTGA_L008_R2_001.fastq'
fq_out = 'Sample1_ACTTGA_R1_R2rc'
FASTQ_R1_R2_merge(in_dir, fq_r1, fq_r2, fq_out)
'''
'''
# FILTER READ 1 FOR ASSEMBLY
fq_in = '/home/antolinlab/Downloads/Sample1_ACTTGA_L008_R1_001.fastq.gz'
q = 25
p = 50
fq_out = '/home/antolinlab/Downloads/Sample1_ACTTGA_L008_R1_001_filtered.fastq.gz'
FASTQ_quality_filter(fq_in, fq_out, q, p)

# FILTER READ 2 FOR DBR PARSING
fq_in = '/home/antolinlab/Downloads/Sample1_ACTTGA_L008_R2_001.fastq.gz'
q = 25
p = 50
fq_out = '/home/antolinlab/Downloads/Sample1_ACTTGA_L008_R2_001_filtered.fastq.gz'
FASTQ_quality_filter(fq_in, fq_out, q, p)

# FILTER MERGED READS 1 AND 2
fq_in = '/home/antolinlab/Downloads/Sample1_ACTTGA_R1_R2rc_untabbed_merged.fq'
q = 25
p = 50
fq_out = '/home/antolinlab/Downloads/Sample1_ACTTGA_R1_R2rc_untabbed_merged_filtered.fq'
FASTQ_quality_filter(fq_in, fq_out, q, p)

# BUILD DBR PARSING DICTIONARIES
R1_dict('/home/antolinlab/Downloads/Sample1_ACTTGA_R1_R2rc_untabbed_merged_filtered.fq', test_dict = True, save_path='/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/initial_qualFilter_R1_dict')
DBR_dict('/home/antolinlab/Downloads/Sample1_ACTTGA_R1_R2rc_untabbed_merged_filtered.fq', dbr_start = 231, dbr_stop = 238, test_dict = True, save_path='/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/initial_qualFilter_dbr_dict')  

# DEMULTIPLEX MERGED READS
in_file = '/home/antolinlab/Downloads/Sample1_ACTTGA_R1_R2rc_untabbed_merged_filtered.fq'
mp_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed/'
prefix = 'pilot_demultiplex_'
barcode_file = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/pilot_barcode_file'
Demultiplex(in_file, barcode_file, mp_dir, prefix)

# TRIM DEMULTIPLEXED READS
in_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed/'
out_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed_Trimmed/'
suffix = '_trimmed.fq' # eventually make cmd arg
first_base = 11
last_base = 238
Trim(in_dir, out_dir, suffix, first_base, last_base)

# INITIAL DENOVO ASSEMBLY
in_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed_Trimmed/'
denovo_path = '/home/antolinlab/Downloads/stacks-1.31/scripts/denovo_map.pl '
stacks_executables = '/home/antolinlab/Downloads/stacks-1.31/scripts'
out_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/'
denovo_Stacks(in_dir, denovo_path, stacks_executables, out_dir, m=10, n=2, b=1, D= 'initial_assembly')

# CREATE PSEUDOREFERENCE
in_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/'
out_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/'
BWA_path = '/home/antolinlab/Downloads/bwa.kit/bwa'
out_name = 'pseudoref.fa'
GeneratePseudoref(in_dir, out_dir, out_name, BWA_path)
'''
# INITIAL PSEUDOREFERENCE-MAPPED ASSEMBLY
out_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Assembled/'
in_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed_Trimmed/'
BWA_path = '/home/antolinlab/Downloads/bwa.kit/bwa'
pseudoref_full_path = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/pseudoref.fa'
refmap_BWA(in_dir, out_dir, BWA_path, pseudoref_full_path)