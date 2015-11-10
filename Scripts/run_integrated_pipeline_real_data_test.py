#!/usr/bin/python

#########################################################################
### De novo pipeline for processing READ 1 only -- PILOT LIBRARY TEST ###
#########################################################################

from integrated_denovo_pipeline import *
from DBR_Parsing import *
from setuptools.extension import Library

### Test functions
# ASSEMBLE SINGLE SET OF PAIRED READS WITH PEAR 
#in_dir = '/home/antolinlab/Downloads/CWD_RADseq/'
#forward = 'Library12_S65_L008_R1_001.fastq.gz'
#reverse = 'Library12_S65_L008_R2_001.fastq.gz'
#out_dir = in_dir
#out_name = 'pear_merged_'
#extra_params = '-m 309 -n 209'
#PEAR_assemble(in_dir, forward, reverse, out_dir, out_name, extra_params)

# POST-PEAR QC CHECKS FOR LIBRARY 12 ONLY
# count up the sequence lengths
# sed -n '2~4p' pear_merged_Library12_S65_L008__001.fastq.assembled.fastq | awk '{print length}' >> pear_merged_Library12_L8_assembled_seq_length.txt
# count up the R1 cutsites
#sed -n '2~4p' pear_merged_Library12_L8.assembled.fastq | cut -c 6-10 | sort | uniq -c | sort -nr -k 1 >> pear_merged_Library12_L8_assembled_R1cut_check.txt
# count up the R2 cutsites
#sed -n '2~4p' pear_merged_Library12_L8.assembled.fastq | rev | cut -c 11-14 | sort | uniq -c | sort -nr -k 1 >> pear_merged_Library12_L8_assembled_R2cut_check.txt

# PATHS TO EXECUTABLES
denovo_path = '/home/antolinlab/Downloads/stacks-1.31/scripts/denovo_map.pl '
stacks_executables = '/home/antolinlab/Downloads/stacks-1.31/scripts'

# PATHS TO INPUTS AND OUTPUTS
# user only needs to specify parent directory; the remaining directories should be automatically generated
parentDir = '/home/antolinlab/Downloads/CWD_RADseq/'
pearInDir = parentDir
pearOutDir = parentDir + '/pear_merged/'
filterInDir = pearOutDir
filterOutDir =  parentDir + '/qual_filtered/'
dbrInDir = filterOutDir
dbrOutDir = parentDir + '/DBR_dir/'
demultiplexInDir = filterOutDir
demultiplexOutDir = parentDir + '/demultiplexed/'
trimInDir = demultiplexOutDir
trimOutDir = parentDir + '/trimmed/'
stacksInDir = trimOutDir
stacksOutDir = parentDir + '/StacksOutput/' # stacks doesn't allow an output to be specified
pseudorefInDir = stacksOutDir
pseudorefOutDir = parentDir
BWAinDir = parentDir
BWAoutDir = parentDir + '/BWA/'

# ASSEMBLE ITERATIVELY WITH PEAR 
out_name = 'pear_merged_'
extra_params = '-m 309 -n 209'
iterative_PEAR_assemble(in_dir = pearInDir, 
                        out_dir = pearOutDir, 
                        out_name = out_name, 
                        extra_params = extra_params,
                        regexR1='R1', regexR2='R2')

# QUALITY FILTER PEAR ASSEMBLED DATA
out_name = '.qual_filtered' # gets appended to input file name
q = 30
p = 50
read = '.assembled.fastq' # extension for pear-assembled reads
iterative_FASTQ_quality_filter(directory = filterInDir, 
                               out_dir = filterOutDir, 
                               out_name = out_name, 
                               q = q, 
                               p = p, 
                               read = read)

# MAKE DBR DICTIONARIES FOR QUAL FILTERED PEAR DATA
seq_type = 'pear'
iterative_DBR_dict(in_dir = dbrInDir, 
                   seqType = seq_type,
                   save = dbrOutDir,
                   dbr_start = -10,
                   dbr_stop = -2)

# DEMULTIPLEX
out_prefix = '/demultiplexed_'
iterative_Demultiplex(in_dir = demultiplexInDir, 
                      barcode_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/BarcodesRound1/', 
                      out_dir = demultiplexOutDir, 
                      out_prefix = out_prefix)

# TRIM TO UNIFORM LENGTH
suffix = '_trimmed.fq'
first_base = 11
last_base = 196
Trim(in_dir = trimInDir, 
     out_dir = trimOutDir, 
     suffix = suffix, 
     first_base = first_base, 
     last_base = last_base)

# RUN STACKS SIMULTANEOUSLY ON ALL LIBRARIES
denovo_Stacks(in_dir = stacksInDir, 
              denovo_path = denovo_path, 
              stacks_executables = stacks_executables, 
              out_dir = stacksOutDir, 
              m = 10, 
              n = 2, 
              b = 1, 
              D = '_initial_assembly')
'''
# GENERATE THE PSEUDOREFERENCE GENOME
GeneratePseudoref(in_dir, out_dir, out_name, BWA_path)

# REFERENCE MAP QUALITY FILTERED/DEMULTIPLEXED MERGED READS TO THE PSEUDOREFERENCE
refmap_BWA(in_dir, out_dir, BWA_path, pseudoref_full_path)

# FILTER THE DBRS -- THIS REQUIRES A REVISION TO CREATE A FUNCTION IN assembled_DBR_filtering
filterDBR(sample_list, dbrDict)

sample_list = the samples from a given Library
dbrDict = the DBR dictionary generated from that Library
'''