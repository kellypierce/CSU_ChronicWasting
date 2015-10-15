#!/usr/bin/python

#########################################################################
### De novo pipeline for processing READ 1 only -- PILOT LIBRARY TEST ###
#########################################################################

from integrated_denovo_pipeline import *
from DBR_Parsing import *

### Test functions
# ASSEMBLE SINGLE SET OF PAIRED READS WITH PEAR 
#in_dir = '/home/antolinlab/Downloads/CWD_RADseq/'
#forward = 'Library12_S65_L008_R1_001.fastq.gz'
#reverse = 'Library12_S65_L008_R2_001.fastq.gz'
#out_dir = in_dir
#out_name = 'pear_merged_'
#extra_params = '-m 280 -n 150'
#PEAR_assemble(in_dir, forward, reverse, out_dir, out_name, extra_params)

# ASSEMBLE ITERATIVELY WITH PEAR 
#in_dir = '/home/antolinlab/Downloads/CWD_RADseq/'
#out_dir = '/home/antolinlab/Downloads/CWD_RADseq/pear_merged/'
#out_name = 'pear_merged_'
#extra_params = '-m 280 -n 150'
#iterative_PEAR_assemble(in_dir, out_dir, out_name, extra_params, r1='R1', r2='R2')

# count up the sequence lengths
# sed -n '2~4p' pear_merged_Library12_L8.assembled.fastq | awk '{print length}' >> pear_merged_Library12_L8_assembled_seq_length.txt

# count up the R1 cutsites
#sed -n '2~4p' pear_merged_Library12_L8.assembled.fastq | cut -c 6-10 | sort | uniq -c | sort -nr -k 1 >> pear_merged_Library12_L8_assembled_R1cut_check.txt

# count up the R2 cutsites
#sed -n '2~4p' pear_merged_Library12_L8.assembled.fastq | rev | cut -c 11-14 | sort | uniq -c | sort -nr -k 1 >> pear_merged_Library12_L8_assembled_R2cut_check.txt

# QUALITY FILTER PEAR ASSEMBLED DATA
directory = '/home/antolinlab/Downloads/CWD_RADseq/'
out_dir = '/home/antolinlab/Downloads/CWD_RADseq/qual_filtered/'
out_name = '.qual_filtered' # gets appended to input file name
q = 30
p = 50
read = '.assembled.fastq' # extension for pear-assembled reads
iterative_FASTQ_quality_filter(directory, out_dir, out_name, q, p, read=read)
'''
# MAKE DBR DICTIONARIES FOR QUAL FILTERED PEAR DATA
fq_path = '/home/antolinlab/Downloads/CWD_RADseq/qual_filtered/'
iterative_DBR_dict(fq_path, seqType, read, save_path)
'''