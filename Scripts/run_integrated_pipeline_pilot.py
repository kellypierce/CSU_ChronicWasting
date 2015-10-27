#!/usr/bin/python

#########################################################################
### De novo pipeline for processing READ 1 only -- PILOT LIBRARY TEST ###
#########################################################################

from integrated_denovo_pipeline import *

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
'''
### PHASE 2: DENOVO ASSEMBLY OF DBR FILTERED DATA
in_file = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/DBR_filtered_sequences.fastq'
mp_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed_Phase2/'
prefix = 'pilot_demultiplex_phase2_'
barcode_file = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/pilot_barcode_file'
#Demultiplex(in_file, barcode_file, mp_dir, prefix)

in_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed_Phase2/'
out_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed_Phase2_Trimmed/'
suffix = '_phase2_trimmed.fq' # eventually make cmd arg
first_base = 6
#Trim(in_dir, out_dir, suffix, first_base)

in_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed_Phase2_Trimmed/'
denovo_path = '/home/antolinlab/Downloads/stacks-1.31/scripts/denovo_map.pl '
stacks_executables = '/home/antolinlab/Downloads/stacks-1.31/scripts'
out_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/'
denovo_Stacks(in_dir, denovo_path, stacks_executables, out_dir, m=10, n=2, b=1, D= 'phase2_assembly')

'''
# CREATE PSEUDOREFERENCE
in_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/'
out_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/'
BWA_path = '/home/antolinlab/Downloads/bwa.kit/bwa'
out_name = 'pseudoref.fa'
GeneratePseudoref(in_dir, out_dir, out_name, BWA_path)

# INITIAL PSEUDOREFERENCE-MAPPED ASSEMBLY
out_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Assembled/'
in_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Demultiplexed_Trimmed/'
BWA_path = '/home/antolinlab/Downloads/bwa.kit/bwa'
pseudoref_full_path = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/pseudoref.fa'
#refmap_BWA(in_dir, out_dir, BWA_path, pseudoref_full_path)
'''