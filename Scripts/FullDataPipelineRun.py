#!/usr/bin/python

#########################################################################
### De novo pipeline for processing complete CWD mule deer dataset    ###
#########################################################################

#### IMPORTS ####
from integrated_denovo_pipeline import *
from DBR_Parsing import *
from assembled_DBR_filtering import *

#########################################################################
### NOTES ON THIS RUN                                                 ###
#########################################################################

# 1. All samples were merged with pear previously
#    - The pear merge function is in this script for posterity's sake, but not executed
#    - The pear merged files were instead copied into the directory for this analysis
#    - The decision not to re-run pear was based on the fact that the output is the same regardless, and it takes ages to run
# 2. The raw .fastq.gz and pear .fastq file names have been changed
#    - The pattern 'LibraryXXX' must be present in the file name for the script to recognize the samples
#    - Since not all data came off the sequencer with that information, names needed to be changed
#    - Other information in the file name was left unchanged
# 3. Libraries 12 and 120 were split across two lanes each
#    - This means individual samples with the same barcode (NOT technical replicates) exist in two .fastq.gz files (and assembled pear files)
#    - This further means that individual read counts in a single library file will be low coverage -- the desired coverage is obtained only when looking across data from both lanes
#    - The pipeline script cannot handle that scenario explicity and will overwrite sample data from the first lane processed with sample data from the second lane processed
#    - I circumvent this problem by concatenating the files, which are renamed accordingly. This was done with the pear merged files, but not with the original data files
#    - The un-concatenated data are removed from the pear merged files directory

#########################################################################
### INPUTS AND OUTPUTS                                                ###
#########################################################################

# user only needs to specify parent directory; the remaining directories should be automatically generated

#### PART 1: INITIAL ASSEMBLY AND DBR FILTERING
parentDir = '/home/pierce/CWD_RADseq/raw_new/'
pearInDir = parentDir
pearOutDir = parentDir + '/pear_merged_parallel/'
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
pseudorefOutDir = parentDir + '/pseudoreference.fastq'
BWAinDir = parentDir
BWAoutDir = parentDir + '/BWA/'
DBRfilteredseqs = parentDir + '/dbrFiltered/'

#### PART 2: REASSEMBLING THE FILTERED SEQUENCES
re_demultiplexInDir = DBRfilteredseqs
re_demultiplexOutDir = parentDir + '/dbrFiltered_demultiplexed/'
re_trimInDir = re_demultiplexOutDir
re_trimOutDir = parentDir + '/dbrFiltered_trimmed/'
re_stacksInDir = re_trimOutDir
re_stacksOutDir = parentDir + '/dbrFiltered_StacksOutput/' # stacks doesn't allow an output to be specified
re_pseudorefInDir = re_stacksOutDir
re_pseudorefOutDir = parentDir + '/dbrFiltered_pseudoreference.fastq'
re_BWAinDir = parentDir
re_BWAoutDir = parentDir + '/dbrFiltered_BWA2/'
finalBCFout = parentDir + '/dbrFiltered_pseudorefMapped_genotypes2.bcf'
finalVCFout = parentDir + '/dbrFiltered_pseudorefMapped_genotypes2.vcf'

#########################################################################
### FUNCTION CALLS TO RUN THE PIPELINE                                ###
#########################################################################

# ASSEMBLE ITERATIVELY WITH PEAR 
out_name = 'pear_merged_'
extra_params = '-m 309 -n 209'
#iterative_PEAR_assemble(in_dir = pearInDir, 
#                        out_dir = pearOutDir, 
#                        out_name = out_name, 
#                        extra_params = extra_params,
#                        regexR1='R1', regexR2='R2')

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
                      barcode_dir = '/home/pierce/CSU_ChronicWasting/RevisedBarcodes', 
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

# GENERATE THE PSEUDOREFERENCE GENOME
GeneratePseudoref(in_dir = pseudorefInDir, 
                  out_file = pseudorefOutDir,  
                  BWA_path = BWA) # imported from integrated_denovo_pipeline.py

# REFERENCE MAP QUALITY FILTERED/DEMULTIPLEXED MERGED READS TO THE PSEUDOREFERENCE
refmap_BWA(in_dir = trimOutDir, # input demultiplexed, trimmed reads
           out_dir = BWAoutDir, 
           BWA_path = BWA, # imported from integrated_denovo_pipeline.py 
           pseudoref_full_path = pseudorefOutDir)


DBR_Filter(assembled_dir = BWAoutDir, # the SAM files for the data mapped to pseudoreference
           out_dir = DBRfilteredseqs, # the output file, full path, ending with .fasta
           n_expected = 2, # the number of differences to be tolerated
           barcode_dir = '/home/pierce/CSU_ChronicWasting/BarcodesRound1/', # the barcodes for individuals in the library referenced in dict_in
           dict_dir = dbrOutDir, # a single dictionary of DBRs (for one library only)
           barcode_file=None, # if just a single library is being used, can directly pass the barcode file
           test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
           phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
           samMapLen=None)

# DEMULTIPLEX
out_prefix = '/re_demultiplexed_'
iterative_Demultiplex(in_dir = re_demultiplexInDir, 
                      barcode_dir = '/home/pierce/CSU_ChronicWasting/RevisedBarcodes', 
                      out_dir = re_demultiplexOutDir, 
                      out_prefix = out_prefix)

# TRIM TO UNIFORM LENGTH
suffix = '_re_trimmed.fq'
new_first_base = 6
Trim(in_dir = re_trimInDir, 
     out_dir = re_trimOutDir, 
     suffix = suffix, 
     first_base = new_first_base)

# RUN STACKS SIMULTANEOUSLY ON ALL LIBRARIES
denovo_Stacks(in_dir = re_stacksInDir, 
              denovo_path = denovo_path, 
              stacks_executables = stacks_executables, 
              out_dir = re_stacksOutDir, 
              m = 10, 
              n = 2, 
              b = 1, 
              D = '_final_assembly')

# GENERATE THE PSEUDOREFERENCE GENOME
GeneratePseudoref(in_dir = re_pseudorefInDir, 
                  out_file = re_pseudorefOutDir,  
                  BWA_path = BWA) # imported from integrated_denovo_pipeline.py

# REFERENCE MAP QUALITY FILTERED/DEMULTIPLEXED MERGED READS TO THE PSEUDOREFERENCE
refmap_BWA(in_dir = re_trimOutDir, # input demultiplexed, trimmed reads
           out_dir = re_BWAoutDir, 
           BWA_path = BWA, # imported from integrated_denovo_pipeline.py 
           pseudoref_full_path = re_pseudorefOutDir)

# CALL THE GENOTYPES USING SAMTOOLS MPILEUP; CONVERT OUTPUT TO VCF FILE
callGeno(sam_in = re_BWAoutDir, 
         pseudoref = re_pseudorefOutDir, 
         BCFout = finalBCFout, 
         VCFout = finalVCFout)
