'''
Created on May 12, 2016

@author: Kelly
'''
#!/usr/bin/python

#############################################################################
### De novo pipeline for processing CWD mule deer dataset WITH DUPLICATES ###
#############################################################################

#### IMPORTS ####
from integrated_denovo_pipeline import *
from DBR_Parsing import *
from assembled_DBR_filtering import *

if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-t', '--threads', help = 'Number of threads to use for parallel operations.')
    parser.add_argument('-p', '--parent_directory', help = 'Path to directory containing the raw sequence data.')
    parser.add_argument('-b', '--barcodes', help = 'Path to directory containing the barcode sequences for each library.')
    opts = parser.parse_args()
    
    threads=int(opts.threads)
    
    #########################################################################
    ### SET UP THREADS TO USE FOR PARALLEL PROCESSES                      ###
    #########################################################################
    
    for i in xrange(threads):
        t = Thread(target=worker)
        t.daemon = True
        t.start()
    
    #########################################################################
    ### INPUTS AND OUTPUTS                                                ###
    #########################################################################
    
    
    #### PART 1: INITIAL ASSEMBLY AND DBR FILTERING
    # user only needs to specify parent directory; the remaining directories should be automatically generated
    parentDir = '/home/pierce/CWD_RADseq/raw_new/'
    
    # The remaining directories are automatically generated from the parent directory
    # The following QC steps have been performed in earlier runs of the script and do not need to be redone:
    #    - PEAR merging of raw sequence reads
    #    - Quality filtering of PEAR merged sequences
    #    - Generation of DBR dictionaries from merged, quality filtered sequences
    # Code for these QC steps is included but commented out and not executed with this run.
    #pearInDir = parentDir
    #pearOutDir = parentDir + '/pear_merged_parallel/'
    #filterInDir = pearOutDir
    filterOutDir =  parentDir + '/qual_filtered/'
    #dbrInDir = filterOutDir
    dbrOutDir = parentDir + '/DBR_dir/'
    demultiplexInDir = filterOutDir
    demultiplexOutDir = parentDir + '/demultiplexed_duplicates_combined/'
    trimInDir = demultiplexOutDir
    trimOutDir = parentDir + '/trimmed_duplicates_combined/'
    stacksInDir = trimOutDir
    stacksOutDir = parentDir + '/StacksOutput_duplicates_combined/' # stacks doesn't allow an output to be specified
    pseudorefInDir = stacksOutDir
    pseudorefOutDir = parentDir + '/pseudoreference_duplicates_combined.fastq'
    BWAoutDir = parentDir + '/BWA_duplicates_combined/'
    DBRfilteredseqs = parentDir + '/dbrFiltered_duplicates_combined/'
    
    #### PART 2: REASSEMBLING THE FILTERED SEQUENCES
    re_demultiplexInDir = DBRfilteredseqs
    re_demultiplexOutDir = parentDir + '/dbrFiltered_demultiplexed_duplicates_combined/'
    re_BWAinDir = DBRfilteredseqs
    re_BWAoutDir = parentDir + '/dbrFiltered_BWA_duplicates_included/'
    finalBCFout = parentDir + '/dbrFiltered_pseudorefMapped_genotypes_duplicates_combined.bcf'
    finalVCFout = parentDir + '/dbrFiltered_pseudorefMapped_genotypes_duplicates_combined.vcf'
    
    #########################################################################
    ### FUNCTION CALLS TO RUN THE PIPELINE                                ###
    #########################################################################
    
    # ASSEMBLE ITERATIVELY WITH PEAR 
    #out_name = 'pear_merged_'
    #extra_params = '-m 309 -n 209'
    #parallel_PEAR_assemble(in_dir = pearInDir, 
    #                        out_dir = pearOutDir, 
    #                        out_name = out_name, 
    #                        extra_params = extra_params,
    #                        regexR1='R1', regexR2='R2')
    
    # QUALITY FILTER PEAR ASSEMBLED DATA
    #out_name = '.qual_filtered' # gets appended to input file name
    #q = 30
    #p = 50
    #read = '.assembled.fastq' # extension for pear-assembled reads
    #parallel_FASTQ_quality_filter(directory = filterInDir, 
    #                               out_dir = filterOutDir, 
    #                               out_name = out_name, 
    #                               q = q, 
    #                               p = p, 
    #                               read = read)
    
    # MAKE DBR DICTIONARIES FOR QUAL FILTERED PEAR DATA
    #seq_type = 'pear'
    #parallel_DBR_dict(in_dir = dbrInDir, 
    #                   seqType = seq_type,
    #                   save = dbrOutDir,
    #                   dbr_start = -10,
    #                   dbr_stop = -2)
    
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
    
    # RUN USTACKS SIMULTANEOUSLY ON ALL LIBRARIES
    denovo_Ustacks(in_dir = stacksInDir, 
                  denovo_path = denovo_path, 
                  stacks_executables = stacks_executables, 
                  out_dir = stacksOutDir, 
                  m = 10, 
                  n = 2, 
                  b = 1, 
                  D = '_initial_assembly')
    
    # RUN CSTACKS SIMULTANEOUSLY ON ALL LIBRARIES (same args as above)
    denovo_Cstacks(in_dir = stacksInDir, 
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
    parallel_refmap_BWA(in_dir = trimOutDir, # input demultiplexed, trimmed reads
               out_dir = BWAoutDir, 
               BWA_path = BWA, # imported from integrated_denovo_pipeline.py 
               pseudoref_full_path = pseudorefOutDir)
        
    # FILTER OUT PCR DUPLICATES USING THE DBR SEQUENCES
    DBR_Filter(assembled_dir = BWAoutDir, # the SAM files for the data mapped to pseudoreference
               out_dir = DBRfilteredseqs, # the output file, full path, ending with .fasta
               n_expected = 2, # the number of differences to be tolerated
               barcode_dir = '/home/pierce/CSU_ChronicWasting/Barcodes_WithDuplicates', # the barcodes for individuals in the library referenced in dict_in
               dict_dir = dbrOutDir, # a single dictionary of DBRs (for one library only)
               sample_regex = '.*_(\d{1,3}T?)_.*',
               barcode_file=None, # if just a single library is being used, can directly pass the barcode file
               test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
               phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
               samMapLen=None)
    
    # REFERENCE MAP DBR FILTERED READS TO THE PSEUDOREFERENCE
    parallel_refmap_BWA(in_dir = re_BWAinDir, # input demultiplexed, trimmed reads
               out_dir = re_BWAoutDir, 
               BWA_path = BWA, # imported from integrated_denovo_pipeline.py 
               pseudoref_full_path = pseudorefOutDir)
    
    # CALL THE GENOTYPES USING SAMTOOLS MPILEUP; CONVERT OUTPUT TO VCF FILE
    callGeno(sam_in = re_BWAoutDir, 
             pseudoref = pseudorefOutDir, 
             BCFout = finalBCFout, 
             VCFout = finalVCFout,
             samtoolsPath = samtoolsPath,
             bcftoolsPath = bcftoolsPath)
