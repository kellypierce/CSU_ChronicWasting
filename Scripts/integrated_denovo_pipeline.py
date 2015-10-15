#!/usr/bin/python

###################################################
### De novo pipeline for processing READ 1 only ###
###################################################

# Note: hard-coded file paths to be replaced with cmd args

import subprocess
from subprocess import call, Popen, PIPE
import os as os
from os import linesep, path, R_OK, X_OK
import re
from DBR_Parsing import *
import string as string
from string import Template, join
import logging as logging
from logging import debug, critical, error, info
import sys as sys
import warnings
import fnmatch

def configureLogging(verbose = False):
    '''
    setup the logger
    '''
    logger = logging.getLogger()
    logger.handlers = []

    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    
    consoleLogHandler = logging.StreamHandler()

    # do not need colors if logging to a non-shell
    if sys.stdout.isatty():
        consoleLogHandler.setFormatter(logging.Formatter("\033[93m%(filename)s:%(lineno)s\033[0m: %(message)s"))
    else:
        consoleLogHandler.setFormatter(logging.Formatter("%(message)s"))

    logger.addHandler(consoleLogHandler)

def checkFile(filename):
    '''
    return true if this is a file and is readable on the current filesystem
    '''
    try:
        if os.path.exists(os.path.abspath(filename)) and os.path.isfile(os.path.abspath(filename)) and os.access(os.path.abspath(filename), R_OK):
            return True
        fullPath = string.join(os.getcwd(), filename[1:])
        return os.path.exists(fullPath) and os.path.isfile(fullPath) and os.access(fullPath, R_OK)
    except IOError:
        return False

def checkExe(filename):
    '''
    return true if this is an executable file on the current filesystem
    '''
    if not isinstance(filename, str):
        raise TypeError("need a string, got a %s" % type(filename))
    return (os.path.exists(filename) and os.path.isfile(filename) and os.access(filename, X_OK))

################################## GLOBALS ####################################

# run configureLogging first
configureLogging(False)

# PEAR assembly
pearPath = '/home/antolinlab//Downloads/PEAR/src/pear'
pearSimpleTemplate = Template('%s -f $f -r $r -o $o' % pearPath)
pearExtraParamsTemplate = Template('%s -f $f -r $r -o $o $e' % pearPath)

# Fastq filtering of assembled data
#qualityFilter = '/Users/Kelly/Downloads/bin/fastq_quality_filter'
qualityFilter = '/home/antolinlab/Downloads/fastx_toolkit-0.0.14/src/fastq_quality_filter/fastq_quality_filter'
fqfStdinTemplate = Template('%s -q $q -p $p -Q33 -z -o $output' % qualityFilter)
fqfFileTemplate = Template('%s -q $q -p $p -Q33 -z -i $input -o $output' % qualityFilter)

# DBR extraction
#dbrStdinTemplate =
#dbrFileTemplate =

# Demultiplexing
#demultiplexStdinTemplate =
#demultiplexFileTemplate =

# Trimming

# Stacks de novo assembly

# Extract reference from Stacks consensus, index

# Reference-based assembly

# DBR filter


###############################################################################

def iterative_PEAR_assemble(in_dir, out_dir, out_name, extra_params, r1='*', r2='*'):
    files = os.listdir(in_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for f in files:
        if re.findall(r1, f):
            forward = f
            reverse = re.sub(r1, r2, f)
            PEAR_assemble(in_dir, forward, reverse, out_dir, out_name, extra_params)

def PEAR_assemble(in_dir, forward, reverse, out_dir, out_name,  extra_params=None):
    if not checkFile(in_dir + forward):
        raise IOError("Where is the forward read file: %s" % forward)
    if not checkFile(in_dir + reverse):
        raise IOError("where is the reverse read file: %s" % reverse)
    
    if not forward.endswith(".fastq.gz"):
        warnings.warn("Expect raw sequence data in .fastq.gz format")
    if not reverse.endswith(".fastq.gz"):
        warnings.warn("Expect raw sequence data in .fastq.gz format")
        
    if not checkExe(pearPath):
        raise Exception("could not find %s in the filesystem for execution, is the environment setup correctly?" % pearPath)
    info('Merging overlapping reads %s and %s with PEAR.' % (forward, reverse))
    
    out_split = os.path.splitext(forward)
    out_suffix = re.sub(r'R\d{1}', '', out_split[0])
    out_full = out_dir + out_name + out_suffix
    
    info('Saving outputs to %s' % out_full)
    
    if extra_params:
        commandLine = pearExtraParamsTemplate.substitute(f = in_dir + forward,
                                                         r = in_dir + reverse,
                                                         o = out_full,
                                                         e = extra_params)
    else:
        commandLine = pearSimpleTemplate.substitute(f = in_dir + forward,
                                                         r = in_dir + reverse,
                                                         o = out_full)
    pearProcess = Popen(commandLine, shell=True)
    pearProcess.wait()

def iterative_FASTQ_quality_filter(directory, out_dir, out_name, q, p, read='*'):
    files = os.listdir(directory)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    info('Quality filtering files containing %s' % read)
    info('Results saved to %s' % out_dir)
    for f in files:
        r2 = re.findall(read, f)
        if r2: 
            fileRoot = os.path.splitext(f)
            out_file = directory + fileRoot[0] + out_name
            in_file = directory + f
            FASTQ_quality_filter(in_file, out_file, q, p)
 
def FASTQ_quality_filter(fq_in, fq_out, q, p, qualityFilter = qualityFilter):
    if not checkFile(fq_in):
        raise IOError("where is the input file: %s" % fq_in)
    if not fq_in.endswith("gz"):
        warnings.warn("prefer to pass compressed files, consider gzipping the inputs")
    if not checkExe(qualityFilter):
        raise Exception("could not find %s in the filesystem for execution, is the environment setup correctly?" % qualityFilter)
    info("Quality filtering with FASTX Toolkit. Output file %s" % fq_out)
    if fq_in.endswith('gz'): # chain a gz decompressor thread to fqf
        commandLine = fqfStdinTemplate.substitute(q = q, p = p, output = fq_out)
        debug(commandLine)
        zcatProcess = Popen('zcat %s' % fq_in, shell = True, stdout = subprocess.PIPE)
        fqfProcess = Popen(commandLine, shell = True, stdin = zcatProcess.stdout)
    else:
        commandLine = fqfFileTemplate.substitute(q = q, p = p, output = fq_out, input = fq_in)
        debug(commandLine)
        fqfProcess = Popen(commandLine,
                           shell = True) 
    fqfProcess.wait() 

def iterative_DBR_dict(in_dir, seqType, read, save_path):
    if seqType == 'pear':
        # read the 8 bases at the end
        dbr_start = -9
        dbr_stop = None
    else:
        # read the 8 bases at the beginning
        dbr_start = 0
        dbr_stop = 9
    files = os.listdir(in_dir)
    for f in files:
        toRead = re.findall(read, f)
        if toRead: 
            DBR_dict(f, dbr_start, dbr_stop, test_dict=True, save_path=save_path)


def Demultiplex(in_file, barcode_file, out_dir, out_prefix):    
    print 'Demultiplexing sequence data with FASTX Toolkit.\n'
    
    prefix_path = out_dir + out_prefix
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    demultiplex_call = 'zcat ' + in_file + ' | /usr/local/bin/fastx_barcode_splitter.pl --bcfile ' + barcode_file + ' --prefix ' + prefix_path + ' --bol'
    
    subprocess.call(demultiplex_call, shell=True)
    return

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
    return

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
    return

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
    return

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
    return

''' Deprecated
def FASTQ_R1_R2_merge(in_dir, fq_r1, fq_r2, fq_out):
    print 'Taking reverse complement of read 2 with FASTX Toolkit.\n'
    rc_call = 'fastx_reverse_complement -i ' + in_dir + fq_r2 + ' -o ' + in_dir + fq_out
    subprocess.call(rc_call, shell=True)
    
    print 'Merging read 1 and read 2 reverse complement into a single FASTQ file.'
    merge_call = 'paste ' + in_dir + fq_r1 + ' ' + in_dir + fq_out + ' > ' + in_dir + fq_out + '_merged.fq'
    subprocess.call(merge_call, shell=True)
    rm_tab_call = "sed 's/[ \t]//g' " + in_dir + fq_out + '_merged.fq' + ' >> ' + in_dir + fq_out + '_untabbed_merged.fq'
    subprocess.call(rm_tab_call, shell=True)
    return
'''
