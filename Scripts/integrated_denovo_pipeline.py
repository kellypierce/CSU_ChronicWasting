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
import string as string
from string import Template, join
import logging as logging
from logging import debug, critical, error, info
import sys as sys
import warnings
import fnmatch
import json
import itertools
import gzip
import multiprocessing as mp
import pdb

################################## GLOBALS ####################################

# paths to executables on cluster (all paths are in ~./bashrc)
pearPath = 'pear-0.9.6-bin-64'
qualityFilter = 'fastq_quality_filter'
trimmer = 'fastx_trimmer'
demultiplexer = 'fastx_barcode_splitter.pl'
denovo_path = 'denovo_map.pl'
stacks_executables = '/home/pierce/bin/stacks-1.35/' # this is where ustacks lives, which should be all we need. other stacks scripts are in /opt/software/stacks-1.26/scripts
BWA = 'bwa'
samtoolsPath = 'samtools'
bcftoolsPath = 'bcftools'

# a note on Stacks:
# The standard installation does not configure the installation directory as expected. Consequently, denovo_map.pl can't find the other executable files.
# Manually editing denovo_map.pl does not seem to help (I'm sure there's a way, but picking apart the code is not something I'm anxious to do).
# The work-around is this:
# - install Stacks in home directory (not system-wide, because subsequent file path edits may cause problems for others)
# - move all the contents of /your/install/location/stacks-X.XX/scripts up one level so they are in .../stacks-X.XX (this is necessary because ustacks et al. are not in scripts, but we really want all those executables in the same place)
# - ensure that all those scripts have execution privileges
# - run denovo_map.pl with the "-e" command that specifies the path to executables as /your/install/location/stacks-X.XX

# PEAR assembly
pearSimpleTemplate = Template('%s -f $f -r $r -o $o' % pearPath)
pearExtraParamsTemplate = Template('%s -f $f -r $r -o $o $e' % pearPath)

# Fastq filtering of assembled data
fqfStdinTemplate = Template('%s -q $q -p $p -Q33 -z -o $output' % qualityFilter)
fqfFileTemplate = Template('%s -q $q -p $p -Q33 -z -i $input -o $output' % qualityFilter)

# Demultiplexing
# fastx_barcode_splitter always takes input as stdin
# files across all libraries should go to the same directory
demultiplexStdinTemplate = Template('%s --bcfile $b --prefix $p --bol' % demultiplexer)

# Trimming
uniformLengthTemplate = Template('%s -f $f -l $l -i $in_path -o $out_path' % trimmer)

# Stacks de novo assembly
# Extract reference from Stacks consensus, index
# Reference-based assembly
# DBR filter

# Genotype calling
samtoolsView = Template('%s view -F 4 -b -S -o $output $input' % samtoolsPath)
samtoolsSort = Template('%s sort -o $output $input' % samtoolsPath)
samtoolsIndex = Template('%s index $input' % samtoolsPath)
samtoolsMpileup = Template('%s mpileup -t DP -C50 -u -I -f $reference -o $bcf_out $input' % samtoolsPath)
bcftoolsView = Template('%s call -v -m $input > $output' % bcftoolsPath)
    

###############################################################################

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

# run configureLogging first
configureLogging(False)

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
        
'''
def checkDir(dirname):
    #return true if this is a directory and is readable on the current filesystem
    try:
        if os.path.exists(os.path.abspath(dirname)) and os.path.isdir(os.path.abspath(dirname)) and os.access(os.path.abspath(dirname), R_OK):
            return True
    except IOError:
        return False
'''

def checkExe(filename):
    '''
    return true if this is an executable file on the current filesystem
    '''
    if not isinstance(filename, str):
        raise TypeError("need a string, got a %s" % type(filename))
    return (os.path.exists(filename) and os.path.isfile(filename) and os.access(filename, X_OK))

def iterative_DBR_dict(in_dir, seqType, save, dbr_start, dbr_stop):
    #if not checkDir(in_dir):
    #    raise IOError("Input is not a directory: %s" % in_dir)
    if seqType == 'read2':
        warnings.warn('Expect directory containing only Read 2 files; any other files present in %s will be incorporated into DBR dictionary.' % in_dir)
        # read the 8 bases at the beginning
        #dbr_start = 0
        #dbr_stop = 9
    elif seqType == 'pear':
        files = os.listdir(in_dir)
        for f in files:
	    inFile = in_dir + f
            DBR_dict(inFile, dbr_start, dbr_stop, test_dict=True, save=save, f=f)
    else:
        raise IOError("Input sequence type specified as %s. Options are 'pear' or 'read2'." % seqType)

def DBR_dict(fq_path, dbr_start, dbr_stop, test_dict = False, save = None, f=None):
    # DBR is in read 2
    # if merged, it will be the last -2 to -9 (inclusive) bases, starting with base 0 and counting from the end
    # if not merged, it will be bases 2 to 9
    if not checkFile(fq_path):
        raise IOError("where is the input file: %s" % fq_path)
    info('Creating {ID: dbr} dictionary from %s.' % fq_path)
    dbr = {}
    fq_line = 1
    if fq_path.endswith('gz'):
        openFxn = gzip.open
    else:
        openFxn = open
    with openFxn(fq_path, 'r') as db:
        for line in db:
            if fq_line == 1:
                ID = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', line)[1]
                fq_line = 2
            elif fq_line == 2:
                seq = list(line) # split the sequence line into a list
                tag = ''.join(seq[dbr_start:dbr_stop])
                dbr[ID] = tag
                fq_line = 3
            elif fq_line == 3:
                fq_line = 4
            elif fq_line == 4:
                fq_line = 1
    if test_dict:
        print 'Checking DBR dictionary format.'
        x = itertools.islice(dbr.iteritems(), 0, 4)
        for key, value in x:
            print key, value
        #print dbr['8:1101:15808:1492'] # this is the first entry in /home/antolinlab/Downloads/CWD_RADseq/pear_merged_Library12_L8.assembled.fastq
    if save:
        if not os.path.exists(save):
		os.makedirs(save)
	fq_name = os.path.splitext(f)[0]
	fq_dbr_out = save + fq_name + '.json'
        print 'Writing dictionary to ' + fq_dbr_out
        with open(fq_dbr_out, 'w') as fp:          
            json.dump(dbr, fp)
            
def iterative_PEAR_assemble(in_dir, out_dir, out_name, extra_params, regexR1='*', regexR2='*'):
    files = os.listdir(in_dir)
    #print(in_dir, files)
    read1 = fnmatch.filter(files, '*'+regexR1+'*')
    mergedProcess = [mp.Process(target=PEAR_assemble, args=(in_dir, 
                                                    r1, 
                                                    re.sub(regexR1, regexR2, r1), 
                                                    out_dir, 
                                                    out_name, 
                                                    extra_params)) for r1 in read1]
    for mP in mergedProcess:
        mP.start()
    for mP in mergedProcess:
        mP.join()
    
def PEAR_assemble(in_dir, forward, reverse, out_dir, out_name,  extra_params=None):
    if not checkFile(in_dir + forward):
        raise IOError("Where is the forward read file: %s" % forward)
    if not checkFile(in_dir + reverse):
        raise IOError("where is the reverse read file: %s" % reverse)
    
    if not forward.endswith(".fastq.gz"):
        warnings.warn("Expect raw sequence data in .fastq.gz format")
    if not reverse.endswith(".fastq.gz"):
        warnings.warn("Expect raw sequence data in .fastq.gz format")
        
    #if not checkExe(pearPath):
    #    raise Exception("could not find %s in the filesystem for execution, is the environment setup correctly?" % pearPath)
    info('Merging overlapping reads %s and %s with PEAR.' % (forward, reverse))
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    out_split = os.path.splitext(forward)
    out_suffix = re.sub(r'R\d{1}', '', out_split[0])
    out_full = out_dir + out_name + out_suffix
    
    info('Identified read 1 as %s' % forward)
    info('Identified read 2 as %s' % reverse)
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
    subprocess.call(commandLine, shell=True)
    #pearProcess = Popen(commandLine, shell=True)
    #pearProcess.wait()

def iterative_FASTQ_quality_filter(directory, out_dir, out_name, q, p, read='*'):
    files = os.listdir(directory)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    filterProcess = []
    for f in files:
        r2 = re.findall(read, f)
        if r2: 
            fileRoot = os.path.splitext(f)
            # untested: file path for outputs. they went to the wrong place in part3, so I fixed what I thought was the bug in the file path, but just moved the data to the expected directory manually isntead of re-running
            out_file = out_dir + fileRoot[0] + out_name
            in_file = directory + f
            info('Quality filtering %s' % f)
            info('Results saved to %s' % out_dir)
	    filterProcess.append(mp.Process(target=FASTQ_quality_filter, args=(in_file, out_file, q, p)))
    for fP in filterProcess:
        fP.start()
    for fP in filterProcess:
        fP.join()

def FASTQ_quality_filter(fq_in, fq_out, q, p, qualityFilter = qualityFilter):
    if not checkFile(fq_in):
        raise IOError("where is the input file: %s" % fq_in)
    if not fq_in.endswith("gz"):
        warnings.warn("prefer to pass compressed files, consider gzipping the inputs")
    #if not checkExe(qualityFilter):
    #    raise Exception("could not find %s in the filesystem for execution, is the environment setup correctly?" % qualityFilter)
    info("Quality filtering with FASTX Toolkit. Output file %s" % fq_out)
    
    # all the files will be saved as gzipped, so they need a .gz extension... but people (me) will probably forget to put gz in the file name
    if 'gz' in fq_out:
        out = fq_out
    else:
        out = fq_out + '.gz'
    
    if fq_in.endswith('gz'): # chain a gz decompressor thread to fqf
        commandLine = fqfStdinTemplate.substitute(q = q, p = p, output = out)
        debug(commandLine)
        zcatProcess = Popen('zcat %s' % fq_in, shell = True, stdout = subprocess.PIPE)
        fqfProcess = Popen(commandLine, shell = True, stdin = zcatProcess.stdout)
    else:
        commandLine = fqfFileTemplate.substitute(q = q, p = p, output = out, input = fq_in)
        debug(commandLine)
        fqfProcess = Popen(commandLine,
                           shell = True) 
    fqfProcess.wait() 

## TRIM R2 END OF MERGED SEQUENCE BEFORE DEMULTIPLEXING TO ENFORCE UNIFORM READ LENGTH?

def Trim(in_dir, out_dir, suffix, first_base, last_base=None):    
    info('Trimming DBR and enzyme cut sites from %s with FASTX Toolkit.' % in_dir)
    
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
        if last_base:
            trim_call = "fastx_trimmer -f " + str(first_base) + ' -l ' + str(last_base) + " -i " + full_path + " -o " + new_path
            subprocess.call(trim_call, shell=True)
        else:
            trim_call = "fastx_trimmer -Q33 -f " + str(first_base) + " -i " + full_path + " -o " + new_path
            subprocess.call(trim_call, shell=True)
    return

def iterative_Demultiplex(in_dir, # directory of un-demultiplexed libraries
                          barcode_dir, #directory containing the barcodes for each library
                          out_dir, # full path for outputs 
                          out_prefix): # text string to add to file names
    #if not checkDir(in_dir):
    #    raise IOError("Input is not a directory: %s" % in_dir)
    #if not checkFile(barcode_file):
    #    raise IOError("Where is the barcode file? %s" % barcode_file)
    #pdb.set_trace()
    files = os.listdir(in_dir) 
    for f in files:
        sampleID_match = re.match(".*(Library\d{2,3}).*", f)
        #sampleID = re.match(".*(\d{3}[a-z]?).*", f).groups()[0]
        if sampleID_match: # if we get a match
            sampleID = sampleID_match.groups()[0] # extract that match
            bcs = os.listdir(barcode_dir)
            out_name = out_prefix + sampleID + '_'
            for b in bcs:
                if sampleID in b:
                    barcode_file = barcode_dir + '/' + b
                    in_f = in_dir + '/' + f
                    Demultiplex(in_f, barcode_file, out_dir, out_name)
                

def Demultiplex(in_file, barcode_file, out_dir, out_prefix): 
#	if not checkFile(in_file):
#        raise IOError("Input is not a file: %s" % in_file)
#	if not checkFile(barcode_file):
#        raise IOError("Where is the barcode file? %s" % barcode_file)   
    
    prefix_path = out_dir + '/' + out_prefix

    info('Demultiplexing %s with FASTX Toolkit; output files saved in %s' % (in_file, prefix_path))
        
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    commandLine = demultiplexStdinTemplate.substitute(b = barcode_file, p = prefix_path)
        
    if in_file.endswith('gz'): # chain a gz decompressor thread to fqf
        zcatProcess = Popen('zcat %s' % in_file, shell = True, stdout = subprocess.PIPE)
        demultiplexProcess = Popen(commandLine, shell = True, stdin = zcatProcess.stdout)
    else:
        catProcess = Popen('cat %s' % in_file, shell = True, stdout = subprocess.PIPE)
        demultiplexProcess = Popen(commandLine, shell = True, stdin = catProcess.stdout)
    demultiplexProcess.wait()

def denovo_Stacks(in_dir, denovo_path, stacks_executables, out_dir, m, n, b, D):    
    print 'Assembling sequences de novo using Stacks denovo_map.pl\n'
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    rm_unmatched = False
    
    ustacks_path = stacks_executables + '/ustacks'
    print ustacks_path
    ustacksMessageTemplate = Template('*****Processing sample $n of $total with USTACKS.*****')
    ntotal = len(os.listdir(in_dir))
    c = 0
    for i in os.listdir(in_dir):
        
        print ustacksMessageTemplate.substitute(n = c, total = ntotal)
        c += 1
        
        # don't proceed if the inputs don't have the proper extension
        assert '.fq' in i, "%s needs extension .fq for Stacks compatibility" % i

        if 'unmatched' in i: # don't process the unmatched sequences
            print 'Skipping sequences with unmatched barcodes.\n'
            rm_unmatched = True
        else:
            
            # full path to trimmed file
            new_path=os.path.join(in_dir, i)            
            # Run ustacks
            # example usage: ustacks -t fastq -f ./samples/f0_male.fq    -o ./stacks -i 1 -d -r -m 3 -p 15
            ustacks_args = [ustacks_path, ' -t fastq ' + ' -f ' + new_path + ' -o ', out_dir, ' -m ', str(m), ' -r -R']
            ustacks_call = ''.join(ustacks_args)
            subprocess.call(ustacks_call, shell=True)
    
    cstacksMessageTemplate = Template('*****Preparing sample $name for CSTACKS.*****')
    segfaultWarningTemplate = Template('### WARNING ### empty .alleles.tsv file for %sample')
    # Generate "-s" arguments to Stacks denovo_map.pl
    s_list=[]
    
    for j in os.listdir(out_dir):
        # a clumsy hack to deal with empty allele.tsv files:
        # the alleles.tsv file is the shortest one; if it only has one line (a header) it will cause a segfault
        if 'alleles' in j:
            alleles_file = os.path.join(out_dir, j)
            with open(alleles_file) as f:
                lineCount = sum(1 for _ in f)
            if lineCount > 1:
                tagsFile = re.sub('alleles', 'tags', alleles_file)
                print cstacksMessageTemplate.substitute(name = tagsFile)
                basej = os.path.splitext(os.path.splitext(tagsFile)[0])[0]
                #tag_path = os.path.join(out_dir, basej)
                # generate a string for the -s cstacks argument
                s_list.append('-s ')
                #s_list.append(tag_path)
                s_list.append(basej)
                s_list.append(' ')
            else:
                basej = os.path.splitext(os.path.splitext(alleles_file)[0])[0]
                print segfaultWarningTemplate.substitute(sample = basej)
        # TODO: find some way to warn when samples are dropped because they don't get tags files after ustacks runs!
        #print j
        #if 'tags' in j: # we only want the tags.tsv files for cstacks
            
            
        
    # join list into a string to pass to denovo_map.pl
    formatted_list = ''.join(s_list)

    # Run denovo_map.pl
    #build_args = [denovo_path, ' -e ', stacks_executables, ' -o ', out_dir, ' -m ', str(m), ' -n ', str(n), ' -t ', '-b ', str(b), ' -D ', D, ' -S ', formatted_list]
    #denovo_call = ''.join(build_args)
    #subprocess.call(denovo_call, shell=True)
    
    # Run cstacks
    # example usage: cstacks -b 1 -o ./stacks -s ./stacks/f0_male -s ./stacks/f0_female -p 15
    cstacks_path = stacks_executables + '/cstacks'
    cstacks_args = [cstacks_path, ' -b 1 -n ', str(n), ' -o ', out_dir, ' ', formatted_list]
    cstacks_call = ''.join(cstacks_args)
    print cstacks_call
    subprocess.call(cstacks_call, shell=True)
    
    ### Add a line to move the files!
    
    return

def GeneratePseudoref(in_dir, out_file, BWA_path):    
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

    #out_path=os.path.join(out_dir, out_name)
    print "Writing pseudoreference genome to file " + out_file + '\n'
    #with open(out_path, 'w') as fq: # need to check if it already exists
    with open(out_file, 'w') as fq:
        fq.write(fastq_interleaved)
    
    # Index the new pseudoreference genome
    print 'Indexing the reference genome using BWA.\n'
    index_call =  BWA_path + ' index ' + out_file # index = BWA function to use
    print index_call
    subprocess.call(index_call, shell = True)
    return

def refmap_BWA(in_dir, out_dir, BWA_path, pseudoref_full_path):    
    
    #### NEED TO CHECK IF LIBRARIES SPLIT ACROSS LANES (AND IN DIFFERENT FASTQ FILES) HAVE SAMPLES OVERWRITTEN HERE
    #### I SUSPECT THIS IS THE CASE; IF SO FASTQ FILES SHOULD BE CONSOLIDATED BY LIBRARY (JUST CAT THE FILES)
    
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

def callGeno(sam_in, pseudoref, BCFout, VCFout):
    #print sam_in, pseudoref, BCFout, VCFout
    # set up the individual files for transfer from sam to bam and bam indexing
    
    print 'Processing sam files into sorted bam files.'
    
    for sam in os.listdir(sam_in):
        # TODO: add if 'sam' in sam (check that we don't already have bam files)
        # samtools view -F 4 will filter OUT reads with bitwise flag 0004 -- these are unmapped reads
        fname = os.path.splitext(sam)[0]
        
        samPath = sam_in + '/' + sam
        bam = sam_in + '/' + fname + '.bam' # for bam output
        sorted = sam_in + '/' + fname + '.sorted.bam' # for sorting output
        
        view_cmd = samtoolsView.substitute(output = bam, input = samPath)
        subprocess.call(view_cmd, shell=True)
        
        sort_cmd = samtoolsSort.substitute(output = sorted, input = bam)
        subprocess.call(sort_cmd, shell=True)
        
        index_cmd = samtoolsIndex.substitute(input = sorted)
        subprocess.call(index_cmd, shell=True)
    
    # take the sorted, indexed bam files and perform the genotype calling with mpileup
    print 'Calling genotypes with samtools mpileup'
    wildcard_in = sam_in + '/*.sorted.bam'
    #print wildcard_in
    mpileup_cmd = samtoolsMpileup.substitute(reference = pseudoref, 
                                             input = wildcard_in, 
                                             bcf_out = BCFout)
                                             
    #print mpileup_cmd
    subprocess.call(mpileup_cmd, shell=True)
    
    # convert the resulting bcf file to a vcf file
    print 'Converting genotypes file to VCF format'
    bcfView_cmd = bcftoolsView.substitute(input = BCFout,
                                          output = VCFout)
    #print bcfView_cmd
    subprocess.call(bcfView_cmd, shell=True)
    print 'RUN COMPLETED.'
        
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
