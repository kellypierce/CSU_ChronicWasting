#!/usr/bin/python

#########################################################################################
#                                                                                       #
#    integrated_denovo_pipeline.py: a denovo pipeline for processing DBR RADseq Data    #
#    Copyright (C) 2016 Kelly Anne Pierce (kellyannepierce@gmail.com)                   #
#                                                                                       #
#########################################################################################

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
from collections import defaultdict
from collections import Counter
import time
import pdb
import multiprocessing
from Queue import Queue
from threading import Thread


################################## GLOBALS ####################################

# Process queue
processQueue = Queue()

# paths to executables on cluster (all paths are in ~./bashrc)
pearPath = 'pear-0.9.6-bin-64'
qualityFilter = 'fastq_quality_filter'
trimmer = 'fastx_trimmer'
demultiplexPath = 'fastx_barcode_splitter.pl'
denovo_path = 'denovo_map.pl'
stacks_executables = '/home/pierce/bin/stacks-1.35/' # this is where ustacks lives, which should be all we need. other stacks scripts are in /opt/software/stacks-1.26/scripts
BWA = 'bwa'
samtoolsPath = 'samtools'
bcftoolsPath = 'bcftools'
phred_dict = {'"':1.0,"#":2.0,"$":3.0,"%":4.0,"&":5.0,"'":6.0,"(":7.0,")":8.0,"*":9.0,"+":10.0,
              ",":11.0,"-":12.0,".":13.0,"/":14.0,"0":15.0,"1":16,"2":17.0,"3":18.0,"4":19.0,"5":20.0,
              "6":21.0,"7":22.0,"8":23.0,"9":24.0,":":25.0,";":26,"<":27.0,"+":28.0,">":29.0,"?":30.0,
              "@":31.0,"A":32.0,"B":33.0,"C":34.0,"D":35.0,"E":36,"F":37.0,"G":38.0,"H":39.0,"I":40.0,
              "J":41.0,"K":42.0}
# conversion reference: http://drive5.com/usearch/manual/quality_score.html

# a note on Stacks:
# The standard installation does not configure the installation directory as expected. Consequently, denovo_map.pl can't find the other executable files.
# Manually editing denovo_map.pl does not seem to help (I'm sure there's a way, but picking apart the code is not something I'm anxious to do).
# The work-around is this:
# - install Stacks in home directory (not system-wide, because subsequent file path edits may cause problems for others)
# - move all the contents of /your/install/location/stacks-X.XX/scripts up one level so they are in .../stacks-X.XX (this is necessary because ustacks et al. are not in scripts, but we really want all those executables in the same place)
# - ensure that all those scripts have execution privileges
# - run denovo_map.pl with the "-e" command that specifies the path to executables as /your/install/location/stacks-X.XX

###############################################################################

# To do
# 1. Check that SAM files contain a map for all the sequences so that FASTQ filtering doesn't leave some bad quality data behind
#    Alternatively, instead of tracking tags to remove, track tags to keep
# 2. Make this work for data that span multiple libraries by calling the function separately for each library following the template below:
#    DBR_filter(assembled_dir = '/path/to/assembly_library1',
#               dict_in = '/path/to/dbr_dict_library1',
#               out_seqs = '/path/to/filtered_library1.fastq',
#               n_expected = 2)

## Parallelization things from Joe
class Work():
    def __init__(self, commandline, shell, cwd = os.getcwd(), libraryPath = None, analysisFile = None):
        self.commandline = commandline
        # probably do not need this libraryPath code... for Joe's dynamically loaded libraries with duplicate names
        if libraryPath is not None:
            self.env = os.environ.copy()
            self.env["LD_LIBRARY_PATH"] = "%s:%s" % (libraryPath, self.env["LD_LIBRARY_PATH"])
        else:
            self.env = os.environ.copy()
            self.shell = shell
            self.cwd = cwd

def worker():
    'run subprocesses in an orderly fashion'
    while True:
        workItem = processQueue.get()
        if workItem.commandline != None:
            p = Popen(workItem.commandline,
                      env = workItem.env,
                      shell = workItem.shell,
                      cwd = workItem.cwd)
            p.wait()
            processQueue.task_done()

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

def qual_median(QUAL, phred_dict):
    
    listQUAL = list(QUAL)
    list_intQUAL =[]
    for q in listQUAL:
        list_intQUAL.append(phred_dict[q])
    #return np.median(list_intQUAL)
    #pdb.set_trace()
    list_intQUAL.sort() # this modifies the list in place -- you can't assign it to a new variable. alt, qsort = sorted(list_intQUAL) if we need the original list preserved for some reason
    qlen = len(list_intQUAL)
    if qlen % 2 == 0: # even length list -- take the average of the two middle values
        median_qual = (list_intQUAL[(qlen/2)-1]+list_intQUAL[(qlen/2)])/2
    else: # odd length list -- take the middle value (list length / 2 will be automatically rounded up by python integer division)
        median_qual = list_intQUAL[qlen/2]
    return median_qual
    
def find_SampleID(filename, regexSample):
    #sampleID_match = re.match(".*(\d{3}[a-z]?).*", filename)
    # this revision is VERY specific to my technical replicates
    # TODO: find a way to pass the regex capture group as an argument so that this (and related functions) are more flexible
    #sampleID_match = re.match(".*(\d{1,3}T?).*", filename)
    sampleID_match = re.match(".*("+regexSample+").*", filename)
    if sampleID_match:
        sampleID = sampleID_match.groups()[0]
        return sampleID
    else:
        return None
    
#def find_LibraryID(filename):
#    #libraryID_match = re.match(".*(Library\d{2,3}).*", filename)
#    libraryID_match = re.match(".*(Library\d{1,3}[A|B]?).*", filename)
#    if libraryID_match: # if we get a match (this allows the script to proceed if a file has a mismatched name)
#        libraryID = libraryID_match.groups()[0] # extract the library ID match
#        return libraryID 
#    else:
#        return None
    
    
def find_LibraryID(filename, regexLibrary):
    #libraryID_match = re.match(".*(Library\d{2,3}).*", filename)
    #libraryID_match = re.match(".*(Library\d{1,3}[A|B]?).*", filename)
    libraryID_match = re.match('.*('+regexLibrary+').*')
    if libraryID_match: # if we get a match (this allows the script to proceed if a file has a mismatched name)
        libraryID = libraryID_match.groups()[0] # extract the library ID match
        return libraryID 
    else:
        return None

def find_BarcodeFile(library, directory):
    if library:
        if os.path.isdir(directory):
            bcs = os.listdir(directory)
            for b in bcs:
                if library in b:
                    bcf = directory + '/' + b
                    return bcf
        else: # if it's just a single file
            if os.path.isfile(directory):
                return directory
            else:
                return None
    else:
        return None

def parallel_PEAR_assemble(regexR1, regexR2, regexLibrary, in_dir, out_dir, pearPath, out_name = 'pear_merged_', extra_params=None):
    files = os.listdir(in_dir)
    
    #print(in_dir, files)
    # find Read 1 file
    read1 = fnmatch.filter(files, '*'+regexR1+'*')
    # which library?
    read1Library = find_LibraryID(read1, regexLibrary)
    # file name for output
    out_file = out_name + read1Library
    
    # assemble parallel processes (calls to function PEAR_assemble)
    mergedProcess = [mp.Process(target=PEAR_assemble, args=(in_dir + r1,
                                                            in_dir + re.sub(regexR1, regexR2, r1),
                                                            out_dir, 
                                                            out_file, 
                                                            extra_params,
                                                            pearPath)) for r1 in read1]
    for mP in mergedProcess:
        mP.start()
    for mP in mergedProcess:
        mP.join()
    
def PEAR_assemble(R1, R2, out_dir, out_file, pearPath, extra_params=None):
    if not checkFile(R1):
        raise IOError("Where is the Read 1 file: %s" % R1)
    if not checkFile(R2):
        raise IOError("Where is the Read 2 read file: %s" % R2)
    
    if not R1.endswith(".fastq.gz"):
        warnings.warn("Expect raw sequence data in .fastq.gz format")
    if not R2.endswith(".fastq.gz"):
        warnings.warn("Expect raw sequence data in .fastq.gz format")
        
    #if not checkExe(pearPath):
    #    raise Exception("could not find %s in the filesystem for execution, is the environment setup correctly?" % pearPath)
    info('Merging overlapping reads %s and %s with PEAR.' % (R1, R2))
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    info('Identified read 1 as %s' % R1)
    info('Identified read 2 as %s' % R2)
    info('Saving outputs to %s' % out_file)
    
    pearSimpleTemplate = Template('%s -f $f -r $r -o $o' % pearPath)
    pearExtraParamsTemplate = Template('%s -f $f -r $r -o $o $e' % pearPath)
    
    if extra_params:
        commandLine = pearExtraParamsTemplate.substitute(f = R1,
                                                         r = R2,
                                                         o = out_file,
                                                         e = extra_params)
    else:
        commandLine = pearSimpleTemplate.substitute(f = R1,
                                                    r = R2,
                                                    o = out_file)
    subprocess.call(commandLine, shell=True)
    #pearProcess = Popen(commandLine, shell=True)
    #pearProcess.wait()

def parallel_FASTQ_quality_filter(in_dir, out_dir, out_name, q, p, qualityFilter, read='*'):
    
    # find all the files in the input directory
    files = os.listdir(in_dir)
    
    # make the output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    # container for parallel processes
    filterProcess = []
    
    # iterate through the input files
    for f in files:
        # don't assume that only merged reads are in directory...
        r2 = re.findall(read, f)
        if r2: 
            # get the base name of the file (without any extension)
            fileRoot = os.path.splitext(f)[0]
            # untested: file path for outputs. they went to the wrong place in part3, so I fixed what I thought was the bug in the file path, but just moved the data to the expected directory manually isntead of re-running
            out_file = out_dir + fileRoot + out_name
            in_file = in_dir + f
            # add the process to the list
            filterProcess.append(mp.Process(target=FASTQ_quality_filter, args=(in_file, out_file, q, p, qualityFilter)))
            
    for fP in filterProcess:
        fP.start()
    for fP in filterProcess:
        fP.join()

def FASTQ_quality_filter(in_file, out_file, q, p, qualityFilter):
    if not checkFile(in_file):
        raise IOError("where is the input file: %s" % in_file)
    if not in_file.endswith("gz"):
        warnings.warn("prefer to pass compressed files, consider gzipping the inputs")
    #if not checkExe(qualityFilter):
    #    raise Exception("could not find %s in the filesystem for execution, is the environment setup correctly?" % qualityFilter)
    
    info("Quality filtering with FASTX Toolkit. Output file %s" % out_file)
    
    # all the files will be saved as gzipped, so they need a .gz extension... but people (me) will probably forget to put gz in the file name
    if 'gz' in out_file:
        out = out_file
    else:
        out = out_file + '.gz'
    
    info('Quality filtering %s' % in_file)
    info('Results saved to %s' % out_file)
    
    fqfStdinTemplate = Template('%s -q $q -p $p -Q33 -z -o $output' % qualityFilter)
    fqfFileTemplate = Template('%s -q $q -p $p -Q33 -z -i $input -o $output' % qualityFilter)
    
    if in_file.endswith('gz'): # chain a gz decompressor thread to fqf
        commandLine = fqfStdinTemplate.substitute(q = q, p = p, output = out)
        debug(commandLine)
        zcatProcess = Popen('zcat %s' % in_file, shell = True, stdout = subprocess.PIPE)
        fqfProcess = Popen(commandLine, shell = True, stdin = zcatProcess.stdout)
    else:
        commandLine = fqfFileTemplate.substitute(q = q, p = p, output = out, input = in_file)
        debug(commandLine)
        fqfProcess = Popen(commandLine, shell = True) 

## TRIM R2 END OF MERGED SEQUENCE BEFORE DEMULTIPLEXING TO ENFORCE UNIFORM READ LENGTH?

def parallel_Trim(in_dir, out_dir, trimPath, first_base, last_base=None, suffix = '_trimmed.fq'):
    
    # new directory for trimmed files
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    trimProcess = []
        
    for i in os.listdir(in_dir):
        # save file as out_file
        out_name = i + suffix
        
        # full path to file for trimming
        in_file=os.path.join(in_dir, i)
        out_file=os.path.join(out_dir, out_name)
        
        trimProcess.append(mp.Process(target=Trim, args=(in_file, out_file, first_base, last_base, trimPath)))
    
    for tP in trimProcess:
        tP.start()
    for tP in trimProcess:
        tP.join()

def Trim(in_file, out_file, first_base, trimPath, last_base=None):    
    
    info('Trimming DBR and enzyme cut sites from %s with FASTX Toolkit.' % in_file)
    
    uniformLengthTemplate = Template('%s -f $f -l $l -i $input -o $output' % trimPath)
    uniformLengthTemplate_Q33 = Template('%s -Q33 -f $f -i $input -o $output' % trimPath)

    # Remove barcodes and R1 enzyme cut site
    if last_base:
        trim_call = uniformLengthTemplate.substitute(f = str(first_base), l = str(last_base), input = in_file, output = out_file)
        #trim_call = "fastx_trimmer -f " + str(first_base) + ' -l ' + str(last_base) + " -i " + full_path + " -o " + new_path
        subprocess.call(trim_call, shell=True)
    else:
        trim_call = uniformLengthTemplate_Q33.substitute(f = str(first_base), input = in_file, output = out_file) #can't recall why this one needs the Q33 arg and the other doesn't... something about whether or not FASTQ Trimmer recognizes the qual encoding of pear merged vs. raw illumina data?
        #trim_call = "fastx_trimmer -Q33 -f " + str(first_base) + " -i " + full_path + " -o " + new_path
        subprocess.call(trim_call, shell=True)
    return

def iterative_Demultiplex(in_dir, # directory of un-demultiplexed libraries
                          barcode_dir, #directory containing the barcodes for each library
                          out_dir, # full path for outputs 
                          regexLibrary,
                          out_prefix = 'demultiplexed_'): # text string to add to file names

    #if not checkDir(in_dir):
    #    raise IOError("Input is not a directory: %s" % in_dir)
    #if not checkFile(barcode_file):
    #    raise IOError("Where is the barcode file? %s" % barcode_file)
    #pdb.set_trace()
    
    files = os.listdir(in_dir) 
    
    demultiplexProcess = []
    
    for f in files:
        # for libraries with format 'Library#' in name (only numbers to distinguish library)
        #sampleID_match = re.match(".*(Library\d{2,3}).*", f)
        
        # for libraries with formats 'Library#' or 'Library#A' in name (letters and numbers to distinguish library)
        # this should also work for libraries with only numbers: '\w?' should capture 0 or more words after the digits
        
        #sampleID_match = re.match(".*(Library\d{1,3}\w?).*", f)
        #sampleID_match = re.match(".*(Library\d{1,3}[A|B]?).*", f)
        sampleID_match = re.match(".*("+regexLibrary+").*", f)
        
        if sampleID_match: # if we get a match
            sampleID = sampleID_match.groups()[0] # extract that match
            bcs = os.listdir(barcode_dir)
            out_name = out_prefix + sampleID + '_'
            for b in bcs:
                if sampleID in b:
                    barcode_file = barcode_dir + '/' + b
                    in_f = in_dir + '/' + f
                    #Demultiplex(in_f, barcode_file, out_dir, out_name)
                    
                    demultiplexProcess.append(mp.Process(target=Trim, args=(in_f, barcode_file, out_dir, out_name)))
    
    for dP in demultiplexProcess:
        dP.start()
    for dP in demultiplexProcess:
        dP.join()
                

def Demultiplex(in_file, barcode_file, out_dir, demultiplexPath, out_prefix = 'demultiplexed_'): 
#	if not checkFile(in_file):
#        raise IOError("Input is not a file: %s" % in_file)
#	if not checkFile(barcode_file):
#        raise IOError("Where is the barcode file? %s" % barcode_file)   
    
    prefix_path = out_dir + '/' + out_prefix

    info('Demultiplexing %s with FASTX Toolkit; output files saved in %s' % (in_file, prefix_path))
        
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # fastx_barcode_splitter always takes input as stdin
    # files across all libraries should go to the same directory
    demultiplexStdinTemplate = Template('%s --bcfile $b --prefix $p --bol' % demultiplexPath)

    commandLine = demultiplexStdinTemplate.substitute(b = barcode_file, p = prefix_path)
        
    if in_file.endswith('gz'): # chain a gz decompressor thread to fqf
        zcatProcess = Popen('zcat %s' % in_file, shell = True, stdout = subprocess.PIPE)
        demultiplexProcess = Popen(commandLine, shell = True, stdin = zcatProcess.stdout)
    else:
        catProcess = Popen('cat %s' % in_file, shell = True, stdout = subprocess.PIPE)
        demultiplexProcess = Popen(commandLine, shell = True, stdin = catProcess.stdout)
    demultiplexProcess.wait()

def denovo_Ustacks(in_dir, denovo_path, stacks_executables, out_dir, m, n, b, D):    
    print 'Assembling sequences de novo using ustacks\n'
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    #rm_unmatched = False
    
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

def denovo_Cstacks(in_dir, denovo_path, stacks_executables, out_dir, m, n, b, D):    
    
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

def parallel_refmap_BWA(in_dir, out_dir, BWA_path, pseudoref_full_path):

    print 'Mapping sequence data to pseudoreference genome using BWA.\n'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)    
    
    #refmapProcess = []
    
    # the regex below also finds unmatched samples    
    rex = re.compile(r'\d+')
    
    for i in os.listdir(in_dir):
        if 'unmatched' not in i: # independently check to remove unmatched files from the files to be refmapped
            if rex.search(i):
                fname, fext = os.path.splitext(i)
                in_file = in_dir + i
                out_file = out_dir + fname + '.sam'
                commandline = refmap_BWA(in_file, fname, out_file, BWA_path, pseudoref_full_path, execute=False)
                #refmapProcess.append(mp.Process(target=refmap_BWA, args=(in_file, fname, out_file, BWA_path, pseudoref_full_path)))
        
                processQueue.put(Work(commandline = commandline, shell = True), True, 360)
    
    processQueue.join()
            
    #for rP in refmapProcess:
    #    rP.start()
    #for rP in refmapProcess:
    #    rP.join()  

def refmap_BWA(in_file, fname, out_file, BWA_path, pseudoref_full_path, execute=True):    
    
    #### NEED TO CHECK IF LIBRARIES SPLIT ACROSS LANES (AND IN DIFFERENT FASTQ FILES) HAVE SAMPLES OVERWRITTEN HERE
    #### I SUSPECT THIS IS THE CASE; IF SO FASTQ FILES SHOULD BE CONSOLIDATED BY LIBRARY (JUST CAT THE FILES)
    
    BWAMemTemplate = Template('%s mem -M -R $rgh $p $input > $out' % BWA_path)
    
    print 'Reference mapping ' + fname + '\n'
    
    read_group_header = '"@RG\\tID:' + fname + '\\tPL:Illumina\\tLB:' + fname + '"' 
    bwa_mem_call = BWAMemTemplate.substitute(rgh = read_group_header, input = in_file, p = pseudoref_full_path, out = out_file)
    #bwa_mem_call = BWA_path + ' mem -M -R ' + read_group_header + " " + pseudoref_full_path + ' ' + in_dir + i + ' > ' + out_dir + fname + '.sam'
    
    print bwa_mem_call
    
    #if this is a single run (default), then call BWA from here
    if execute:
        subprocess.call(bwa_mem_call, shell=True)
        return
    
    #if this is a parallel run (multiple individuals to map), pass the cmd line back to parallel_refmap_BWA
    else:
        return bwa_mem_call

def callGeno(sam_in, pseudoref, BCFout, VCFout, samtoolsPath, bcftoolsPath):
    #print sam_in, pseudoref, BCFout, VCFout
    # set up the individual files for transfer from sam to bam and bam indexing
    
    print 'Processing sam files into sorted bam files.'
    
    samtoolsView = Template('%s view -F 4 -b -S -o $output $input' % samtoolsPath)
    samtoolsSort = Template('%s sort -o $output $input' % samtoolsPath)
    samtoolsIndex = Template('%s index $input' % samtoolsPath)
    samtoolsMpileup = Template('%s mpileup -t DP -C50 -u -I -f $reference -o $bcf_out $input' % samtoolsPath)
    bcftoolsView = Template('%s call -v -m $input > $output' % bcftoolsPath)
    
    for sam in os.listdir(sam_in):
        # TODO: add if 'sam' in sam (check that we don't already have bam files)
        # samtools view -F 4 will filter OUT reads with bitwise flag 0004 -- these are unmapped reads
        fname = os.path.splitext(sam)[0]
        
        samPath = sam_in + '/' + sam
        bam = sam_in + '/' + fname + '.bam' # for bam output
        sorted_sam = sam_in + '/' + fname + '.sorted.bam' # for sorting output
        
        view_cmd = samtoolsView.substitute(output = bam, input = samPath)
        subprocess.call(view_cmd, shell=True)
        
        sort_cmd = samtoolsSort.substitute(output = sorted_sam, input = bam)
        subprocess.call(sort_cmd, shell=True)
        
        index_cmd = samtoolsIndex.substitute(input = sorted_sam)
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
