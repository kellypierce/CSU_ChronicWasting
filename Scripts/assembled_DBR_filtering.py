#!/usr/bin/python

############################################################
##                                                        ##
##  Compare DBRs in assembled data to DBRs from raw data  ##
##                                                        ##
############################################################

from collections import defaultdict
from collections import Counter
from subprocess import call, Popen, PIPE
import subprocess
import os as os
from os import linesep, path, R_OK, X_OK
import sys
from logging import debug, critical, error, info
import warnings
import json
import re
import itertools
import fnmatch
#import numpy as np
import time
import pdb
import multiprocessing as mp
from Queue import Queue
from threading import Thread
import warnings
import string as string
from string import Template, join
import logging as logging
import gzip
from collections import defaultdict
from collections import Counter

# To do
# 1. Check that SAM files contain a map for all the sequences so that FASTQ filtering doesn't leave some bad quality data behind
#    Alternatively, instead of tracking tags to remove, track tags to keep
# 2. Make this work for data that span multiple libraries by calling the function separately for each library following the template below:
#    DBR_filter(assembled_dir = '/path/to/assembly_library1',
#               dict_in = '/path/to/dbr_dict_library1',
#               out_seqs = '/path/to/filtered_library1.fastq',
#               n_expected = 2)

################################## GLOBALS ####################################

# Process queue
processQueue = Queue()

phred_dict = {'"':1.0,"#":2.0,"$":3.0,"%":4.0,"&":5.0,"'":6.0,"(":7.0,")":8.0,"*":9.0,"+":10.0,
              ",":11.0,"-":12.0,".":13.0,"/":14.0,"0":15.0,"1":16,"2":17.0,"3":18.0,"4":19.0,"5":20.0,
              "6":21.0,"7":22.0,"8":23.0,"9":24.0,":":25.0,";":26,"<":27.0,"+":28.0,">":29.0,"?":30.0,
              "@":31.0,"A":32.0,"B":33.0,"C":34.0,"D":35.0,"E":36,"F":37.0,"G":38.0,"H":39.0,"I":40.0,
              "J":41.0,"K":42.0}
# conversion reference: http://drive5.com/usearch/manual/quality_score.html

################################# FUNCTIONS ####################################

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
        if isinstance(workItem, filterGroup):
            #raise Exception("correct type")
            workItem.DBR_Filter()
        elif isinstance(workItem, Work):
            p = Popen(workItem.commandline,
                      env = workItem.env,
                      shell = workItem.shell,
                      cwd = workItem.cwd)
            p.wait()
        else:
            raise Exception(type(workItem))
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
    
def find_SampleID(filename, r):
    #sampleID_match = re.match(".*(\d{3}[a-z]?).*", filename)
    # this revision is VERY specific to my technical replicates
    # TODO: find a way to pass the regex capture group as an argument so that this (and related functions) are more flexible
    #sampleID_match = re.match(".*(\d{1,3}T?).*", filename)
    sampleID_match = re.match(r, filename)
    if sampleID_match:
        sampleID = sampleID_match.groups()[0]
        return sampleID
    else:
        return None
    
def find_LibraryID(filename, r):
    #libraryID_match = re.match(".*(Library\d{2,3}).*", filename)
    #libraryID_match = re.match(".*(Library\d{1,3}[A|B]?).*", filename)
    libraryID_match = re.match(r, filename)
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
    
def find_DBRdictionary(library, directory):
    if library: # library can also be returned as 'None' for files with improper naming
        if os.path.isdir(directory):
            dcs = os.listdir(directory)
            for d in dcs:
                if library in d:
                    dcf = directory + '/' + d
                    return dcf
        else: # if it's just a single file
            if os.path.isfile(directory):
                return directory
            else:
                return None
    else:
        return None

#processQueue.put(Work(commandline = commandline, shell = True), True, 360)
    
#processQueue.join()

# def barcodeDict(barcode_file):
#     
#     bc_dict = {} # define empty container for barcode dictionary
#             
#     with open(barcode_file, 'r') as bc:
#         for line in bc:
#             row=line.split()
#             bc_dict[row[0]]=row[1]
#     
#     return bc_dict

def assemblyDict(DBR_dict, sampleID, assembled_file, samMapLen):#, barcode_dict):
    
    dbr = json.load(DBR_dict)
    #original_barcode = barcode_dict[sampleID]
    
    # initialize an empty dictionary with each iteration of the for-loop
    assembly_dict_2 = {}
    assembly_dict_3 = defaultdict(list)
    
    # print some info to track progress
    #path=os.path.join(assembled_dir, i)
    print 'Creating filtering dictionaries from ' + assembled_file
    
    # get the sample number for use later
    #number = re.split('(\d)', i)[1] # enclose the regex in parentheses to keep it in the output
    
    # start counter for the number of primary reads
    n_primary = 0
    
    # open the sam file and process its contents
    with open(assembled_file, 'r') as inFile:
        for line in inFile:
            if not line.startswith("@"): # ignore the header lines
                
                fields = line.split("\t")
    
                # extract the info for the dictionary for each line
                QNAME = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', fields[0])[1] # FASTQ ID = QNAME column in sam file
                FLAG = fields[1] # bitwise flag with map info; == 0 if primary read
                RNAME = fields[2] # ref sequence name -- where did the sequence map?
                POS = fields[3] # position of map
                MAPQ = fields[4] # mapping quality
                CIGAR = fields[5] # additional mapping info
                SEQ = fields[9] # the actual sequence
                QUAL = fields[10] # sequence quality score
                
                # extract the DBR corresponding to the QNAME for each row
                dbr_value = dbr.get(QNAME) 
                
                # bitwise FLAG == 0 means that the read is the PRIMARY READ. There will only be one of these per sequence, so only mapped primary reads should be considered.
                if FLAG == '0':
                    if samMapLen:
                        if len(SEQ) == samMapLen: #if we specify an expected sequence length in the samfile
                            # tally the new primary read
                            #n_primary += 1
                        
                            # WE NEED TWO DICTIONARIES TO REPRESENT ALL THE RELATIONSHIP BETWEEN RNAME, QNAME, dbr_value, QUAL, AND count
                            # build a dictionary with structure {DBR: (locus: count)}                    
                            if RNAME in assembly_dict_2:
                                if dbr_value in assembly_dict_2.get(RNAME):
                                    assembly_dict_2[RNAME][dbr_value]=assembly_dict_2[RNAME][dbr_value]+1
                                else:
                                    assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1
                            else:
                                assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1 # add the new DBR and its associated locus and count
                
                            # build a dictionary with structure {RNAME: {DBR:[[QNAME, QUAL]]}}        
                            if RNAME in assembly_dict_3:
                                if dbr_value in assembly_dict_3.get(RNAME):
                                    assembly_dict_3[RNAME][dbr_value].append([QNAME, QUAL, SEQ])
                                else:
                                    assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL, SEQ]]
                            else:
                                assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL, SEQ]]
                            # tally the new primary read
                            n_primary += 1
                    else: #if we're not using stacks to re-assemble and we don't care about expected lengths...
                        # WE NEED TWO DICTIONARIES TO REPRESENT ALL THE RELATIONSHIP BETWEEN RNAME, QNAME, dbr_value, QUAL, AND count
                        # build a dictionary with structure {DBR: (locus: count)}                    
                        if RNAME in assembly_dict_2:
                            if dbr_value in assembly_dict_2.get(RNAME):
                                assembly_dict_2[RNAME][dbr_value]=assembly_dict_2[RNAME][dbr_value]+1
                            else:
                                assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1
                        else:
                            assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1 # add the new DBR and its associated locus and count
                
                        # build a dictionary with structure {RNAME: {DBR:[[QNAME, QUAL]]}}        
                        if RNAME in assembly_dict_3:
                            if dbr_value in assembly_dict_3.get(RNAME):
                                assembly_dict_3[RNAME][dbr_value].append([QNAME, QUAL, SEQ])
                            else:
                                assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL, SEQ]]
                        else:
                            assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL, SEQ]]
                        # tally the new primary read
                        n_primary += 1
                        
    return assembly_dict_2, assembly_dict_3, n_primary

class filterGroup():
    def __init__(self,
                 assembled_file, # the SAM files for the data mapped to pseudoreference
                 out_file, # the output file, full path, ending with .fasta
                 n_expected, # the number of differences to be tolerated
                 dict_file, # a single dictionary of DBRs (for one library only) 
                 sampleID,
                 logfile,
                 test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
                 phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
                 samMapLen=None):
        
        self.assembled_file = assembled_file
        self.out_file = out_file
        self.n_expected = n_expected
        self.dict_file = dict_file
        self.sampleID = sampleID
        self.logfile = logfile
        self.test_dict = test_dict
        self.phred_dict = phred_dict
        self.samMapLen = samMapLen

#     def DBR_Filter(assembled_file, # the SAM files for the data mapped to pseudoreference
#                    out_file, # the output file, full path, ending with .fasta
#                    n_expected, # the number of differences to be tolerated
#                    dict_file, # a single dictionary of DBRs (for one library only)
#                    sampleID,
#                    logfile,
#                    test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
#                    phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
#                    samMapLen=None): # expected sequence length will help when primary reads are still not perfectly aligned with reference
     
    def DBR_filter(self):
            
        print 'Opening DBR dictionary ' + self.dict_file  
        
        # open a connection to an output file
        with open(self.out_file, 'a') as of:
            
            # probably unnecessary ...
            #bc_dict = barcodeDict(barcode_file)
            
            # return two dictionaries for use in filtering
            assembly_dictionaries = assemblyDict(self.dict_file, self.sampleID, self.assembled_file, self.samMapLen)
            
            # separate the two dictionaries
            assembly_dict_2 = assembly_dictionaries[0]
            assembly_dict_3 = assembly_dictionaries[1]
            n_primary = assembly_dictionaries[2]
            
            print 'Checking DBR counts against expectations.'
            
            total_removed = 0
            delete_list = []
            keep_list = []
            
            for RNAME, value in assembly_dict_2.iteritems():
                #print 'RNAME', RNAME
                # ignore the data where the reference is "unmapped" -- RNAME = '*'
                if RNAME != '*':
                    # get all the DBRs and counts that went into that locus in that sample
                    for subvalue in value.iteritems():
                        #print 'Subvalue', subvalue
                        dbr_value = subvalue[0] 
                        count = subvalue[1]
                        if count > self.n_expected:
                            ##################################################
                            ## THIS IS WHERE THE FILTERING HAPPENS           #
                            ##################################################
                            #print 'count', count, 'n exp', n_expected
                            # the other dictionary contains the full quality information for each RNAME:DBR pair (this will be multiple entries of sequence IDs and qualities
                            qname_qual = assembly_dict_3[RNAME][dbr_value] #this is a list of lists: [[QNAME, QUAL, SEQ], [QNAME, QUAL, SEQ], ...]
                            #print RNAME, dbr_value, len(qname_qual)
                            #print qname_qual
                            ID_quals = {} # we'll make yet another dictionary to store the QNAME and the median QUAL
                            for i in qname_qual:
                                id_val=i[0] # this is the QNAME (Illumina ID)
                                id_seq=i[2] # the full sequence
                                id_qual=i[1] # the full quality
                                ID_quals[id_val] = (qual_median(i[1], phred_dict), id_seq, id_qual)
                            n_remove = count - self.n_expected
                            total_removed += n_remove
                            t = 0
                            k = 1
                            while k <= self.n_expected:
                                to_keep = max(ID_quals, key=lambda x:ID_quals[x]) 
                                #print(to_keep)
                                keep = ID_quals[to_keep]
                                keep_list.append(to_keep)
                                #write out the data to keep, NOT appending the original barcode to the beginning of the sequence
                                of.write('@'+to_keep+'\n'+ keep[1]+'\n+\n'+ keep[2]+'\n')
                                #out_file.write([keep.split('\n', 1)[0] for i in keep])
                                k += 1
            with open(self.logfile,'a') as log:
                log.write(self.sampleID+','+str(total_removed)+','+str(n_primary)+','+time.strftime("%d/%m/%Y")+','+(time.strftime("%H:%M:%S"))+'\n')
                print 'Removed ' + str(total_removed) + ' PCR duplicates out of ' + str(n_primary) + ' primary mapped reads.'                                                    
                    
            if self.test_dict: # check construction by printing first entries to screen
                print 'Checking dictionary format (version 3).'
                x = itertools.islice(assembly_dict_3.iteritems(), 0, 4)
                for keyX, valueX in x:
                    print keyX, valueX
                print 'Checking dictionary format (version 2).'
                y = itertools.islice(assembly_dict_2.iteritems(), 0, 4)
                for keyY, valueY in y:
                    print keyY, valueY
               
def parallel_DBR_Filter(assembled_dir, # the SAM files for the data mapped to pseudoreference
               out_dir, # the output file, full path, ending with .fasta
               n_expected, # the number of differences to be tolerated
               dict_dir, # a single dictionary of DBRs (for one library only)
               regexSample,
               regexLibrary,
               test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
               phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
               samMapLen=None): # expected sequence length will help when primary reads are still not perfectly aligned with reference
    
    logfile = out_dir + '/DBR_filtered_sequences_logfile.csv' 
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    for i in os.listdir(assembled_dir):
            
        if 'unmatched' not in i: # skip the SAM files with sequences that didn't match
            
            # extract the sample ID with a regex
            sampleID = find_SampleID(i, regexSample)
            
            # extract the library ID with a regex
            libraryID = find_LibraryID(i, regexLibrary)
            
            # use the library ID to find the right DBR dictionary
            dict_in = find_DBRdictionary(libraryID, dict_dir)
            
            if sampleID and libraryID and dict_in: # and bcf ... (deprecated) # if all of these != None
            
                print 'sample', sampleID 
                print 'library', libraryID
                print 'dictionary file', dict_in
                
                inFile = os.path.join(assembled_dir, i)
                # one output fastq file per sample    
                out_seqs_final = out_dir + '/DBR_filtered_sequences_' + libraryID + '_' + sampleID + '.fastq'
                
                processQueue.put(filterGroup(assembled_file = inFile,
                                         out_file = out_seqs_final,
                                         n_expected = n_expected,
                                         dict_file = dict_in,
                                         sampleID = sampleID,
                                         logfile = logfile))
    processQueue.join()

    return

def parallel_DBR_dict(in_dir, seqType, dbr_start, dbr_stop, test_dict = False, save = None):
    #if not checkDir(in_dir):
    #    raise IOError("Input is not a directory: %s" % in_dir)
    if seqType == 'read2':
        warnings.warn('Expect directory containing only Read 2 files; any other files present in %s will be incorporated into DBR dictionary.' % in_dir)
    elif seqType == 'pear':
        warnings.warn('Expect directory containing only merged Read 1 and Read 2 files; any other files present in %s will be incorporated into DBR directory' % in_dir)
    else:
        raise IOError("Input sequence type specified as %s. Options are 'pear' or 'read2'." % seqType)
    
    dbrProcess = [mp.Process(target=DBR_dict, args=(in_dir+in_file, 
                                                    seqType,
                                                    dbr_start,
                                                    dbr_stop,
                                                    test_dict,
                                                    save)) for in_file in in_dir]
     
    for dP in dbrProcess:
        dP.start()
    for dP in dbrProcess:
        dP.join()
    

def DBR_dict(in_file, dbr_start, dbr_stop, test_dict = False, save = None):
    # DBR is in read 2
    # if merged, it will be the last -2 to -9 (inclusive) bases, starting with base 0 and counting from the end
    # if not merged, it will be bases 2 to 9
    if not checkFile(in_file):
        raise IOError("where is the input file: %s" % in_file)
    info('Creating {ID: dbr} dictionary from %s.' % in_file)
    dbr = {}
    fq_line = 1
    if in_file.endswith('gz'):
        openFxn = gzip.open
    else:
        openFxn = open
    with openFxn(in_file, 'r') as db:
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
        fq_name = os.path.splitext(in_file)[0]
        fq_dbr_out = fq_name + save + '.json'
        print 'Writing dictionary to ' + fq_dbr_out
        with open(fq_dbr_out, 'w') as fp:          
            json.dump(dbr, fp)
            
'''
def DBR_Filter(assembled_dir, # the SAM files for the data mapped to pseudoreference
               out_dir, # the output file, full path, ending with .fasta
               n_expected, # the number of differences to be tolerated
               barcode_dir, # the barcodes for individuals in the library referenced in dict_in
               dict_dir, # a single dictionary of DBRs (for one library only)
               sample_regex, # regular expression to find the sample ID
               barcode_file=None, # if just a single library is being used, can directly pass the barcode file
               test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
               phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
               samMapLen=None): # expected sequence length will help when primary reads are still not perfectly aligned with reference
    
    # addition: grep the name of the file in the assembled_dir for a number in column 1 of the barcode file
    # get the barcode and append it to the front of the sequence
    # if we don't do this, the demultiplexing won't work for phase 2
    # alternatively we can just give a new barcode, but then we'd need to track a secondary barcode file for all the individuals
    # that sounds not-fun
    
    #pdb.set_trace()
    #logfile = os.path.splitext(out_seqs)[0] + '_logfile.csv'
    logfile = out_dir + '/DBR_filtered_sequences_logfile.csv'
    
    # for each sample -- each file in assembled_dir is a sam file for a single sample
    for i in os.listdir(assembled_dir):
            
        if 'unmatched' not in i: # skip the SAM files with sequences that didn't match
            
            print i
            
            # extract the sample ID with a regex
            sampleID = find_SampleID(i, sample_regex) # find the sample ID, potentially with some extra characters to distinguish from library ID
            
            # extract the library ID with a regex
            libraryID = find_LibraryID(i)
            
            # use the library ID to find the right barcode file
            bcf = find_BarcodeFile(libraryID, barcode_dir)
            
            # use the library ID to find the right DBR dictionary
            dict_in = find_DBRdictionary(libraryID, dict_dir)
            
            if sampleID and libraryID and bcf and dict_in: # if all of these != None
            
                print 'sample', sampleID 
                print 'library', libraryID
                print 'barcode file', bcf 
                print 'dictionary file', dict_in
            
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
            
                out_seqs_final = out_dir + '/DBR_filtered_sequences_' + libraryID + '.fastq'
            
                #os.path.isfile(fname)
    
                with open(out_seqs_final, 'a') as out_file:
                    
                    bc_dict = {} # define empty container for barcode dictionary
            
                    with open(bcf, 'r') as bc:
                        for line in bc:
                            row=line.split()
                            bc_dict[row[0]]=row[1]
            
                        print 'Opening DBR dictionary ' + dict_in  
                        with open(dict_in, 'r') as f:
                            dbr = json.load(f)
                            
                            original_barcode = bc_dict[sampleID]
                            #print original_barcode
                            # suggestion on error checking: 
                            # normally i capture the .match() value and say 'if object:"
                            # "else: print('did not match correctly')
                            
                            # initialize an empty dictionary with each iteration of the for-loop
                            assembly_dict_2 = {}
                            assembly_dict_3 = defaultdict(list)
                            
                            # print some info to track progress
                            path=os.path.join(assembled_dir, i)
                            print 'Creating filtering dictionaries from ' + path
                            
                            # get the sample number for use later
                            #number = re.split('(\d)', i)[1] # enclose the regex in parentheses to keep it in the output
                            
                            # start counter for the number of primary reads
                            n_primary = 0
                            
                            delete_list = []
                            keep_list = []
                            
                            # open the sam file and process its contents
                            with open(path, 'r') as inFile:
                                for line in inFile:
                                    if not line.startswith("@"): # ignore the header lines
                                        fields = line.split("\t")
                                        
                                        # extract the info for the dictionary for each line
                                        QNAME = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', fields[0])[1] # FASTQ ID = QNAME column in sam file
                                        FLAG = fields[1] # bitwise flag with map info; == 0 if primary read
                                        RNAME = fields[2] # ref sequence name -- where did the sequence map?
                                        POS = fields[3] # position of map
                                        MAPQ = fields[4] # mapping quality
                                        CIGAR = fields[5] # additional mapping info
                                        SEQ = fields[9] # the actual sequence
                                        QUAL = fields[10] # sequence quality score
                                        
                                        # extract the DBR corresponding to the QNAME for each row
                                        dbr_value = dbr.get(QNAME) 
                                        
                                        # after trimming, the sequences are 116 bases long
                                        # the POS for all matched sequences should be 1 because we used a pseudoreference build from our RADtags
                                        # if all 116 bases match, POS = 1 and CIGAR = 116M (meaning 116 bases match the reference)
                                        # repetitive regions may map to multiple regions but not exactly -- this would cause the same sequence ID to be present > 1x in the dictionary
                                        # sequences present >1x in the dictionary break the count of DBRs, and cause problems with filtering downstream
                                        # the inexact matches of repetitive regions can be detected by POS != 1 or CIGAR != 116M
                                        # a more flexible/broader use way of filtering would be to do a REGEX search on the CIGAR score and report the number of digits preceeding the M
                                        # then you could only keep the entry that has the largest number of matches (but M can show up multiple times in the CIGAR score if separated by an insertion, so think about this more)
                                        # I had thought that filtering on POS == 1 will keep only the good matches, but it is possible to have multiple matches to POS == 1 with different levels of clipping
                                        
                                        # A MORE GENERAL SOLUTION THAT SHOULD WORK FOR A VARIETY OF CIRCUMSTANCES:
                                        # bitwise FLAG == 0 means that the read is the PRIMARY READ. There will only be one of these per sequence, so only mapped primary reads should be considered.
                                        if FLAG == '0':
                                            if samMapLen:
                                                if len(SEQ) == samMapLen: #if we specify an expected sequence length in the samfile
                                                    # tally the new primary read
                                                    #n_primary += 1
                                                
                                                    # WE NEED TWO DICTIONARIES TO REPRESENT ALL THE RELATIONSHIP BETWEEN RNAME, QNAME, dbr_value, QUAL, AND count
                                                    # build a dictionary with structure {DBR: (locus: count)}                    
                                                    if RNAME in assembly_dict_2:
                                                        if dbr_value in assembly_dict_2.get(RNAME):
                                                            assembly_dict_2[RNAME][dbr_value]=assembly_dict_2[RNAME][dbr_value]+1
                                                        else:
                                                            assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1
                                                    else:
                                                        assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1 # add the new DBR and its associated locus and count
                                    
                                                    # build a dictionary with structure {RNAME: {DBR:[[QNAME, QUAL]]}}        
                                                    if RNAME in assembly_dict_3:
                                                        if dbr_value in assembly_dict_3.get(RNAME):
                                                            assembly_dict_3[RNAME][dbr_value].append([QNAME, QUAL, SEQ])
                                                        else:
                                                            assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL, SEQ]]
                                                    else:
                                                        assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL, SEQ]]
                                                    # tally the new primary read
                                                    n_primary += 1
                                            else: #if we're not using stacks to re-assemble and we don't care about expected lengths...
                                                # WE NEED TWO DICTIONARIES TO REPRESENT ALL THE RELATIONSHIP BETWEEN RNAME, QNAME, dbr_value, QUAL, AND count
                                                # build a dictionary with structure {DBR: (locus: count)}                    
                                                if RNAME in assembly_dict_2:
                                                    if dbr_value in assembly_dict_2.get(RNAME):
                                                        assembly_dict_2[RNAME][dbr_value]=assembly_dict_2[RNAME][dbr_value]+1
                                                    else:
                                                        assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1
                                                else:
                                                    assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1 # add the new DBR and its associated locus and count
                                
                                                # build a dictionary with structure {RNAME: {DBR:[[QNAME, QUAL]]}}        
                                                if RNAME in assembly_dict_3:
                                                    if dbr_value in assembly_dict_3.get(RNAME):
                                                        assembly_dict_3[RNAME][dbr_value].append([QNAME, QUAL, SEQ])
                                                    else:
                                                        assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL, SEQ]]
                                                else:
                                                    assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL, SEQ]]
                                                # tally the new primary read
                                                n_primary += 1
                                    
                                # NOW THAT DICTIONARIES ARE MADE, REMOVE DUPLICATE SEQUENCES BASED ON DBR COUNTS
                                # for each assembled locus, get the associated dbr_value and count
                                print 'Checking DBR counts against expectations.'
                                total_removed = 0
                                for RNAME, value in assembly_dict_2.iteritems():
                                    #print 'RNAME', RNAME
                                    # ignore the data where the reference is "unmapped" -- RNAME = '*'
                                    if RNAME != '*':
                                        # get all the DBRs and counts that went into that locus in that sample
                                        for subvalue in value.iteritems():
                                            #print 'Subvalue', subvalue
                                            dbr_value = subvalue[0] 
                                            count = subvalue[1]
                                            if count > n_expected:
                                                ##################################################
                                                ## THIS IS WHERE THE FILTERING HAPPENS           #
                                                ##################################################
                                                #print 'count', count, 'n exp', n_expected
                                                # the other dictionary contains the full quality information for each RNAME:DBR pair (this will be multiple entries of sequence IDs and qualities
                                                qname_qual = assembly_dict_3[RNAME][dbr_value] #this is a list of lists: [[QNAME, QUAL, SEQ], [QNAME, QUAL, SEQ], ...]
                                                #print RNAME, dbr_value, len(qname_qual)
                                                #print qname_qual
                                                ID_quals = {} # we'll make yet another dictionary to store the QNAME and the median QUAL
                                                for i in qname_qual:
                                                    id_val=i[0] # this is the QNAME (Illumina ID)
                                                    id_seq=i[2] # the full sequence
                                                    id_qual=i[1] # the full quality
                                                    ID_quals[id_val] = (qual_median(i[1], phred_dict), id_seq, id_qual)
                                                n_remove = count - n_expected
                                                total_removed += n_remove
                                                t = 0
                                                k = 1
                                                while k <= n_expected:
                                                    to_keep = max(ID_quals, key=lambda x:ID_quals[x]) 
                                                    #print(to_keep)
                                                    keep = ID_quals[to_keep]
                                                    keep_list.append(to_keep)
                                                    #write out the data to keep, appending the original barcode to the beginning of the sequence
                                                    out_file.write('@'+to_keep+'\n'+ original_barcode+keep[1]+'\n+\n'+ "KKKKK"+keep[2]+'\n')
                                                    #out_file.write([keep.split('\n', 1)[0] for i in keep])
                                                    k += 1
                                with open(logfile,'a') as log:
                                    log.write(sampleID+','+str(total_removed)+','+str(n_primary)+','+time.strftime("%d/%m/%Y")+','+(time.strftime("%H:%M:%S"))+'\n')
                                    print 'Removed ' + str(total_removed) + ' PCR duplicates out of ' + str(n_primary) + ' primary mapped reads.'                                                    
                                        
                                if test_dict: # check construction by printing first entries to screen
                                    print 'Checking dictionary format (version 3).'
                                    x = itertools.islice(assembly_dict_3.iteritems(), 0, 4)
                                    for keyX, valueX in x:
                                        print keyX, valueX
                                    print 'Checking dictionary format (version 2).'
                                    y = itertools.islice(assembly_dict_2.iteritems(), 0, 4)
                                    for keyY, valueY in y:
                                        print keyY, valueY
        
'''
#TODO: why does DBR_filter need to write out a single fastq file -- why redo all that demultiplexing??
'''
# OTHER METRICS FOR DESCRIBING OVERALL SEQUENCE QUALITY (QUAL = ASCII character string)

# counter object (from collections import Counter)
countQUAL = Counter(QUAL)

# frequency table of ASCII characters
freqQUAL = countQUAL.most_common()

# most frequently observed ASCII character
modeQUAL = countQUAL.most_common(1)

# most frequently observed integer score
intQUAL = phred_dict[modeQUAL[0][0]]

# list of integer qualities
listQUAL = list(QUAL) # split the ASCII string into a list
list_intQUAL = []
for q in listQUAL:
    list_intQUAL.append(phred_dict[q])
    
# median quality
medQUAL = np.median(list_intQUAL)
'''