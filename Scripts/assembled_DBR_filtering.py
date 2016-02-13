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
import os, os.path
import sys
import json
import re
import itertools
#import numpy as np
import time

# To do
# 1. Check that SAM files contain a map for all the sequences so that FASTQ filtering doesn't leave some bad quality data behind
#    Alternatively, instead of tracking tags to remove, track tags to keep
# 2. Make this work for data that span multiple libraries by calling the function separately for each library following the template below:
#    DBR_filter(assembled_dir = '/path/to/assembly_library1',
#               dict_in = '/path/to/dbr_dict_library1',
#               out_seqs = '/path/to/filtered_library1.fastq',
#               n_expected = 2)


phred_dict = {'"':1.0,"#":2.0,"$":3.0,"%":4.0,"&":5.0,"'":6.0,"(":7.0,")":8.0,"*":9.0,"+":10.0,
              ",":11.0,"-":12.0,".":13.0,"/":14.0,"0":15.0,"1":16,"2":17.0,"3":18.0,"4":19.0,"5":20.0,
              "6":21.0,"7":22.0,"8":23.0,"9":24.0,":":25.0,";":26,"<":27.0,"+":28.0,">":29.0,"?":30.0,
              "@":31.0,"A":32.0,"B":33.0,"C":34.0,"D":35.0,"E":36,"F":37.0,"G":38.0,"H":39.0,"I":40.0,
              "J":41.0,"K":42.0}
# conversion reference: http://drive5.com/usearch/manual/quality_score.html

''' stuff from the pilot library
dict_in = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/initial_qualFilter_dbr_dict'
assembled_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Assembled'
out_seqs = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/DBR_filtered_sequences.fastq'
barcode_file = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/pilot_barcode_file'
'''

def qual_mode(QUAL, phred_dict):
    listQUAL = list(QUAL)
    list_intQUAL =[]
    for q in listQUAL:
        list_intQUAL.append(phred_dict[q])
    #return np.median(list_intQUAL)
    qsort = list_intQUAL.sort()
    qlen = len(qsort)
    if qlen % 2 == 0: # even length list -- take the average of the two middle values
        median_qual = (qsort[(qlen/2)-1]+qsort[(qlen/2)])/2
    else: # odd length list -- take the middle value (list length / 2 will be automatically rounded up by python integer division)
        median_qual = qsort[qlen/2]
    return median_qual

'''
def iterative_DBR_filter(in_dir, # directory of ref-mapped sequence data
                         barcode_dir,
                         out_dir, # full path for output fastq file (ALL filtered sequences go to this single file
                         out_prefix, # text string to add to output file (e.g., 'dbr_filtered.fastq')
                         dict_dir, # a directory of DBR dictionaries
                         out_seqs, # the output file, full path, ending with .fasta
                         n_expected, # the number of differences to be tolerated
                         test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
                         phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
                         samMapLen=None): # expected sequence length will help when primary reads are still not perfectly aligned with reference
    #(assembled_dir, dict_in, out_seqs, n_expected, barcode_file, test_dict=True, phred_dict=phred_dict, samMapLen=None)
    #if not checkDir(in_dir):
    #    raise IOError("Input is not a directory: %s" % in_dir)
    #if not checkFile(barcode_file):
    #    raise IOError("Where is the barcode file? %s" % barcode_file)
    #pdb.set_trace()
    files = os.listdir(in_dir) 
    for f in files: # for each individual sam file
        sampleID_match = re.match(".*(Library\d{2,3}).*", f) # get the bit of the file name that links it back up to the barcode file
        #sampleID = re.match(".*(\d{3}[a-z]?).*", f).groups()[0]
        if sampleID_match: # if we get a match
            sampleID = sampleID_match.groups()[0] # extract that match
            bcs = os.listdir(barcode_dir)
            dicts = os.listdir(dict_dir)
            for b in bcs:
                if sampleID in b: # find the corresponding barcode file
                    for d in dicts:
                        if sampleID in d: # find the corresponding dictionary file
                            barcode_file = barcode_dir + '/' + b
                            in_f = in_dir + '/' + f
                            dict_in = dict_dir + '/' + d
                            out_name = out_seqs + '/' + sampleID + '.fq' # output the filtered sequences to a unique file by library name
                            DBR_Filter(assembled_dir = in_f, 
                                       dict_in = dict_in , # a single dictionary of DBRs (for one library only)
                                       out_seqs = out_seqs, # the output file, full path, ending with .fasta
                                       n_expected = 2, # the number of differences to be tolerated
                                       barcode_file = barcode_file, # the barcodes for individuals in the library referenced in dict_in
                                       test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
                                       phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
                                       samMapLen=None)
'''

def DBR_Filter(assembled_dir, # the SAM files for the data mapped to pseudoreference
               out_seqs, # the output file, full path, ending with .fasta
               n_expected, # the number of differences to be tolerated
               barcode_dir, # the barcodes for individuals in the library referenced in dict_in
               dict_dir, # a single dictionary of DBRs (for one library only)
               barcode_file=None, # if just a single library is being used, can directly pass the barcode file
               test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
               phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
               samMapLen=None): # expected sequence length will help when primary reads are still not perfectly aligned with reference
    
    # addition: grep the name of the file in the assembled_dir for a number in column 1 of the barcode file
    # get the barcode and append it to the front of the sequence
    # if we don't do this, the demultiplexing won't work for phase 2
    # alternatively we can just give a new barcode, but then we'd need to track a secondary barcode file for all the individuals
    # that sounds not-fun
    
    logfile = os.path.splitext(out_seqs)[0] + '_logfile.csv'
    
    with open(out_seqs, 'w') as out_file: 
        # for each sample -- each file in assembled_dir is a sam file for a single sample
        for i in os.listdir(assembled_dir):
             
            bcf = '' # define empty container for barcode text file name
            dict_in = '' # define an empty container for the DBR dictionary file name
            bc_dict = {} # define empty container for barcode dictionary
            
            if barcode_file != None: # if we just have a single barcode file, then we don't need to match file names
                
                bcf = barcode_file
                dict_in = dict_dir # add some error checking here eventually: if there's just one barcode file there MUST be only one DBR dictionary
                
            elif barcode_file == None: # if we're processing data from several libraries
                # find which library the sample was in
                if 'unmatched' not in i: # don't process the unmatched sequences
                    
                    libraryID_match = re.match(".*(Library\d{2,3}).*", i)
                
                    if libraryID_match: # if we get a match (this allows the script to proceed if a file has a mismatched name)
                        
                        libraryID = libraryID_match.groups()[0] # extract the library ID match
                                    
                        # list all the files in the barcode directory and find the match
                        bcs = os.listdir(barcode_dir)
                        for b in bcs:
                            if libraryID in b:
                                bcf = barcode_dir + '/' + b
                        
                        dcs = os.listdir(dict_dir)
                        for d in dcs:
                            if libraryID in d:
                                dict_in = dict_dir + '/' + d
                        
            # make a dictionary of the barcodes in the barcode file
            # this will be re-done for each sample
            # in the future it may be faster to batch by library ID somehow (search for unique library IDs in the files and then process by library?)      
            with open(bcf, 'r') as bc:
                for line in bc:
                    row=line.split()
                    bc_dict[row[0]]=row[1]

            print bc_dict
            # find the sample ID
            sampleID = re.match(".*(\d{3}[a-z]?).*", i).groups()[0]
            print libraryID, sampleID
            
            print 'Opening DBR dictionary ' + dict_in  
            with open(dict_in, 'r') as f:
                dbr = json.load(f)
                '''
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
                                        ID_quals[id_val] = (qual_mode(i[1], phred_dict), id_seq, id_qual)
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
        
'''
DBR_Filter(assembled_dir=assembled_dir, 
           dict_in=dict_in, 
           out_seqs=out_seqs, 
           n_expected = 2, 
           barcode_file=barcode_file)
'''

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