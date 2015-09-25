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
import numpy as np

# To do
# 1. Check that SAM files contain a map for all the sequences so that FASTQ filtering doesn't leave some bad quality data behine
#    Alternatively, instead of tracking tags to remove, track tags to keep
# 2. Make the FASTQ filtering function
#    line = 1
#    with open(filtered_undemultiplexed_fastq, 'r') as raw:
#        with open(dbr_filtered_output_file, 'a') as out:
#            for line in raw:
#                if line == 1:
#                    ID = regex for ID
#                    if ID in long_list_of_IDs_to_save:
#                        out.write(ID)
#        ... and so on... this needs some work, but should move through lines in sets of 4
# 2. Make this work for data that span multiple libraries
#    Something like 
#       for inner_dir in outer_dir: # directory holding directories of demultiplexed data from multiple libraries
#           for file in inner_dir:  # the demultiplexed data themselves
#                do the de novo assembly to get one big catalogs.tsv file for the pseudoreference
# 3. Do we need to make the big "dbr_dict" from the qual-fitered, un-demultiplexed FASTQ?
#    It seems that all the data we need are in the dictionaries built as part of filtering.
#    As long as we have some external representation of the DBR frequencies, we don't really need a separate dictionary
# 4. How should this be executed? Many cmd line params is tedious, and that would usurp the modular nature of this code.
#    Instead I think that each execution should be its own python script, e.g. 'pilot_assemly.py'
#    These python scripts would import DBR_Parsing, integrated_denovo_pipeline, and assembled_DBR_filtering python modules.
#    After the imports, there would be a series of calls to the relevant functions.
#    An added benefit is that each "job script" is saved -- you'll know exactly what you ran and on what date.
# 5. On that note, make this stuff functions!


def qual_mode(QUAL, phred_dict):
    listQUAL = list(QUAL)
    list_intQUAL =[]
    for q in listQUAL:
        list_intQUAL.append(phred_dict[q])
    return np.median(list_intQUAL)

phred_dict = {'"':1,"#":2,"$":3,"%":4,"&":5,"'":6,"(":7,")":8,"*":9,"+":10,
              ",":11,"-":12,".":13,"/":14,"0":15,"1":16,"2":17,"3":18,"4":19,"5":20,
              "6":21,"7":22,"8":23,"9":24,":":25,";":26,"<":27,"+":28,">":29,"?":30,
              "@":31,"A":32,"B":33,"C":34,"D":35,"E":36,"F":37,"G":38,"H":39,"I":40,
              "J":41}

# conversion reference: http://drive5.com/usearch/manual/quality_score.html

test_dict=True
write_dict=False
save_reduced_dict=False

# load Dictionary

dict_in = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/initial_qualFilter_dbr_dict'
#seq_dict_in = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/initial_qualFilter_R1_dict'

print 'Opening DBR dictionary ' + dict_in  
with open(dict_in, 'r') as f:
    dbr = json.load(f)

#print 'Opening Read 1 dictionary ' + seq_dict_in
#with open(seq_dict_in, 'r') as s:
#    R1 = json.load(s)

assembled_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Assembled'
for i in os.listdir(assembled_dir):
    
    # initialize an empty dictionary with each iteration of the for-loop
    assembly_dict_2 = {}
    assembly_dict_3 = defaultdict(list)
    
    # print some info to track progress
    path=os.path.join(assembled_dir, i)
    print 'Creating filtering dictionaries from ' + path
    
    # get the sample number for use later
    number = re.split('(\d+)', i)[1] # enclose the regex in parentheses to keep it in the output
    
    # start counter for the number of primary reads
    n_primary = 0
    
    delete_list = []
    
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
                
                    # tally the new primary read
                    n_primary += 1
                
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
                            assembly_dict_3[RNAME][dbr_value].append([QNAME, QUAL])
                        else:
                            assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL]]
                    else:
                        assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL]]

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
                # we know the DBRs with count = 1 are okay to save, so we don't need to calculate expected occurrences for them
                if count > 1:
                    ##################################################
                    ## THIS IS WHERE THE LIKELIHOOD FUNCTION WILL GO #
                    ## FOR NOW THE EXPECTATION IS RANDOM             #
                    ##################################################
                    #n_expected = np.random.uniform(2, 10)
                    n_expected = 5
                    # if we count more DBR occurrences than we expect, we need to evaluate the quality of the sequences associated with that DBR
                    # calculating median QUAL from the full ASCII QUAL score is computationally intensive enough that we only want to do it when absolutely necessary -- only if we observe a DBR more often than expected
                    if count > n_expected:
                        #print 'count', count, 'n exp', n_expected
                        # the other dictionary contains the full quality information for each RNAME:DBR pair (this will be multiple entries of sequence IDs and qualities
                        qname_qual = assembly_dict_3[RNAME][dbr_value] #this is a list of lists: [[QNAME, QUAL], [QNAME, QUAL], ...]
                        #print RNAME, dbr_value, len(qname_qual)
                        #print qname_qual
                        ID_quals = {} # we'll make yet another dictionary to store the QNAME and the median QUAL
                        for i in qname_qual:
                            id_val=i[0]
                            ID_quals[id_val] = qual_mode(i[1], phred_dict)
                        # determine how many sequences will need to be removed and start a counter, t
                        #print 'start length', len(ID_quals)
                        #if len(ID_quals) < count:
                            #print qname_qual
                        n_remove = count - n_expected
                        total_removed += n_remove
                        t = 0
                        #print 'count', count, 'expected', n_expected, 'remove', n_remove, 'tally', t
                        # with the full {ID: qual} dictionary available, we can now determine which sequences should be removed based on their score                        
                        while t < n_remove:
                            t += 1 
                            to_remove = min(ID_quals, key=lambda x:ID_quals[x])
                            # remove the sequences by their IDs from the master R1 sequence ID file
                            del ID_quals[to_remove]
                            # exit with error if a key is missing
                            delete_list.append(to_remove)
                            #assert to_remove in R1, "%s was not in R1" % to_remove
                            #del R1[to_remove] #8:1311:18936:26586
    print 'Removed ' + str(total_removed) + ' PCR duplicates out of ' + str(n_primary) + ' primary mapped reads.'                                                    
    if write_dict:
        dict_out = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/' + number + '_assembly_dict.json'
        print 'WARNING: which of the 3 dictionaries to write?'
        #print 'Writing dictionary to ' + dict_out
        #with open(dict_out, 'w') as fp:          
        #    json.dump(assembly_dict, fp)
                
    if test_dict: # check construction by printing first entries to screen
        print 'Checking dictionary format (version 3).'
        x = itertools.islice(assembly_dict_3.iteritems(), 0, 4)
        for keyX, valueX in x:
            print keyX, valueX
        print 'Checking dictionary format (version 2).'
        y = itertools.islice(assembly_dict_2.iteritems(), 0, 4)
        for keyY, valueY in y:
            print keyY, valueY

#if save_reduced_dict:
#    R1_reduced_dict_out = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_reduced_dict.json'
#    print 'Writing dictionary to ' + R1_reduced_dict_out
#    with open('/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_reduced_dict', 'w') as fp:          
#        json.dump(R1, fp)

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