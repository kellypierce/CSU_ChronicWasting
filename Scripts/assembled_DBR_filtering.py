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
import sqlite3
import numpy as np

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

test_dict=False
write_dict=False
save_reduced_dict=False

# load Dictionary

dict_in = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/dbr_dict'
seq_dict_in = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_dict'
  
with open(dict_in, 'r') as f:
    dbr = json.load(f)

with open(seq_dict_in, 'r') as s:
    R1 = json.load(s)
    
print 'Checking dictionary format (version 3).'
r = itertools.islice(R1.iteritems(), 0, 10)
for keyr, valuer in r:
    print keyr, valuer
print len(R1)

assembled_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Assembled'
for i in os.listdir(assembled_dir):
    
    # initialize an empty dictionary with each iteration of the for-loop
    assembly_dict_2 = {}
    assembly_dict_3 = defaultdict(list)
    
    # print some info to track progress
    path=os.path.join(assembled_dir, i)
    print 'Creating dictionary from ' + path
    
    # get the sample number for use later
    number = re.split('(\d+)', i)[1] # enclose the regex in parentheses to keep it in the output
    
    # open the sam file and process its contents
    with open(path, 'r') as inFile:
        for line in inFile:
            if not line.startswith("@"): # ignore the header lines
                fields = line.split("\t")
                
                # extract the info for the dictionary for each line
                QNAME = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', fields[0])[1] # FASTQ ID = QNAME column in sam file
                RNAME = fields[2] # ref sequence name -- where did the sequence map?
                QUAL = fields[10] # sequence quality score
                
                # extract the DBR corresponding to the QNAME for each row
                dbr_value = dbr.get(QNAME)[0] 
                
                # WE NEED TWO DICTIONARIES TO REPRESENT ALL THE RELATIONSHIP BETWEEN RNAME, QNAME, dbr_value, QUAL, AND count
                # build a dictionary with structure {DBR: (locus: count)}                    
                if RNAME not in assembly_dict_2:
                    assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1 # add the new DBR and its associated locus and count
                else:
                    if dbr_value not in assembly_dict_2.get(RNAME):
                        assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1
                    else:
                        assembly_dict_2[RNAME][dbr_value]=assembly_dict_2[RNAME][dbr_value]+1    
                
                # build a dictionary with structure {RNAME: {DBR:[[QNAME, QUAL]]}}        
                if RNAME not in assembly_dict_3:
                    assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL]]
                else:
                    if dbr_value not in assembly_dict_3.get(RNAME):
                        assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL]]
                    else:
                        assembly_dict_3[RNAME][dbr_value].append([QNAME, QUAL])

    # NOW THAT DICTIONARIES ARE MADE, REMOVE DUPLICATE SEQUENCES BASED ON DBR COUNTS
    # for each assembled locus, get the associated dbr_value and count
    for RNAME, value in assembly_dict_2.iteritems():
        print 'RNAME', RNAME
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
                    n_expected = 10
                    # if we count more DBR occurrences than we expect, we need to evaluate the quality of the sequences associated with that DBR
                    # calculating median QUAL from the full ASCII QUAL score is computationally intensive enough that we only want to do it when absolutely necessary -- only if we observe a DBR more often than expected
                    if count > n_expected:
                        #print 'count', count, 'n exp', n_expected
                        # the other dictionary contains the full quality information for each RNAME:DBR pair (this will be multiple entries of sequence IDs and qualities
                        qname_qual = assembly_dict_3[RNAME][dbr_value] #this is a list of lists: [[QNAME, QUAL], [QNAME, QUAL], ...]
                        ID_quals = {} # we'll make yet another dictionary to store the QNAME and the median QUAL
                        for i in qname_qual:
                            id=i[0]
                            ID_quals[id] = qual_mode(i[1], phred_dict)
                        # determine how many sequences will need to be removed and start a counter, t
                        n_remove = count - n_expected
                        t = 0
                        print 'count', count, 'expected', n_expected, 'remove', n_remove, 'tally', t
                        # with the full {ID: qual} dictionary available, we can now determine which sequences should be removed based on their score                        
                        while t < n_remove:
                            t += 1 
                            print ID_quals
                            print 'tally', t
                            print 'number of sequences for the given sample/locus/DBR combo', len(ID_quals)
                            if len(ID_quals)>0:
                                to_remove = min(ID_quals, key=ID_quals.get)
                                if to_remove == '8:1311:18936:26586':
                                    print '########################################### Goddamn duplicate'
                                #print 'removal ID', to_remove
                                #print 'sub dict value', ID_quals[to_remove]
                                #print 'full dict value', R1[to_remove]
                                # remove the sequences by their IDs from the master R1 sequence ID file
                                del ID_quals[to_remove]
                                del R1[to_remove] #8:1311:18936:26586
                                print len(R1)
                            else:
                                print "They're all gone!"
                            #try:
                            #    del R1[to_remove]
                            #except KeyError as ke:
                            #    print ke 
                            #    print '\n'.join(R1.keys())
                            #    print '\n' + to_remove + '\n'
                            
                            #if to_remove in R1:
                            #    print 'key present'
                            #    del R1[to_remove]
                            #else:
                            #    print '##################################'
                            
                            
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

if save_reduced_dict:
    R1_reduced_dict_out = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_reduced_dict.json'
    print 'Writing dictionary to ' + R1_reduced_dict_out
    with open('/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_reduced_dict', 'w') as fp:          
        json.dump(R1, fp)


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