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

# load Dictionary

dict_in = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/dbr_dict'
seq_dict_in = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/R1_dict.json'
  
with open(dict_in, 'r') as f:
    dbr = json.load(f)

with open(seq_dict_in, 'r') as s:
    R1 = json.load(s)

assembled_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Assembled'
for i in os.listdir(assembled_dir):
    
    # initialize an empty dictionary with each iteration of the for-loop
    assembly_dict_0 = {}
    assembly_dict_1 = {}
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
                RNAME = fields[2] # ref sequence name ('locus', if you will)
                QUAL = fields[10] # sequence quality score
                
                dbr_value = dbr.get(QNAME)[0] 
                
                
                #assembly_dict_0.setdefault(RNAME, {})[QNAME]=QUAL # dict with ID and quality scores
                #assembly_dict_1.setdefault(RNAME, {})[QNAME]=dbr_value # dict with ID and DBR
                #assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[QNAME, QUAL] # dict with ID and DBR
    
                #print locus, dbr_value
                
                # FUNCTIONAL: build a dictionary with structure {DBR: (locus: count)}                    
                if RNAME not in assembly_dict_2:
                    assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1 # add the new DBR and its associated locus and count
                else:
                    if dbr_value not in assembly_dict_2.get(RNAME):
                        assembly_dict_2.setdefault(RNAME, {})[dbr_value]=1
                    else:
                        assembly_dict_2[RNAME][dbr_value]=assembly_dict_2[RNAME][dbr_value]+1    
                        
                if RNAME not in assembly_dict_3:
                    assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL]]
                else:
                    if dbr_value not in assembly_dict_3.get(RNAME):
                        assembly_dict_3.setdefault(RNAME, {})[dbr_value]=[[QNAME, QUAL]]
                    else:
                        assembly_dict_3[RNAME][dbr_value].append([QNAME, QUAL])

    # now each sample has 2 dictionaries describing the relationships between RNAME, QNAME, dbr_value, QUAL, and count
    for RNAME, value in assembly_dict_2.iteritems():
        print 'RNAME', RNAME
        if RNAME != '*':
            for subvalue in value.iteritems():
                #print 'Subvalue', subvalue
                dbr_value = subvalue[0] 
                count = subvalue[1]
                if count > 1:
                    # we don't have a likelihood function yet, so draw from a random distribution to see if we need to calculate quality score:
                    n_expected = np.random.uniform(2, 10)
                    if count > n_expected:
                        print 'count', count, 'n exp', n_expected
                        qname_qual = assembly_dict_3[RNAME][dbr_value]
                        ID_quals = {}
                        for i in qname_qual:
                            id=i[0]
                            ID_quals[id] = qual_mode(i[1], phred_dict)
                        n_remove = count - n_expected
                        t = 0
                        while t <= n_remove:
                            #print ID_quals
                            to_remove = min(ID_quals, key=ID_quals.get)
                            #print to_remove
                            del seq_dict_in[to_remove]
                            t += 1
                            
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



'''
# we only need to calculate the quality score for some loci
                countQUAL = Counter(QUAL)
                freqQUAL = countQUAL.most_common() # counts of all unique ASCII characters
                modeQUAL = countQUAL.most_common(1) # most frequent occurrence
                listQUAL = list(QUAL) # split the string into a list
                list_intQUAL = []
                for q in listQUAL:
                    list_intQUAL.append(phred_dict[q])
                intQUAL = phred_dict[modeQUAL[0][0]]
                medQUAL = np.median(list_intQUAL)
                print 'mode', intQUAL, 'median', medQUAL
'''