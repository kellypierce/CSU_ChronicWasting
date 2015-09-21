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

test_dict=True
write_dict=False

# load Dictionary

dict_in = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/dbr_dict'
   
with open(dict_in, 'r') as f:
    dbr = json.load(f)
    #x = itertools.islice(dbr.iteritems(), 0, 4)
    #for key, value in x:
    #       print key, value
    #print dbr.get('8:2205:9340:77988')
    #if dbr.get('8:2205:9340:77988') == u'CTCACGGG':
    #    print 'match'
    #8:2205:9340:77988 [u'CTCACGGG']

assembled_dir = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/Assembled'
for i in os.listdir(assembled_dir):
    
    # initialize an empty dictionary with each iteration of the for-loop
    assembly_dict_0 = {}
    assembly_dict_1 = {}
    assembly_dict_2 = {}
    
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
                ID = re.split('(\d[:|_]\d+[:|_]\d+[:|_]\d+)', fields[0])[1]
                locus = fields[2]
                
                assembly_dict_0.setdefault(fields[2], {})[ID]=fields[10] # dict with ID and quality scores
                assembly_dict_1.setdefault(fields[2], {})[ID]=dbr.get(ID) # dict with ID and DBR
                dbr_value = dbr.get(ID)[0]
                #print locus, dbr_value
                
                # FUNCTIONAL: build a dictionary with structure {locus: {DBR: count}}
                if locus in assembly_dict_2:
                    if dbr_value not in assembly_dict_2.get(locus):
                        #print 'new'
                        assembly_dict_2.setdefault(locus, {})[dbr_value]=1 # dict with only DBR
                    else:
                        old_count=assembly_dict_2[locus][dbr_value]
                        #print old_count
                        assembly_dict_2[locus][dbr_value]=old_count+1
                else:
                    assembly_dict_2.setdefault(locus, {})
    
    if write_dict:
        dict_out = '/home/antolinlab/Desktop/CSU_ChronicWasting/PilotAnalysis/' + number + '_assembly_dict.json'
        print 'Writing dictionary to ' + dict_out
        with open(dict_out, 'w') as fp:          
            json.dump(assembly_dict, fp)
                
    if test_dict: # check construction by printing first entries to screen
        print 'Checking dictionary format.'
        x = itertools.islice(assembly_dict.iteritems(), 0, 4)
        for key, value in x:
            print key, value
    #c = Counter()
    #for key in assembly_dict.iterkeys():
    #    print key
    #    for subkey in key.iterkeys():
    #        c.update(set(subkey))
    #Counter(c.values())
        