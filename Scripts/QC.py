#!/usr/bin/python

import os
import re
import subprocess

# run QC on the library-demultiplexed ddRADseq data with DBRs:

# check for barcodes, read 1
def barcode_r1(directory, out_name):
	files = os.listdir(directory)
	out_dir = directory + '/barcodes/'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	for f in files:
		r1 = re.findall(r'R1', f)
		if r1:
			barcode_check_call = 'zcat ' + directory + f + " | sed -n '2~4p' | cut -c 1-5 | sort | uniq -c | sort -nr -k 1 > " + out_dir + f + '_' + out_name
			subprocess.call(barcode_check_call, shell=True)

# check for cutsite, read 1
def cutsite_r1(directory, out_name):
	files = os.listdir(directory)
	out_dir = directory + '/cutsites/'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	for f in files:
		r1 = re.findall(r'R1', f)
		if r1:
			cutsite_check_call = 'zcat ' + directory + f + " | sed -n '2~4p' | cut -c 6-10 | sort | uniq -c | sort -nr -k 1 > " + out_dir + f + '_' + out_name
			subprocess.call(cutsite_check_call, shell=True)

def cutsite_r2(directory, out_name):
	files = os.listdir(directory)
	out_dir = directory + '/cutsites_read2/'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	for f in files:
		r1 = re.findall(r'R2', f)
		if r1:
			cutsite_check_call = 'zcat ' + directory + f + " | sed -n '2~4p' | cut -c 9-14 | sort | uniq -c | sort -nr -k 1 > " + out_dir + f + '_' + out_name
			subprocess.call(cutsite_check_call, shell=True)

# check for degeneracy, read 2
def degeneracy_r2(directory, out_name):
	files = os.listdir(directory)
	out_dir = directory + '/degeneracy'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	for f in files:
		r1 = re.findall(r'R2', f)
		if r1: 
			degeneracy_check_call = 'zcat ' + directory + f + " | sed -n 2~4p | cut -c 1-8 | sort | uniq -c | sort -nr -k 1 > " + out_dir + f + '_' + out_name
			subprocess.call(degeneracy_check_call, shell=True)
			
def run_FastQC(directory, out_name):
	files = os.listdir(directory)
	out_dir = directory + out_name
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	for f in files:
		fastqc_call = '/opt/software/FastQC/fastqc ' + directory + f + " -o " + out_dir
		subprocess.call(fastqc_call, shell=True)

barcode_r1('/home/antolinlab/Downloads/CWD_RADseq/', 'barcode_check')
cutsite_r1('/home/antolinlab/Downloads/CWD_RADseq/', 'cutsite_check')
cutsite_r2('/home/antolinlab/Downloads/CWD_RADseq/', 'cutsite_r2_check')
degeneracy_r2('/home/antolinlab/Downloads/CWD_RADseq/', 'degeneracy_check')
run_FastQC('/home/antolinlab/Downloads/CWD_RADseq/', 'fastQC')
