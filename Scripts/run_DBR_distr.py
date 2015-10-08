#!/usr/bin/PYTHON

print 'running quality filtering and DBR counting for DBR distribution generation'

import integrated_denovo_pipeline as pipeline
import QC as qc
import os
import re

def iterative_FASTQC(directory, out_name, q, p):
    files = os.listdir(directory)
    out_dir = directory + '/qual_filtered_R2_for_DBR_distr/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print out_dir
    for f in files:
        r2 = re.findall(r'R2', f)
        if r2: 
            in_file = directory + f
            out_file = out_dir + f + out_name
            pipeline.FASTQ_quality_filter(in_file, out_file, q, p)

iterative_FASTQC('/home/antolinlab/Downloads/CWD_RADseq/', 'qual_filtered_30.fastq.gz', 30, 50)
qc.degeneracy_r2('/home/antolinlab/Downloads/CWD_RADseq/', 'qual_filtered_30_degeneracy_check')
