#!/usr/bin/PYTHON

import integrated_denovo_pipeline as pipeline
import QC as qc
import os
import re

# run the FASTX toolkit quality filters on Read 2 only (lots of data; will be slow)
# get the degeneracy stats for the quality filtered Read 2 files
pipeline.iterative_FASTQ_quality_filter(directory = '~/CWD_RADseq/', out_dir = '/qual_filtered_R2_for_DBR_distr/', out_name = 'qual_filtered_30.fastq.gz', q = 30, p = 50, read = 'R2')
qc.degeneracy_r2(directory = '~/CWD_RADseq/', out_name = 'qual_filtered_30_degeneracy_check')
