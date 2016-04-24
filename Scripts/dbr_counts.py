#!/usr/bin/python

from string import Template, join
import re
import os
import multiprocessing as mp
from subprocess import call, Popen, PIPE
import subprocess
from assembled_DBR_filtering import *

parallel_DBR_count(in_dir = '/home/pierce/CWD_RADseq/raw_new/qual_filtered/',
                   dbr_start = -9,
                   dbr_stop = -2,
                   save = '/home/pierce/CWD_RADseq/raw_new/degeneracy_checks/',
                   saveType = 'text')


