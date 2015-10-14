from string import Template
import logging as logging
from logging import debug, critical, error, info
from os import linesep, path
import os as os
import sys as sys

def configureLogging(verbose = False):
    '''
    setup the logger
    '''
    logger = logging.getLogger()
    logger.handlers = []

    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    
    consoleLogHandler = logging.StreamHandler()

    # do not need colors if logging to a non-shell
    if sys.stdout.isatty():
        consoleLogHandler.setFormatter(logging.Formatter("\033[93m%(filename)s:%(lineno)s\033[0m: %(message)s"))
    else:
        consoleLogHandler.setFormatter(logging.Formatter("%(message)s"))

    logger.addHandler(consoleLogHandler)

def checkFile(filename):
    '''
    return true if this is a file and is readable on the current filesystem
    '''
    try:
        if os.path.exists(os.path.abspath(filename)) and os.path.isfile(os.path.abspath(filename)) and os.access(os.path.abspath(filename), R_OK):
            return True
        fullPath = join(os.getcwd(), filename[1:])
        return os.path.exists(fullPath) and os.path.isfile(fullPath) and os.path.access(fullPath, R_OK)
    except IOError:
        return False

def checkExe(filename):
    '''
    return true if this is an executable file on the current filesystem
    '''
    if not isinstance(filename, str):
        raise TypeError("need a string, got a %s" % type(filename))
    return (exists(filename) and isfile(filename) and access(filename, X_OK))

################################## GLOBALS ####################################
# run configureLogging first
configureLogging(False)
qualityFilter = '/Users/Kelly/Downloads/bin/fastq_quality_filter'
fqfStdinTemplate = Template('%s -q $q -p $p -Q33 -z -o $output' % qualityFilter)
fqfFileTemplate = Template('%s -q $q -p $p -Q33 -z -i $input -o $output' % qualityFilter)
###############################################################################

# function to iteratively call "fastq_quality_filter"
def FASTQ_quality_filter(fq_in, fq_out, q, p, qualityFilter = qualityFilter):
    '''
    fq_in: full path to the in-file, a zipped text file with sequence and quality data
    fq_out: full path to the output directory
    q additional param to fastq_quality filter
    p additional param to fastq_quality filter
    '''
    # note that the above is a PEP0257 docstring, it's standardized
    
    if not checkFile(fq_in):
        raise FileNotFoundError("where is the input file: %s" % fq_in)
    
    if not fq_in.endswith("gz"):
        warn("prefer to pass compressed files, consider gzipping the inputs")

    if not checkExe(qualityFilter):
        raise Exception("could not find %s in the filesystem for execution, is the environment setup correctly?" % qualityFilter)
	    
    info("Quality filtering with FASTX Toolkit.")
    
    # files seem to be handled differently when they are compressed
    # I gather "fastq_quality_filter" doesn't really expect a compressed file, but the below works okay in the shell
    
    # chain a gz decompressor thread to fqf
    if fq_in.endswith('gz'):
        commandLine = fqfStdinTemplate.substitute(q = q,
                                                  p = p,
                                                  output = fq_out)
        
        debug(commandLine)
        
        zcatProcess = Popen('zcat %s' % fq_in,
                            shell = True,
                            stdout = subprocess.PIPE)
        
        fqfProcess = Popen(commandLine,
                           shell = True,
                           stdin = zcatProcess.stdout)
        
    #if 'gz' in fq_in:
    #    fqc_call = "zcat " + fq_in + " | fastq_quality_filter -q " + str(q) + " -p " + str(p) + "-Q33 -z -o " + fq_out 
    
    # the uncompressed files seem to behave better, but I'd like to avoid uncompressing all the files at this stage
    # if that's not possible, then I can go ahead and uncompress
    # but is there a way to pass the compressed file to fastq_quality_filter in a pipe-style buffered fashion?
    
    else:
        #fqc_call = "fastq_quality_filter -q " + str(q) + " -p " + str(p) + " -Q33 -z -i " + fq_in + " -o " + fq_out
        commandLine = fqfFileTemplate.substitute(q = q, 
                                                 p = p,
                                                 output = fq_out,
                                                 input = fq_in)
        
        debug(commandLine)
        
        fqfProcess = Popen(commandLine,
                           shell = True)
    
    # make this synchronous, it will wait until fqf finishes before leaving the function
    # alternatively you can pass in a [], which you can .append(fqfProcess) and have all of them running 
    # in parallel, probably you want to also pass a counter so that it can self-limit how many 
    # processes are running at a time, this would be done via a check, followed by time.sleep(#) to tell it 
    # to go take a nap and check back later   
    fqfProcess.wait() 
    
    #print fqc_call
    #subprocess.call(fqc_call, shell=True)

    #return
    
    # This structure of passing some arguments to a python function that strings them together
    # and makes subprocess calls to external software is a common theme in the pipeline I've developed.
    # Many of the little helper-functions I've made actually iterate through all the files in a directory;
    # this one just happens to be a bit simpler... this function actually gets called by a separate
    # function I wrote that does the iteration (otherwise this function seems silly! Can't remember
    # why I structured it this way, but I'm sure there was a reason... But I'll apply the lessons learned here 
    # to the other functions and hopefully improve their performance a bit. Thanks for your help!
    
FASTQ_quality_filter('/Users/Kelly/Desktop/CSU_ChronicWasting/Library12_S65_L008_R2_001.fastq.gz', '/Users/Kelly/Desktop/CSU_ChronicWasting/test_filter/', q = 30, p = 50)