#!/usr/bin/python

# import python libraries
import argparse
import string as str
from string import Template, join
import subprocess
from subprocess import call, Popen, PIPE
from itertools import izip, izip_longest
import os as os
from os import linesep, path, R_OK, X_OK

# python main module
if __name__ == '__main__':
    
    # command line arguments
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-r', '--runVCFtools', help = 'Logical. Run VCFtools to get 012 output?', default = False)
    parser.add_argument('-v', '--vcftoolsPath', help = 'Path to vcftools executable.', default = 'vcftools')
    parser.add_argument('-i', '--inFile', help = 'Quality filtered VCF File.', default = None)
    parser.add_argument('-o', '--File012', help = 'Name for 012 file output.')
    parser.add_argument('-m', '--missingData', help = 'Character used to specify missing data.', default = '-1')
    parser.add_argument('-p', '--popFile', help = 'File containing populations of individuals in space delimited format on a single line, in an order corresponding to the order in the .indv file.', default = None)
    parser.add_argument('-d', '--phenotype', help = 'File containing phenotype of individuals in space delimited format on a single line, in an order corresponding to the order in the .indv file.', default = None)
    parser.add_argument('-n', '--numberSamples', help = 'Number of samples in file.')
    parser.add_argument('-a', '--adegenetFile', help = 'Name for adegenet-compatible output file.')
    opts = parser.parse_args()
    
    ##############################################
    ### Short function to check file inputs    ###
    ##############################################
    def checkFile(filename):
        '''
        return true if this is a file and is readable on the current filesystem
        '''
        try:
            if os.path.exists(os.path.abspath(filename)) and os.path.isfile(os.path.abspath(filename)) and os.access(os.path.abspath(filename), R_OK):
                return True
            fullPath = str.join(os.getcwd(), filename[1:])
            return os.path.exists(fullPath) and os.path.isfile(fullPath) and os.access(fullPath, R_OK)
        except IOError:
            return False
        
    ##############################################
    ### Parse arguments for inputs and outputs ###
    ##############################################
    vcfRun = opts.runVCFtools
    vcftoolsPath = opts.vcftoolsPath
    
    # are you working directly with a VCF file?
    if vcfRun:
        inVCF = opts.inFile
        checkFile(inVCF)    
        out012 = opts.File012
        # does a 012 output of the same name as you have assigned the new output already exist?
        if os.path.exists(out012):
            raise IOError('VCF 012 output file already exists. Choose a different name or skip VCF 012 generation.')
        
        # when VCF outputs your 012 file it will add extensions to the name you have chosen
        # these lines add extensions to your name strings so you can load the newly created 012 file
        out012_final = str.join([out012, '.012'], sep = '')
        if os.path.exists(out012_final):
            raise IOError('VCF 012 output file already exists. Choose a different name or skip VCF 012 generation.')
        
        out012_indv = str.join([out012, '.012.indv'], sep = '')
        if os.path.exists(out012_indv):
            raise IOError('VCF 012 output file already exists. Choose a different name or skip VCF 012 generation.')
        
    # if you are not working directly with a VCF file, you should already have a 012 file
    else:
        # you would have given the name with the extension as an input, so no extension needs to be added
        out012_final = opts.File012 
        checkFile(out012_final)
        # you will still need a file with the sample IDs, so tell python where to find the .indv based on the .012 output
        out012_indv = str.join([out012_final, '.indv'], sep = '')
        
    # load in the relevant phenotype and population info 
    pheno = opts.phenotype
    popFile = opts.popFile
    m = opts.missingData
    n = int(opts.numberSamples)
    
    # check that the output name is appropriate and the output doesn't already exist
    adegenet = opts.adegenetFile
    if adegenet.endswith('.snp'):
        if os.path.exists(adegenet):
            raise IOError('Parsed VCF 012 file for Adegenet already exists. Choose a different adegenet file name.')
    else:
        raise IOError('Output file name should use extension .snp for Adegenet compatibility.')
    
    # temp name for intermediates in the adegenet parsing process
    out012_tmp = str.join([out012_final, '.tmp'], sep = '')
    adegenet_tmp = str.join([adegenet, '.tmp'], sep = '')
    if os.path.exists(out012_tmp) or os.path.exists(adegenet_tmp):
        raise IOError('Temporary files from previous run still exist. Delete before proceeding.')
    
    ##############################################
    ### Convert the VCF file to 012 output     ###
    ##############################################
    
    if vcfRun:
        vcfTo012_Template = Template('%s --vcf $inFile --012 --out $outFile' % vcftoolsPath)
        commandLine = vcfTo012_Template.substitute(inFile = inVCF, outFile = out012)
        subprocess.call(commandLine, shell=True)
    
    ###############################################
    ### Replace missing data character with "-" ###
    ###############################################
    
    with open(out012_final, 'r') as iv:
        data = iv.readlines()
        with open(out012_tmp, 'a') as ov:
            for line in data:
                # replace m with - to recode missing data
                newLine = line.replace(m, '-')
                # expect tab-delimited 012 file; strip whitespace (this simultaneously takes the line and makes it a list
                newerLine = newLine.split('\t')
                newestLine = filter(None, newerLine)
                # expect first column to contain row numbers; remove
                newestLine.pop(0)
                # coerce from list back to a string
                finalLine = ''.join(newestLine)
                # write the line to a temporary file
                ov.write(finalLine)
    
    ###############################################
    ### interleave sample IDs and 012 genotpyes ###
    ###############################################
    
    # simultaneously open 3 files: 1 to write (adegenet_tmp) and two to read
    with open(adegenet_tmp, 'w') as a, open(out012_indv) as f1, open(out012_tmp) as f2:
        for line1, line2 in izip(f1, f2):
            # substitute lines from f1 and f2 into a template where they are separated by newlines
            newline = '>%s\n%s\n' % (line1.strip(), line2.strip())
            # write the line to a temporary file
            a.write(newline)     

    ##############################################################
    ### add a population ID for each sample to the file header ###
    ##############################################################
    
    # if there is a user-supplied population file
    if popFile:
        with open(popFile, 'r') as pf:
            finalPopline = pf.readlines()
            #print(finalPopline[0])
    # if no user-supplied population file, assign all samples to "pop1"
    else: 
        popline = ['pop1' for i in range(n)]
        finalPopline = ' '.join(popline)
    
    #######################################################
    ### add a phenotype for each sample to slot "other" ###
    #######################################################
    
    # if there is a user supplied phenotype file
    if pheno:
        with open(pheno, 'r') as pf:
            finalPheno = pf.readlines()
        header = '>>>> begin comments - do not remove this line <<<<\n>>>> end comments - do not remove this line <<<<\n>> population\n' + finalPopline + '\n>> ploidy\n2\n>> phenotype\n' + ''.join(finalPheno) + '\n'
    # if no user-supplied phenotype file, do nothing
    else:      
        header = '>>>> begin comments - do not remove this line <<<<\n>>>> end comments - do not remove this line <<<<\n>> population\n' + finalPopline + '\n>> ploidy\n2\n'

    ###################################
    ### write the final output file ###
    ###################################
    
    with open(adegenet_tmp, 'r') as a_original:
        data = a_original.read()
        with open(adegenet, 'w') as a_final:
            a_final.write(header)
            a_final.write(data)
            
