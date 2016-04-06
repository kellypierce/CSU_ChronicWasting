#!/usr/bin/python

import argparse
import string as str
from string import Template, join
import subprocess
from subprocess import call, Popen, PIPE
from itertools import izip, izip_longest
import os as os
from os import linesep, path, R_OK, X_OK

if __name__ == '__main__':
        
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-r', '--runVCFtools', help = 'Logical. Run VCFtools to get 012 output?', default = False)
    parser.add_argument('-v', '--vcftoolsPath', help = 'Path to vcftools executable.', default = 'vcftools')
    parser.add_argument('-i', '--inFile', help = 'Quality filtered VCF File.', default = None)
    parser.add_argument('-o', '--outFile', help = 'Name for 012 file output.')
    parser.add_argument('-p', '--popFile', help = 'File containing populations of individuals in space delimited format on a single line, in an order corresponding to the order in the .indv file.', default = None)
    parser.add_argument('-d', '--phenotype', help = 'File containing phenotype of individuals in space delimited format on a single line, in an order corresponding to the order in the .indv file.', default = None)
    parser.add_argument('-n', '--numberSamples', help = 'Number of samples in file.')
    parser.add_argument('-a', '--adegenetFile', help = 'Name for adegenet-compatible output file.')
    opts = parser.parse_args()
    
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
    
    # parse arguments for inputs and outputs
    vcfRun = opts.runVCFtools
    vcftoolsPath = opts.vcftoolsPath
    if vcfRun:
        inVCF = opts.inFile
        checkFile(inVCF)    
        out012 = opts.outFile
        if os.path.exists(out012):
            raise IOError('VCF 012 output file already exists. Choose a different name or skip VCF 012 generation.')
        out012_final = str.join([out012, '.012'], sep = '')
        out012_indv = str.join([out012, '.012.indv'], sep = '')
    else:
        out012_final = opts.outFile
        checkFile(out012_final)
        out012_indv = str.join([out012_final, '.indv'], sep = '')
    pheno = opts.phenotype
    popFile = opts.popFile
    n = int(opts.numberSamples)
    adegenet = opts.adegenetFile
    if os.path.exists(adegenet):
        raise IOError('Parsed VCF 012 file for Adegenet already exists. Choose a different adegenet file name.')
    
    # temp name for intermediates in the adegenet parsing process
    out012_tmp = str.join([out012_final, '.tmp'], sep = '')
    adegenet_tmp = str.join([adegenet, '.tmp'], sep = '')
    if os.path.exists(out012_tmp) or os.path.exists(adegenet_tmp):
        raise IOError('Temporary files from previous run still exist. Delete before proceeding.')
    
    # Convert the VCF file to 012 output if required
    if vcfRun:
        vcfTo012_Template = Template('%s --vcf $inFile --012 --out $outFile' % vcftoolsPath)
        commandLine = vcfTo012_Template.substitute(inFile = inVCF, outFile = out012)
        subprocess.call(commandLine, shell=True)
    
    # Replace "-1" with "-" as the character for missing data
    with open(out012_final, 'r') as iv:
        data = iv.readlines()
        with open(out012_tmp, 'a') as ov:
            for line in data:
                # replace -1 with - to recode missing data
                #missingRecode = re.sub('-1', '-', iv)
                newLine = line.replace('-1', '-')
                newerLine = newLine.split('\t')
                newestLine = filter(None, newerLine)
                newestLine.pop(0)
                finalLine = ''.join(newestLine)
                #newLine = re.sub(' ', '', missingRecode)
                ov.write(finalLine)
    
    # interleave sample IDs and 012 genotpyes
    with open(adegenet_tmp, 'w') as a, open(out012_indv) as f1, open(out012_tmp) as f2:
        for line1, line2 in izip(f1, f2):
            #a.write(">{}\n{}\n".format(line1.rstrip(), line2.rstrip()))  
            newline = '>%s\n%s\n' % (line1.strip(), line2.strip())
            a.write(newline)     

    # untested:
    if popFile:
        with open(popFile, 'r') as pf:
            finalPopline = pf.readlines()
            #print(finalPopline[0])
    else: 
        popline = ['pop1' for i in range(n)]
        finalPopline = ' '.join(popline)
    if pheno:
        with open(popFile, 'r') as pf:
            finalPheno = pf.readlines()
        header = '>>>> begin comments - do not remove this line <<<<\n>>>> end comments - do not remove this line <<<<\n>> population\n' + ''.join(finalPopline[0]) + '\n>> ploidy\n2\n>> other\n' + ''.join(finalPheno) + '\n'
    else:      
        header = '>>>> begin comments - do not remove this line <<<<\n>>>> end comments - do not remove this line <<<<\n>> population\n' + ''.join(finalPopline[0]) + '\n>> ploidy\n2\n'

    with open(adegenet_tmp, 'r') as a_original:
        data = a_original.read()
        with open(adegenet, 'w') as a_final:
            a_final.write(header)
            a_final.write(data)
            
    #os.remove("*.tmp")