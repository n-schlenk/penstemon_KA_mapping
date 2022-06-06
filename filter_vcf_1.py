#! /usr/bin/env python3
# code by Carrie, edited by Noelle
import os
import re

invcf = open('KA_called_nodup.vcf', 'r')                                # define input file, output file, number F2s and min individuals to call a genotype at a locus
outfile = open('output_filter1.vcf', 'w')
nF2s = 298
minIndiv = 50

def skipHeader(line):                                                   # function 'skipHeader'
    cols = line.replace('\n','').split('\t')                            # split line based on tabs
    if len(cols) < 2:                                                   # if there are less than 2 of them, return 'header'
        return 'header'
    elif cols[0] == '#CHROM':                                           # if the first one is the chromosome header, return 'header'
        return 'header'
    else:
        return cols                                                     # otherwise, everything else stays the same


def extractVCFfields(sampleData):                                       # function 'extractVCFfields' takes in only sample data, not entire line
    fields = sampleData.split(":")                                      # split GT, PL, AD
    dataLength = len(fields)
    if (dataLength > 2):
        alleleDepths = fields[2].split(',')                             # split allele depth section based on comma 
        phreds = fields[1].split(',')                                   # split phred score section based on comma
        return [alleleDepths, phreds]
    else:
        return 'missing'

def findMapQual(infoField):                                             # function 'findMapQual' takes in the information column and splits it to get the MQscore
    info = infoField.split(';')                                         # split field based on semicolon
    listpos = len(info) - 1                                             # list position starts at the second-to-last chunk and moves left until it finds MQ
    MQscore = 0
    while listpos > 0: # start searching for MQ from end of info field
        if (info[listpos].split('='))[0] == 'MQ': # found it!
            MQscore = float((info[listpos].split('='))[1])
            listpos = 0
        else:
            listpos -= 1 # keep searching
    return MQscore                                                      # return MQscore as float variable
###################################################################

nsamples = nF2s + 2
minMQ = 30
for line in invcf:
    cols = skipHeader(line)                                              # use 'skipHeader' function to return only informative lines
    if cols != 'header':   
        MQscore = findMapQual(cols[7])                                   # MQscore is MQ score for each individual locus
        if MQscore >= minMQ:                                             # if MQscore is greater than minimum allowed value...
            altBase = cols[4]                  
            if len(altBase) > 1:                                         # and if there is only one alt allele...
                pass
            else:                                                        # count number of calls made based on presence of 'missing' applied via 'extractVCFfields' function
                calls = 0
                for j in range(9, 9+nsamples):
                    annot = extractVCFfields(cols[j])
                    if annot != "missing":
                        calls += 1
                amph = cols[nsamples-2 + 9]
                kunth = cols[nsamples-1 + 9]
                amphfields = amph.split(':')                             # split genotype data to extract genotype and eliminate loci where either parent is a het
                kunthfields = kunth.split(':')
                if amphfields[0] != './.' and kunthfields[0] != './.':
                    if amphfields[0] == '0/1' or kunthfields[0] == '0/1':
                        pass
                    elif amphfields[0] != kunthfields[0]:                # if the remaining locus is different between parents (and passed filters), keep it
                        if calls >= minIndiv:
                            outfile.write(line)                          # write out lines that passed filters into outfile
                        else:
                            pass
outfile.close()
