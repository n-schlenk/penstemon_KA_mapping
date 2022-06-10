#! /usr/bin/env python3
# written by N. Schlenk
# creates file with coordinates for plotting Lep-MAP3 output

import re
import os

vcf=open('./LG8/LG_8.vcf', 'r')                                         # vcf of scaffold
map=open('./LG8/ordered_lod3_LG8.txt', 'r')                             # ordered map file of scaffold
out=open('./LG8/coords_LG8.txt', 'a')                                   # output file

enum=1
bp={}
for line in vcf:                                                        # populate bp dictionary, keys = marker num, vals = bp
    if line[0] != '#':                                                  
        line_splat = line.split('\t')
        bp[str(enum)] = str(line_splat[1])                             
        enum += 1

cm={}
for line in map:                                                        # populate cm dictionary, keys = marker num, vals = cm
    if line[0] != '#':
        line_splat = line.split('\t')
        cm[str(line_splat[0])] = str(line_splat[1])

for marker_num in cm.keys():                                            # write output to file, base number followed by centimorgan value
    output = str(bp[marker_num]+'\t'+cm[marker_num]+'\n')
    out.write(output)
