#! /usr/bin/env 
# created by N. Schlenk
# generate final map coordinates for each scaffold

map_filepath = 'sc'             # prefix of each map file (to be followed by <number>.txt)
out_filepath = 'coords_sc'      # prefix of each coord file (to be followed by <number>.txt)
vcf_filepath = <filepath for VCF with SNP locations>
first_scaffold = 0
last_scaffold = 7

import re
import os
for i in range(first_scaffold, last_scaffold + 1):
    inmap=open(map_filepath + str(i) + '.txt', 'r')
    invcf=open(vcf_filepath, 'r')
    outfile=open(out_filepath+ str(i) + '.txt', 'w')
    bp = {}
    cm = {}
    n = 0
    for line in invcf:
        if line[0] != '#' and line[0:2] != 'KA' and line[0:3] != 'POS':
            pos = line.split('\t')[1]
            n += 1
            bp[str(n)] = pos
    for line in inmap:
        if line.split()[0] in bp.keys():
            x = bp[line.split()[0]]
            y = (line.split()[1])
            outline = str(x) + '\t' + str(y)
            outfile.write(outline + '\n')


