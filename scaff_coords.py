#! /usr/bin/env 
# created by N. Schlenk
# generate final map coordinates for a single scaffold

import re
import os

inmap=open('[scaffold map].txt', 'r')
invcf=open('[un-LMfiltered scaffold].vcf', 'r')
outfile=open('coords.txt', 'w')

bp = {}
cm = {}

n = 0
for line in invcf:
    if line[0] != '#':
        pos = line.split('\t')[1]
        n += 1
	bp[str(n)] = pos

for line in inmap:
    if line.split()[0] in bp.keys():
        x = bp[line.split()[0]]
        y = (line.split()[1])
        outline = str(x) + '\t' + str(y)
        outfile.write(outline + '\n')


