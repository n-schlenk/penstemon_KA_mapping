#! /usr/bin/env 
# created by N. Schlenk
# generate final map coordinates for a single scaffold

import re
import os

inmap=open('sc0_order.txt', 'r')
invcf=open('scaffold_0.vcf', 'r')
outfile=open('coords0_v2.txt', 'w')

bp = {}
cm = {}

n = 0
for line in invcf:
    if line[0] != '#':
        pos = line.split('\t')[1]
        bp[str(n)] = pos
        n += 1

for line in inmap:
    if line.split()[0] in bp.keys():
        x = bp[line.split()[0]]
        y = (line.split()[1])
        outline = str(x) + '\t' + str(y)
        outfile.write(outline + '\n')

