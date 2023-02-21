#! /usr/bin/env 
# created by N. Schlenk
# generate final map coordinates for each scaffold



#! /usr/bin/env
# created by N. Schlenk
# generate final map coordinates
import re
import os

first_scaffold = 0
last_scaffold = 7
scaffold_index = 12
map_prefix = 'ordered_scaffold_'
coords_prefix = 'coords/coords'
invcf=open('big_vcfs/minminor1_pc.vcf', 'r')
dict = {}
for line in invcf:
    if line[0] != '#' and line[0:2] != 'KA' and line[0:3] != 'POS':
        pos = line.split('\t')[1]
        chr = line[scaffold_index]
        dict.setdefault(chr, [])
        dict[chr].append(pos)

for i in range(first_scaffold, last_scaffold + 1):
    n = 0
    inmap=open(map_prefix + str(i) + '.txt', 'r')
    outfile=open(coords_prefix + str(i) + '.txt', 'w')
    for line in inmap:
        if line[0] != '#' and len(line.split()) > 1:
            print(line)
            outline = str(dict[str(i)][n]) + '\t' + str(line.split()[1])
            outfile.write(outline + '\n')
            n += 1



