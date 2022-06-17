#! /usr/bin/env python3
# written by N. Schlenk
# takes log file with ID in first column and recombination event count in second column
# returns VCF with IDs removed that had more than 4 rec. events in a LG
# also removes IDs that passed the first filter but total rec. count was in top 10% of all 'passed' rec. counts


import os
import re


invcf=open('filter6/KA_called_filter6.vcf', 'r')
log=open('order5_recombination_log.txt', 'r')
outvcf = open('KA_called_filter8.vcf', 'w')

bad = []
good = {}

for line in log:
    line_splat = line.split('\t')
    if int(line_splat[1]) > 4:
        bad.append(line_splat[0])
    elif int(line_splat[1]) < 4:
        if line_splat[0] not in good.keys():
            good[line_splat[0]] = int(line_splat[1])
        else:
            good[line_splat[0]] += int(line_splat[1])

perc90 = round(len(good.keys()) - (len(good.keys()) / 10))
for n,val in enumerate(sorted(good.values())):
    if n < perc90:
        max = val

for i in good.keys():
    if good[i] >= max:
        bad.append(i)

for n,line in enumerate(invcf):
    if n < 36:
        outvcf.write(line)
    else:
        line_splat = line.split('\t')
        if line_splat[0] == '#CHROM':
            for m,id in enumerate(line_splat):
                if id in bad:
                    bad.append(m)
    if n >= 36:
        line_splat = line.replace('\n', '').split('\t')
        outline = ''
        for m,i in enumerate(line_splat):
            if m in bad:
                pass
            else:
                best = line_splat[m]
                outline += best + '\t'
        outline += '\n'
        outvcf.write(outline)     

outvcf.close()
