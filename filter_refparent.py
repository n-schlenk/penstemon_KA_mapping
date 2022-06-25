#! /usr/bin/env python3
# written by N. Schlenk
# removes markers where the last ID (P2, in this case) is not homozygous for ref allele

import re
import os

vcf = open([VCF], 'r')
outvcf = open([filtered VCF], 'w')

nF2s = 298

for line in vcf:
    if line[0] == '#':
        outvcf.write(line)
    elif line[0] != '#':
        line_splat = line.split('\t')
        P1 = line_splat[298+9].split(':')[0]
        P2 = line_splat[298+10].split(':')[0]
        if P2 == '0/0':
            outvcf.write(line)

outvcf.close()
