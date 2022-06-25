#! /usr/bin/env python3
# created by JKK and Carrie
# modified by N. Schlenk

in2a = open([output findbestsnp], 'r')
vcf = open([VCF], 'r')
outfile = open([filtered VCF], 'w')

bestlist = []
for line in in2a:
    cols = line.replace('\n', '').split('\t')
    bestlist.append([cols[0], cols[1]])
for line in vcf:
    cols = line.replace('\n', '').split('\t')
    tig = cols[0]
    pos = cols[1]
    for x in bestlist:
        if tig == x[0] and pos == x[1]:
            outfile.write(line)
outfile.close()
