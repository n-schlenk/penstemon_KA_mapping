#! /usr/bin/env python3
# written by N. Schlenk
# creates csv of genotypes for mapping in Excel

import re
import os

vcf=open('KA_called_filter2.vcf', 'r')
outfile=open('map_all/geno_plot_all.csv', 'w')

nF2s = 298

for line in vcf:
    if line[0] != '#':
        line_splat = line.split('\t')
        out = str('')
        for i in range(9, nF2s+11):
            geno = line_splat[i].split(':')[0]
            if geno == './.':
                out += '-'
            elif geno == '1/1':
                out += 'A'
            elif geno == '0/0':
                out += 'R'
            elif geno == '1/0' or geno == '0/1':
                out += 'H'
            out += ',' 
        out += '\n'
        outfile.write(out)

outfile.close()
             
