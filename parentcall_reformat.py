#! /usr/bin/env python3
# created by N Schlenk


pc = open(<parentcalled VCF>, 'r')      # parentcalled VCF
ref = open(<pre-parentcall VCF>, 'r')   # VCF from pre-parentcalling
outfile = open(<outfile path>, 'w')     # outfile

snp_list = []
for line in pc:
    if line[0:3] == 'PGA':
        snp_list.append(line.split('\t')[1])

for line in ref:
    if line[0:3] == 'PGA' and line.split('\t')[1] in snp_list:
        outfile.write(line)



