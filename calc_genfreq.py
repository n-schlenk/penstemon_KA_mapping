#! usr/bin/env python3
# created by N. Schlenk
# input VCF to calculate genotype frequencies 

import re
import os
from numpy import mean
from numpy import std

invcf = open('filteredvcf.vcf', 'r')

tot_RR = []
tot_RA = []
tot_AA = []
for n,line in enumerate(invcf):
    RR = 0
    RA = 0
    AA = 0
    if n > 37:
	line_splat = line.split('\t')
        for i in range(187):
            id = line_splat[i]
            geno = id.split(':')
            if geno[0] == '0/0':
                RR += 1
            elif geno[0] == '0/1':
                RA += 1
            elif geno[0] == '1/1':
                AA += 1
        tot_RR.append(RR)
        tot_RA.append(RA)
        tot_AA.append(AA)
print(mean(tot_RR), std(tot_RR))
print(mean(tot_RA), std(tot_RA))
print(mean(tot_AA), std(tot_AA))

