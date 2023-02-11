#! usr/bin/env python3
# created by N. Schlenk
# input VCF to calculate genotype frequencies 

import re
import os
from numpy import mean
from numpy import std

invcf = open(<filepath to vcf>, 'r')

tot_RR = 0
tot_RA = 0
tot_AA = 0
for n,line in enumerate(invcf):
    if line[0] != '#' and line[0:3] != 'CHR' and line[0:3] == 'PGA':    # this is specific to our assembly, change to ignore columns beginning with name of assembly
        RR = 0
        RA = 0
        AA = 0
        line_splat = line.split('\t')
        for i in range(len(line_splat) - 2):
            id = line_splat[i]
            geno = id.split(':')
            if geno[0] == '0/0':
                RR += 1
            elif geno[0] == '0/1':
                RA += 1
            elif geno[0] == '1/1':
                AA += 1
        tot_RR += RR
        tot_RA += RA
        tot_AA += AA

tot_snps = tot_RR + tot_RA + tot_AA
tot_RR = round((100 * tot_RR / tot_snps), 3)
tot_RA = round((100 * tot_RA / tot_snps), 3)
tot_AA = round((100 * tot_AA / tot_snps), 3)
print('TOTAL genotype frequencies:')
print('AA : ' + str(tot_RR))
print('AB : ' + str(tot_RA))
print('BB : ' + str(tot_AA))


