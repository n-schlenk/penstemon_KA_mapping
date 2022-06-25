#! /usr/bin/env python3
# created by N. Schlenk
# returns p-value for Chi2 test (df = 1) for HWE per SNP


from scipy.stats import chi2
import re
import os

vcf = open('[VCF].vcf', 'r')
outfile = open('[HWE-filtered VCF].vcf', 'w')

nF2s = 298
minHWE = 0.001
minq = 0.3
maxq = 0.7

for line in vcf:
    n_alleles = 0
    R_allelecount = 0
    A_allelecount = 0
    R_obs = 0
    A_obs = 0
    H_obs = 0
    calls = 0
    if line[0] != '#':
        line_splat = line.split('\t')
        for i in range(9, nF2s+3):
            info = line_splat[i].split(':')
            gen = info[0]
            pl = info[1].split(',')
            ad = info[2].split(',')
            if gen != './.':
                calls += 1
#    optional replacement for using PL, not sure why it gives different results
#    if gen == '0/0':
#        R_obs += 1
#    elif gen == '1/0' or gen == '0/1':
#        H_obs += 1
#    elif gen == '1/1':
#        A_obs += 1
                if pl[0] == '0':
                    R_obs += 1
                if pl[1] == '0':
                    H_obs += 1
                if pl[2] == '0':
                    A_obs += 1
                R_allelecount += int(ad[0])
                A_allelecount += int(ad[1])
        n_alleles = R_allelecount + A_allelecount
        p = ( R_obs / calls ) + (0.5 * H_obs / calls)
        q = ( A_obs / calls ) + (0.5 * H_obs / calls)
        R_exp = p**2 * calls
        A_exp = q**2 * calls
        H_exp = 2 * p * q * calls
        chi2_stat = (((R_obs - R_exp)**2) / R_exp) + (((H_obs - H_exp)**2) / H_exp) + (((A_obs - A_exp)**2) / A_exp)
        chi2_sf = chi2.sf(chi2_stat, 1)
        if chi2_sf >= minHWE and (minq <= q <= maxq):
            outfile.write(line)
outfile.close()


