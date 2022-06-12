#! /usr/bin/env python3
# created by N. Schlenk
# returns p-value for Chi2 test (df = 1) for HWE per SNP


from scipy.stats import chi2
import re
import os

vcf = open('KA_called_filter2.vcf', 'r')
outfile = open('KA_called_filter4.txt', 'w')


nF2s = 298

p_dict = {}
count = 0
for line in vcf:
    n_alleles = 0
    R_allelecount = 0
    A_allelecount = 0
    R_obs = 0
    A_obs = 0
    H_obs = 0
    if line[0] != '#':
        line_splat = line.split('\t')
        for i in range(9, nF2s+11):
            gen = line_splat[i].split(':')[0]
            if gen == '0/0':
                R_allelecount += 2
                R_obs += 1
            elif gen == '1/0' or gen == '0/1':
                R_allelecount += 1
                A_allelecount += 1
                H_obs += 1
            elif gen == '1/1':
                A_allelecount += 2
                A_obs += 1
        n_alleles = R_allelecount + A_allelecount
        if n_alleles != 0:
            R_exp = ((R_allelecount / n_alleles)**2) * n_alleles
            A_exp = ((A_allelecount / n_alleles)**2) * n_alleles
            H_exp = 2 * (R_allelecount / n_alleles) * (A_allelecount / n_alleles)
            if R_exp != 0 and A_exp != 0 and H_exp != 0:
                chi2_stat = (((R_obs - R_exp)**2) / R_exp) + (((H_obs - H_exp)**2) / H_exp) + (((A_obs - A_exp)**2) / A_exp)
                p_val = 1 - (chi2.sf(chi2_stat, 1))
                count += 1
                p_dict[count] = p_val
print(p_dict.values())

vcf.close()
