#! /usr/bin/env python3

from scipy.stats import chi2
import matplotlib.pyplot as plt
import numpy as np

vcf = open('vcfs/no_hwe_v2_postfil.vcf', 'r')
nF2s = 238
minHWE = 0.01

rr = 0
ra = 0
aa = 0

for line in vcf:
    if line[0] != '#':
        freq = [0, 0, 0]
        for i in line.split()[9:9 + nF2s]:
            geno = i.split(':')[0]
            if geno == '0/0':
                freq[0] += 1
            if geno == '0/1':
                freq[1] += 1
            if geno == '1/1':
                freq[2] += 1
        p = float((2 * freq[0]) + freq[1])/(2 * (freq[0] + freq[1] +  freq[2]))
        q = float((2 * freq[2]) + freq[1])/(2 * (freq[0] + freq[1] +  freq[2]))
        exp_r = (p * p) * (freq[0] + freq[1] + freq[2])
        exp_h = (p * q) * (freq[0] + freq[1] + freq[2])
        exp_a = (q * q) * (freq[0] + freq[1] + freq[2])
        r_term = ((freq[0] - exp_r) * (freq[0] - exp_r)) / exp_r
        h_term = ((freq[1] - exp_h) * (freq[1] - exp_h)) / exp_h
        a_term = ((freq[2] - exp_a) * (freq[2] - exp_a)) / exp_a
        rr += freq[0]
        ra += freq[1]
        aa += freq[2]
#        chisq_term = r_term + h_term + a_term
#        pvalue = (1 - chi2.cdf(chisq_term, 1))
#        if pvalue >= minHWE:
#            print(freq)

tot_calls = rr + ra + aa
hwe_hom = tot_calls / 4
hwe_het = tot_calls / 2
hwe_data = [0.25, 0.50, 0.25]
obs_data = [rr/tot_calls, ra/tot_calls, aa/tot_calls]

w = 0.3
plt.bar(np.arange(len(hwe_data)), hwe_data, width = w, color = 'thistle')
plt.bar(np.arange(len(obs_data))+ w, obs_data, width = w, color = 'indigo')
plt.xticks([0.15, 1.15, 2.15], ['K/K', 'K/A', 'A/A'])
plt.legend(['Expected', 'Observed'])
plt.xlabel('Genotype')
plt.ylabel('SNP Frequency')
plt.show()
