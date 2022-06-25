#! /usr/bin/env python3
# created by Carrie on 10/7/17
# edited by N. Schlenk

from numpy import median
from numpy import mean
from numpy import sum
from scipy.stats import chi2

in1 = open('[VCF]', 'r')
outfile = open('[HWE-filtered VCF]', 'w')

nF2s = 298
nsamples = nF2s + 2
minq = 0.3
maxq = 0.7
minHWE = 0.0001

for line in in1:
    cols = line.replace('\n', '').split('\t')
    if len(cols) < 2:
        pass
    elif cols[0] == '#CHROM':
        pass
    else:
        tig = cols[0]
        pos = cols[1]
        totREF = 0
        totALT = 0
        calls = 0  
        reads = []
        RRcount = 0
        RAcount = 0
        AAcount = 0
        for j in range(9, 9+nsamples):
            fields = cols[j].split(':')    
            if fields[0] != './.':
                calls += 1
                alleleDepths = fields[2].split(',')
                refDepth = int(alleleDepths[0])
                altDepth = int(alleleDepths[1])
                reads.append(float(refDepth + altDepth))
                phreds = fields[1].split(',')
                if phreds[0] == '0':
                    RRcount += 1
                elif phreds[1] == '0':
                    RAcount += 1
                elif phreds[2] == '0':
                    AAcount += 1
        meanRPI = mean(reads)
        sumRPI = sum(reads)
        freqRR = float(RRcount) / calls
        freqRA = float(RAcount) / calls
        freqAA = float(AAcount) / calls  
        freqR = freqRR + 0.5 * freqRA
        freqA = freqAA + 0.5 * freqRA
        expRR = calls * freqR * freqR
        expRA = calls * 2 * freqA * freqR
        expAA = calls * freqA * freqA
        if expRR == 0:
            termRR = 0
        else:
            termRR = (((RRcount - expRR) ** 2) / expRR)
        if expRA == 0:
            termRA = 0
        else:
            termRA = (((RAcount - expRA) ** 2) / expRA)
        if expAA == 0:
            termRA = 0
        else:
            termAA = (((AAcount - expAA) ** 2) / expAA)
        chisq = termRR + termRA + termAA
        pvalue = 1.0 - chi2.cdf(chisq, 1)
        if pvalue >= minHWE and (minq <= freqR <= maxq):
            outfile.write(line)
outfile.close()
