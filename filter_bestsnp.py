#! /usr/bin/env python3
# This program thins to best snp per radtag
# written by JKK
# modified by Carrie
# assumes you've already filtered for quality
# modified by N. Schlenk

MinMinor = 8            # min individuals that will have minor allele

vcf = open([VCF], 'r')
out2a = open([TXT output], 'w')


plants = 300
last_scaff = ''
bestsnp = ''
lastpos = 0
cp = 0
best_list = []

for line_idx, line in enumerate(vcf):
    cols = line.replace('\n', '').split('\t')
    scaff = cols[0]
    position = int(cols[1])
    ref_base = cols[3]
    alt_base = cols[4]
    if line_idx % 10000 == 0:
        print(scaff, position)
    mincc = 0
    if len(alt_base) == 1:
        datums = [0, 0, 0]
    for j in range(9, 9 + plants):
        if cols[j] != "./.":
            geno = cols[j].split(':')[0]
            if geno == '0/0':
                datums[0] += 1
            elif geno == '0/1':
                datums[1] += 1
            elif geno == '1/1':
                datums[2] += 1
            else:
                print('wtf2 ', geno)
    mincc = min(datums[0] + datums[1], datums[2], datums[1])
    if scaff != last_scaff or (position - lastpos) > 150:
        if cp >= MinMinor:
            out2a.write(bestsnp)
        cp = mincc
        bestsnp = scaff + '\t' + cols[1] + '\n'
        last_scaff = scaff
        lastpos = position
    elif mincc > cp and position != lastpos:
        cp = mincc
        bestsnp = scaff + '\t' + cols[1] + '\n'
if cp >= MinMinor:
    out2a.write(bestsnp)
out2a.close()
