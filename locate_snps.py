#! /usr/bin/env python3

infile = open('extra_super_filtered.vcf', 'r')

scaffolds = [0, 0, 0, 0, 0, 0, 0, 0]
for line in infile:
    if line[0:3] == 'PGA':
        for i in range(8):
            if line.split('_')[1].replace('scaffold','') == str(i):
                scaffolds[i] += 1
print(scaffolds)
