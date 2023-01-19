#! /usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import numpy as np

awk0 = open('csv_stk_txt/awk_output_sc0.txt', 'r')
co0 = open('coords/sc0_v4.txt', 'r')
vcf0 = open('vcfs/scaffold_0.vcf', 'r')
sc0 = open('fastas/scaffold0.fasta', 'r')
awk2 = open('csv_stk_txt/awk_output_sc2.txt', 'r')
co2 = open('coords/sc2_v4.txt', 'r')
vcf2 = open('vcfs/scaffold_2.vcf', 'r')
sc2 = open('fastas/scaffold2.fasta', 'r')
awk3 = open('csv_stk_txt/awk_output_sc3.txt', 'r')
co3 = open('coords/sc3_v4.txt', 'r')
vcf3 = open('vcfs/scaffold_3.vcf', 'r')
sc3 = open('fastas/scaffold3.fasta', 'r')
awk7 = open('csv_stk_txt/awk_output_sc7.txt', 'r')
co7 = open('coords/sc7_v4.txt', 'r')
vcf7 = open('vcfs/scaffold_7.vcf', 'r')
sc7 = open('fastas/scaffold7.fasta', 'r')

awks = [awk0, awk2]
cos = [co0, co2]
vcfs = [vcf0, vcf2]
scs = [sc0, sc2]

titles = ['Linkage Group 0', 'Linkage Group 2']
colors = ['darkblue', 'mediumvioletred']

bin = 1000000

recombs = {}
for n, awk_file in enumerate(awks):
    recombs.setdefault(n, []) 
    geno_dict = {}
    for line in awk_file:
        if line[0] != '#':
            geno = line.split()[4:len(line.split())]
            for id in range(0, len(geno), 2):
                plant_index = id / 2
                geno_dict.setdefault(plant_index, [])
                g1 = geno[id]
                g2 = geno[id + 1]
                if g1 == '1' and g2 == '1':
                    gt = 'A'
                if g1 == '2' and g2 == '2':
                    gt = 'B'
                else:
                    gt = 'H'
                geno_dict[plant_index].append(gt)
    for id in geno_dict.keys():
        for i in range(3, len(geno_dict[id]) - 3):
            if geno_dict[id][i] != geno_dict[id][i+1]:
                recombs[n].append(i)                                # these are indexed and need to be given physical bp locations

rec_bp = {}
snp_den = {}
for n, vcf in enumerate(vcfs):
    rec_bp.setdefault(n, [])
    snp_den.setdefault(n, [])
    bp = []
    for line in vcf:
        if line[0] == 'P':
            snp_den[n].append(int(line.split()[1]))
    for num, pos in enumerate(snp_den[n]):
        if num in recombs[n]:
            rec_bp[n].append(pos)                                   # now, they are no longer indexed

x_dict = {}
y_dict = {}
for n, coord_file in enumerate(cos):
    x_dict.setdefault(n, [])
    y_dict.setdefault(n, [])
    for line in coord_file:
        x = line.split()[0]
        y = line.split()[1].replace('\n','')
        x_dict[n].append(int(x)/bin)
        y_dict[n].append(float(y))


gc_dict = {}
for num, fasta in enumerate(scs):
    gc_dict.setdefault(num, [])
    gc = 0
    tot = 0
    x = 0
    for line in fasta:
        if line[0] != '>' and line.replace('\n','') != '':
            for nuc in line.replace('\n',''):
                x += 1
                if nuc == 'G' or nuc == 'C':
                    gc += 1
                    tot += 1
                elif nuc == 'A' or nuc == 'T':
                    tot += 1
                if x % bin == 0 and x != 0:
                    gc_dict[num].append(gc / tot)
                    gc = 0
                    tot = 0
    gc_dict[num].append(round((gc / tot), 4))

##### binning step #####
snp_heat = {}
rec_heat = {}
snp_hist = {}
rec_hist = {}

for i in range(0, 70000000, bin):
    for chr, val in snp_den.items():
        snp_heat.setdefault(chr, [])
        snp_hist.setdefault(chr, [])
        x = 0
        for snp in val:
            if i < snp <= i + bin:
                snp_hist[chr].append(i / bin)
                x += 1   
        snp_heat[chr].append(x)
    for chr, val in rec_bp.items():
        rec_heat.setdefault(chr, [])
        rec_hist.setdefault(chr, [])
        x = 0
        for snp in val:
            if i < snp <= i + bin:
                rec_hist[chr].append(i)
                x += 1
        rec_heat[chr].append(x)

patches = []
steps = int(10000000/bin)
fig, ax = plt.subplots(4, sharex = True, gridspec_kw = {'height_ratios': [1, 0.25, 0.1, 0.1]})
for chr in x_dict.keys():
    ax[0].scatter(x_dict[chr], y_dict[chr], c = colors[chr], alpha = 0.2, s = 1)
    ax[0].yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax[0].set_ylabel('Genetic Distance (cM)')
    ax[0].set_title('Metacentric Centromeres (Predicted)', fontsize = 12)
    ax[1].hist(snp_hist[chr], histtype = 'step', color = colors[chr], alpha = 0.7)
    ax[1].set_xticks(*[range(0, (7 * steps), steps)], ['0', '10', '20', '30', '40', '50', '60'])
    ax[1].set_xlim(0, 7 * steps)
    ax[1].set_ylabel('SNP Count')
    patches.append(mpatches.Patch(color = colors[chr], label = titles[chr]))
ax[2].imshow([gc_dict[0]], cmap = 'plasma', aspect = 'auto', vmin = 0.3000, vmax = 0.4000)    
ax[2].set_yticks([0], ['LG0'])
ax[3].imshow([gc_dict[1]], cmap = 'plasma', aspect = 'auto', vmin = 0.3000, vmax = 0.4000)    
ax[3].set_yticks([0], ['LG2'])
ax[3].set_xlabel('Position (Mbp)')
ax[0].legend(handles = patches)
plt.suptitle('Linkage Maps and SNP Coverage', fontsize = 16)
plt.show()
