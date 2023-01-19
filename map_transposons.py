#! /usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

co0 = open('coords/no_hwe_v2_sc0_coords_edit.txt', 'r')
co1 = open('coords/no_hwe_v2_sc1_coords_edit.txt', 'r')
co2 = open('coords/no_hwe_v2_sc2_coords_edit.txt', 'r')
co3 = open('coords/no_hwe_v2_sc3_coords_edit.txt', 'r')
co4 = open('coords/no_hwe_v2_sc4_coords_edit.txt', 'r')
co5 = open('coords/no_hwe_v2_sc5_coords_edit.txt', 'r')
co6 = open('coords/no_hwe_v2_sc6_coords_edit.txt', 'r')
co7 = open('coords/no_hwe_v2_sc7_coords_edit.txt', 'r')
infile = open('csv_stk_txt/kunthii-families.stk', 'r')
fasta0 = open('fastas/scaffold0.fasta', 'r')
fasta1 = open('fastas/scaffold1_v2.fasta', 'r')
fasta2 = open('fastas/scaffold2_v2.fasta', 'r')
fasta3 = open('fastas/scaffold3_v2.fasta', 'r')
fasta4 = open('fastas/scaffold4_v2.fasta', 'r')
fasta5 = open('fastas/scaffold5_v2.fasta', 'r')
fasta6 = open('fastas/scaffold6_v2.fasta', 'r')
fasta7 = open('fastas/scaffold7_v2.fasta', 'r')
copia_out = open('plots/copia_heatmap_perc.txt', 'w')
gypsy_out = open('plots/gypsy_heatmap_perc.txt', 'w')
gc_out = open('plots/gc_heatmap.txt', 'w')

coords = [co0, co1, co2, co3, co4, co5, co6, co7]
fastas = [fasta0, fasta1, fasta2, fasta3, fasta4, fasta5, fasta6, fasta7]
scaffold_id = [0, 1, 2, 3, 4, 5, 6, 7]
colors = ['darkblue', 'mediumvioletred', 'indigo', 'hotpink', 'green', 'coral', 'crimson', 'yellow']
targets = ['Ty1-Copia', 'Gypsy']

by_scaff_repeats = {}
x = 0
for line in infile:
    if line[0] == '#' and line.split()[1] == 'TP':
        x = 0
        for target in targets:
            if target in line.split()[2].split(';'):
                x = 1
                t = target
    if line[0:3] == 'PGA' and x == 1:
        if int(line[12]) in scaffold_id:
            by_scaff_repeats.setdefault(int(line[12]), {})
            by_scaff_repeats[int(line[12])].setdefault(t, [])
            by_scaff_repeats[int(line[12])][t].append(int(line.split()[0].split(':')[1].split('-')[0]) / 1000000)

print('repeats done')

seqs = {}
for n, id in enumerate(scaffold_id):
    seqs.setdefault(id, '')
    fasta = fastas[n]
    for line in fasta:
        if line[0] != '>':
            seqs[id] += line.replace('\n','')

print('sequence done')

gc_content = {}
for id in seqs.keys():
    gc = 0
    at = 0
    gc_content.setdefault(id, [])
    for n, nuc in enumerate(seqs[id]):
        if nuc == 'G' or nuc == 'C' or nuc == 'g' or nuc == 'c':
            gc += 1
        elif nuc == 'A' or nuc == 'T' or nuc == 'a' or nuc == 't':
            at += 1
        if n != 0 and n % 1000000 == 0:
            gc_content[id].append(gc / (at + gc))
            gc = 0
            at = 0
    gc_content[id].append(gc / (at + gc))

print('GC done')

x_dict = {}
y_dict = {}
for n, id in enumerate(scaffold_id):
    x_dict.setdefault(id, [])
    y_dict.setdefault(id, [])
    for line in coords[n]:
        x_dict[n].append(int(line.split()[0]) / 1000000)
        y_dict[n].append(float(line.split()[1]))

print('coords done')

heatmap = {}
for id in scaffold_id:
    heatmap.setdefault(id, {})
    k = 0
    for rep_ele in by_scaff_repeats[id].keys():
        temp = []
        for i in range(len(gc_content[id])):
            j = 0
            for loc in by_scaff_repeats[id][rep_ele]:
                if i < loc < i + 1:
                    j += 1
                    k += 1
            temp.append(j)
        for count in temp:
            heatmap[id].setdefault(rep_ele, [])
            heatmap[id][rep_ele].append(count / k)

print('heatmap done')

a = (len(targets) + 1) * len(scaffold_id) + 1
heights = [5]
for i in range(len(scaffold_id)):
    heights.append(0.2)
    heights.append(0.2)
    heights.append(0.2)







fig, ax = plt.subplots(a, sharex = True, gridspec_kw = {'height_ratios': heights})
for n, scaff in enumerate(scaffold_id):
    ax[0].scatter(x_dict[scaff], y_dict[scaff], c = colors[n], s = 0.5, alpha = 0.3)
    ax[0].set_title('Linkage Map for Chromosomes 1-8')
    for m, rep in enumerate(heatmap[scaff].keys()):
        ax[(3 * n) + 1 + m].imshow([heatmap[scaff][rep]], aspect = 'auto', cmap = 'plasma', vmin = 0.01, vmax = 0.06)
    ax[(3 * n) + 3].imshow([gc_content[scaff]], aspect = 'auto', cmap = 'viridis', vmax = 0.4)
    ax[(3 * n) + 1].set_yticks([])
    ax[(3 * n) + 2].set_yticks([])
    ax[(3 * n) + 3].set_yticks([])
    ax[(3 * n) + 3].set_xticks([0, 10, 20, 30, 40, 50, 60, 70])
plt.show()
