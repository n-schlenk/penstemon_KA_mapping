#! /usr/bin/env python3
# written by N. Schlenk
# input map coordinates, output plot of 8 scaffolds

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

co0 = open('plot2/coords0_v2.txt', 'r')
co1 = open('plot2/coords1_v2.txt', 'r')
co2 = open('plot2/coords2_v2.txt', 'r')
co3 = open('plot2/coords3_v2.txt', 'r')
co4 = open('plot2/coords4_v2.txt', 'r')
co5 = open('plot2/coords5_v2.txt', 'r')
co6 = open('plot2/coords6_v2.txt', 'r')
co7 = open('plot2/coords7_v2.txt', 'r')

co = [co0, co1, co2, co3, co4, co5, co6, co7]

dict = {}
for i in range(8):
    dict.setdefault(i, ([], []))
    for n, scaffold in enumerate(co):
        if n == i:
            for line2 in scaffold:
                dict[i][0].append(int(line2.split()[0]))
                dict[i][1].append(float(line2.split()[1]))





fig, ax = plt.subplots(4, 2, sharex = 'col', sharey = 'row', constrained_layout = True)
ax[0, 0].scatter(dict[0][0], dict[0][1], alpha = 0.3, c = 'darkblue')
ax[1, 0].scatter(dict[1][0], dict[1][1], alpha = 0.3, c = 'darkblue')
ax[2, 0].scatter(dict[2][0], dict[2][1], alpha = 0.3, c = 'darkblue')
ax[3, 0].scatter(dict[3][0], dict[3][1], alpha = 0.3, c = 'darkblue')
ax[3, 0].set_xticks([0, 10000000, 20000000, 30000000, 40000000, 50000000, 60000000, 70000000])
ax[3, 0].set_xticklabels(['0', '10', '20', '30', '40', '50', '60', '70'])
ax[3, 0].xaxis.set_minor_locator(ticker.MultipleLocator(5000000))
ax[0, 1].scatter(dict[4][0], dict[4][1], alpha = 0.3, c = 'darkblue')
ax[1, 1].scatter(dict[5][0], dict[5][1], alpha = 0.3, c = 'darkblue')
ax[2, 1].scatter(dict[6][0], dict[6][1], alpha = 0.3, c = 'darkblue')
ax[3, 1].scatter(dict[7][0], dict[7][1], alpha = 0.3, c = 'darkblue')
ax[3, 1].set_xticks([0, 10000000, 20000000, 30000000, 40000000, 50000000])
ax[3, 1].set_xticklabels(['0', '10', '20', '30', '40', '50'])
ax[3, 1].xaxis.set_minor_locator(ticker.MultipleLocator(5000000))
ax[0, 0].yaxis.set_major_locator(ticker.MultipleLocator(50))
ax[1, 0].yaxis.set_major_locator(ticker.MultipleLocator(50))
ax[2, 0].set_yticks([0, 50, 100, 150], ['0', '', '100', ''])
ax[3, 0].yaxis.set_major_locator(ticker.MultipleLocator(50))
ax[0, 0].yaxis.set_minor_locator(ticker.MultipleLocator(25))
ax[1, 0].yaxis.set_minor_locator(ticker.MultipleLocator(25))
ax[2, 0].yaxis.set_minor_locator(ticker.MultipleLocator(25))
ax[3, 0].yaxis.set_minor_locator(ticker.MultipleLocator(25))
fig.supxlabel('Distance (Mbp)')
fig.supylabel('Distance (cM)')
fig.suptitle('Linkage Map by Chromosome')
ax[0, 0].set_title('Chromosome 1', fontsize = 10)
ax[1, 0].set_title('Chromosome 2', fontsize = 10)
ax[2, 0].set_title('Chromosome 3', fontsize = 10)
ax[3, 0].set_title('Chromosome 4', fontsize = 10)
ax[0, 1].set_title('Chromosome 5', fontsize = 10)
ax[1, 1].set_title('Chromosome 6', fontsize = 10)
ax[2, 1].set_title('Chromosome 7', fontsize = 10)
ax[3, 1].set_title('Chromosome 8', fontsize = 10)
plt.show()

