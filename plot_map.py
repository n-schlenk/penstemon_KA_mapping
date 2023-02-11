#! /usr/bin/env python3
# written by N. Schlenk
# input map coordinates, output plot of 8 scaffolds
# you can plot another number of scaffolds but you will need to adjust the plot code

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

first_scaffold = 0
last_scaffold = 7
filepath_prefix = 'coords'
filepath_suffix = '.txt'
graph_color = 'darkblue'
titles = ['LG0', 'LG1', 'LG2', 'LG3', 'LG4', 'LG5', 'LG6', 'LG7']

dict = {}
for i in range(first_scaffold, last_scaffold + 1):
    dict.setdefault(i, ([],[]))
    with open(filepath_prefix + str(i) + filepath_suffix, 'r') as infile:
        for line in infile:
            dict[i][0].append(int(line.split()[0]))
            dict[i][1].append(float(line.split()[1]))

fig, ax = plt.subplots(4, 2, sharex = 'col', sharey = 'row', constrained_layout = True)
for i in range(0, 4):
    ax[i, 0].scatter(dict[i][0], dict[i][1], alpha = 0.3, c = graph_color)
    ax[i, 0].yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax[i, 0].yaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax[i, 0].set_title(titles[i], fontsize = 10)
for i in range(4, 8):
    ax[i-4, 1].scatter(dict[i][0], dict[i][1], alpha = 0.3, c = graph_color)
    ax[i-4, 1].yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax[i-4, 1].yaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax[i-4, 1].set_title(titles[i], fontsize = 10)
for i in range(0, 2):
    ax[3, i].set_xticks([0, 10000000, 20000000, 30000000, 40000000, 50000000, 60000000, 70000000])
    ax[3, i].set_xticklabels(['0', '10', '20', '30', '40', '50', '60', '70'])
    ax[3, i].xaxis.set_minor_locator(ticker.MultipleLocator(5000000))
fig.supxlabel('Distance (Mbp)')
fig.supylabel('Distance (cM)')
fig.suptitle('Linkage Map by Chromosome')
plt.show()


