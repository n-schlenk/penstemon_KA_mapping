#! /usr/bin/env python3
# Created by N. Schlenk
# Extracts ID and recombination count from log file after looping through OrderMarkers2 [loop] times in each linkage group

import os
import re

recom_log=open('[log file]', 'r')
outfile=open('[recombination log]', 'w')

order = 0
lg = 1
loop = 4                                                            # number of times OrderMarkers2 was run
start = []

for n,line in enumerate(recom_log):
    outline = ''
    line_splat = line.split(' ')
    if line_splat[0] == 'number' and order < loop:
        order += 1
    elif order == loop:
        if line_splat[0].split('\t')[0] == 'Individual':
            outline = line_splat[0].split('\t')[2] + '\t' + str(line_splat[0].split('\t')[4]) + '\t' + str(lg) + '\n'
            outfile.write(outline)
        else:
            order = 0
            lg += 1
outfile.close()
