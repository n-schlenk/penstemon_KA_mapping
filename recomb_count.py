#! /usr/bin/env python3
# Created by N. Schlenk
# Extracts ID and recombination count from log file after looping through OrderMarkers2 4 times in each of the 8 linkage groups

import os
import re

recom_log=open('serial_test_35520326.log', 'r')
outfile=open('order5_recombination_log.txt', 'w')

order = 0
lg = 1
start = []

for n,line in enumerate(recom_log):
    outline = ''
    line_splat = line.split(' ')
    if line_splat[0] == 'number' and order < 4:
        order += 1
    elif order == 4:
        if line_splat[0].split('\t')[0] == 'Individual':
            outline = line_splat[0].split('\t')[2] + '\t' + str(line_splat[0].split('\t')[4]) + '\t' + str(lg) + '\n'
            outfile.write(outline)
        else:
            order = 0
            lg += 1
outfile.close()
