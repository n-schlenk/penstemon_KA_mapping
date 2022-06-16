#! /usr/bin/env python3
# written by N. Schlenk
# returns file with marker number, scaffold, position, and linkage group

import re
import os

vcf=open('[filtered VCF]', 'r')
map=open('[ordered map TXT]', 'r')  # same as map2
map2=open('[ordered map TXT]', 'r') # same as map
outfile=open('[coords TXT]', 'w')

LG_start={}
LG_end={}
pos={}
scaff={}
cm={}
lg={}
outline = ''


def split_map(map_file):
    header_lines = {}
    for n,line in enumerate(map_file):
        line_splat = line.split()
        if line_splat[0] == '#***':
            header_lines[int(line_splat[3])] = n
        last_line = n
    for i in range(1, len(header_lines.keys()) + 1):
        start_pos = header_lines[i] + 2
        if i < len(header_lines.keys()):
            end_pos = header_lines[i + 1] - 1
        elif i == len(header_lines.keys()):
            end_pos = n
        LG_start[i] = start_pos
        LG_end[i] = end_pos
    return(LG_start)
    return(LG_end)

split_map(map)

n = 1
for line in vcf:		
    if line[0] != '#':
        line_splat = line.split('\t')
        if line_splat[0] != 'CHR':
            pos[n] = line_splat[1]         		        	# populate pos dict as: keys = line (marker) number, values = position in bp
            scaff[n] = line_splat[0]	    			    # populate scaff dict as: keys = line (marker) number, values = scaffold
            n += 1

n = 0
for line in map2:
    if line[0] != '#':
        line_splat = line.split('\t')
        cm[int(line_splat[0])] = line_splat[1]  	    	# populate cm dict as: keys = marker number, values = position in cM
    for i in LG_start.keys():
        if n >= LG_start[i] and n <= LG_end[i]:
            lg[int(line_splat[0])] = i	            	    # populate lg dict as: keys = marker number, values = LG assigned in Lep-MAP3    
    n += 1

for marker_num in cm.keys():
    # outline = str(lg[marker_num]) + '\t' + scaff[marker_num] + '\t' + pos[marker_num] + '\t' + cm[marker_num] + '\n'
    outline = str(lg[marker_num]) + '\n'
    outfile.write(outline)

outfile.close()
