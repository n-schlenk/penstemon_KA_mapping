#! /usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import re

x1 = []
y1 = []
mn = []
ch_bp = []
mn2 = []
mn_cm = []

with open("lod9_ordered_sexavg.txt", "r") as ordered:           # open file with marker number and position in cM
    for line in ordered:
        x = re.findall("^[0-9]*\w", line)                       # add marker number to mn
        if x != []:
            mn.append(x)                       
with open("edit_fil.vcf", "r") as filtered:                     # open file with marker positions in db
    for line in filtered:
        y = line.split()
        ch_bp.append(y[0:2])                                    # add marker chromosome and position to ch_bp

for i in mn:                                                    # fix mn list to only contain integers
    a = str(i)
    b = len(a)
    mn2.append(a[2:b-2])

with open("lod9_ordered_sexavg.txt", "r") as ordered:           # open file with marker number and position in cM
    for line in ordered:
        z = line.split()
        mn_cm.append(z[0:2])                                    # add marker number and position to mn_cm

for i in mn_cm:                                                 # create x variable from marker position (bp)
    j = int(i[0])                                               # create y variable from marker position (M)
    k = ch_bp[j-1]
    l = float(i[1])
    m = int(l)
    x1.append(k[1])
    y1.append(m)
    
# fig, ax = plt.subplots()
# ax.set_title("Linkage Group 1")
# ax.set_xlabel("Marker Position (bp)")
# ax.set_ylabel("Marker Position (M)")
# ax.scatter(x1, y1)
# plt.show()

for index, x in enumerate(x1):
    print(x, y1[index])
