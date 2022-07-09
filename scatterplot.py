#! /usr/bin/env python3
# created by N. Schlenk
# produces custom scatterplot from coordinate data
# redundant code --- will fix later

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import os
import re

x0 = [0]
y0 = [0]
x1 = [0]
y1 = [0]
x2 = [0]
y2 = [0]
x3 = [0]
y3 = [0]
x4 = [0]
y4 = [0]
x5 = [0]
y5 = [0]
x6 = [0]
y6 = [0]
x7 = [0]
y7 = [0]

file0=open('coords0.txt', 'r')
file1=open('coords1.txt', 'r')
file2=open('coords2.txt', 'r')
file3=open('coords3.txt', 'r')
file4=open('coords4.txt', 'r')
file5=open('coords5.txt', 'r')
file6=open('coords6.txt', 'r')
file7=open('coords7.txt', 'r')

for line in file0:
    bp = int(line.split()[0])
    cm = float(line.split()[1])
    x0.append(bp)
    y0.append(cm)
for line in file1:
    bp = int(line.split()[0])
    cm = float(line.split()[1])
    x1.append(bp)
    y1.append(cm)
for line in file2:
    bp = int(line.split()[0])
    cm = float(line.split()[1])
    x2.append(bp)
    y2.append(cm)
for line in file3:
    bp = int(line.split()[0])
    cm = float(line.split()[1])
    x3.append(bp)
    y3.append(cm)
for line in file4:
    bp = int(line.split()[0])
    cm = float(line.split()[1])
    x4.append(bp)
    y4.append(cm)
for line in file5:
    bp = int(line.split()[0])
    cm = float(line.split()[1])
    x5.append(bp)
    y5.append(cm)
for line in file6:
    bp = int(line.split()[0])
    cm = float(line.split()[1])
    x6.append(bp)
    y6.append(cm)
for line in file7:
    bp = int(line.split()[0])
    cm = float(line.split()[1])
    x7.append(bp)
    y7.append(cm)

print(pearsonr(x0, y0))
print(pearsonr(x1, y1))
print(pearsonr(x2, y2))
print(pearsonr(x3, y3))
print(pearsonr(x4, y4))
print(pearsonr(x5, y5))
print(pearsonr(x6, y6))
print(pearsonr(x7, y7))
#plt.scatter(x1, y1, color ='slateblue', alpha=0.5)
#plt.scatter(x5, y5, color ='slateblue', alpha=0.5)
#plt.scatter(x6, y6, color ='slateblue', alpha=0.5)
#plt.scatter(x7, y7, color ='slateblue', alpha=0.5)
plt.scatter(x0, y0, color ='slateblue', alpha=0.5)
plt.scatter(x2, y2, color ='mediumvioletred', alpha=0.5)
plt.scatter(x3, y3, color ='plum', alpha=0.5)
plt.scatter(x4, y4, color ='pink', alpha=0.5)
plt.xlabel('Physical Position (bp)')
plt.ylabel('Genetic Distance (cM)')
plt.title('Genetic Map of 4 Chromosomes')
# plt.legend(['Scaffold 0', 'Scaffold 2', 'Scaffold 3', 'Scaffold 4'])
plt.show()

