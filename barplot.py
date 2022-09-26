#! /usr/bin/env python3
# created by N. Schlenk
# words
# use values from calc_genfreq.py to plot exp vs obs gen freq

import re
import os
import matplotlib.pyplot as plt
import numpy as np


vals, sd = (32.057, 58.091, 38.541), (11.518, 14.624, 13.331)
geno = ['K/K', 'K/A', 'A/A']
hwe = (25, 50, 25)

x_axis = np.arange(3)
fig, ax = plt.subplots()
obs = ax.bar(x_axis - 0.2, vals, yerr = sd, width = 0.4, label = 'Observed', color = 'mediumvioletred', alpha = 0.8)
exp = ax.bar(x_axis + 0.2, hwe, label = 'Expected', width = 0.4, color = 'darkgrey')
plt.xticks(x_axis, geno)
plt.legend()
plt.show()
