# SPDX-FileCopyrightText: 2023 CERN
# SPDX-License-Identifier: Apache-2.0

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from cycler import cycler
import math

if len(sys.argv) < 5:
	print("Usage: python3 plot_ratio.py output_file x_label y_label data.csv [data_2.csv data_3.csv ...]")
	exit()

output_file = sys.argv[1]
x_label = sys.argv[2]
y_label = sys.argv[3]
data_files = sys.argv[4:]

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
#Uncomment this line for custom color cycle
#plt.rc('axes', prop_cycle=(cycler('color', ['c', 'orange', 'r', 'g', 'yellow', 'violet'])))

f1 = data_files[0]
f2 = data_files[1]
data1 = pd.read_csv(f1)
data2 = pd.read_csv(f2)

ratio = data1/data2

ratio_mean = ratio.mean()

ratio_error = np.sqrt(data1.std()**2 + data2.std()**2) * ratio_mean

print(ratio_error)

width=0.2

plt.figure(figsize=(18, 10))

x = np.arange(len(ratio.columns))
#Draw grid below other figures
plt.gca().set_axisbelow(True)
plt.grid(True, axis='y', color='black', linestyle='dotted')
#Plot the data in a bar chart
plt.errorbar(x=x, y=ratio_mean, yerr=ratio_error, linewidth=1, marker="s", elinewidth=1, label="Ratio")
plt.xticks(x, ratio.columns)
plt.ylabel(y_label)
plt.legend()

#plt.show()
plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.5)


