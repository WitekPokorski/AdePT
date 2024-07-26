# SPDX-FileCopyrightText: 2023 CERN
# SPDX-License-Identifier: Apache-2.0

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from cycler import cycler


def plot_lines(output_file, x_label, y_label, data_files):

	prop_cycle = plt.rcParams['axes.prop_cycle']
	colors = prop_cycle.by_key()['color']
	#Uncomment this line for custom color cycle
	#plt.rc('axes', prop_cycle=(cycler('color', ['c', 'orange', 'r', 'g', 'yellow', 'violet'])))

	width=0.2

	plt.figure(figsize=(18, 10))

	for i in range(len(data_files)):
		file = data_files[i]
		data = pd.read_csv(file)

		#Sort the columns in ascending order
		data = data.reindex(data.mean().sort_values().index, axis=1)
		
		#Get the mean of each timing
		means = data.mean()
		#Get the standard error for each timing.
		errors = [np.std(data[column]) for column in data.columns]

		x = [2**i for i in np.arange(len(data.columns))]
		
		#Draw grid below other figures
		plt.gca().set_axisbelow(True)
		plt.grid(True, axis='y', color='black', linestyle='dotted')
		
		#Plot the data
		plt.errorbar(x=x, y=means, yerr=errors, label=file, linewidth=1, marker="s", elinewidth=1)
		plt.xticks(x, data.columns)

		plt.ylim(bottom=0, top=1)
		start, end = plt.gca().get_ylim()
		plt.gca().yaxis.set_ticks(np.arange(start, end, 0.1))

		plt.ylabel(y_label, fontsize=18)
		plt.xlabel(x_label, fontsize=18)
		plt.legend(fontsize=14)
	
    #Plot hardware limit
	#plt.axvline(x=32, color="black", linestyle="--", linewidth=0.8)
	#plt.text(31, 0.42, "Hardware Limit", rotation=90, fontsize=14)

	#Execution info
	plt.text(x[-1]/2, 0.05, 
		  "Machine:   2 x AMD EPYC 7552 (96 cores) + Nvidia A100\n\
Test:      1 Run, 1-96 CPU threads, 4 TTBar per thread\n\
Geometry:  CMS Run 2",
		  fontname="monospace",
		  fontsize=14,
		  bbox={"boxstyle" : "Square",
		  		"facecolor" : "white"}
		  )

	#plt.show()
	plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.5)


if __name__ == '__main__':
	if len(sys.argv) < 5:
		print("Usage: python3 plot_lines.py output_file x_label y_label data.csv [data_2.csv data_3.csv ...]")
		exit()

	output_file = sys.argv[1]
	x_label = sys.argv[2]
	y_label = sys.argv[3]
	data_files = sys.argv[4:]

	plot_lines(output_file, x_label, y_label, data_files)