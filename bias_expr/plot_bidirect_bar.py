import sys
import matplotlib as mpl
import matplotlib.pyplot as plt

from collections import defaultdict

mpl.use('agg')

'''
	input file format:
		pair	NaCl.0.3M	NaCl.0.8M	NaCl.0M	poplar
		JP19_08484@JP19_26196	0.6757231048097043	0.10742657783666666	1.4020984435713462	0.6205864104518777
		JP19_08483@JP19_26195	-1.1283240969755395	-1.0377196659491588	-0.8888933274526509	-1.7032040654153415
		JP19_08482@JP19_26194	0.01969619895860117	0.29034238688458813	0.5363142619977517	0.24115741744769623
		JP19_08481@JP19_26193	-2.577139795142226	-2.9496816275689017	-2.582442508742638	-2.915265425125856
		JP19_08480@JP19_26192	0.9509573770634813	0.8573546592663078	0.1938113202427457	0.47064196153900417
		JP19_08479@JP19_26191	0.01155104778554315	-0.21418463942971094	0.0	-0.011568809129889323
		JP19_08478@JP19_26190	-0.30079298952988165	-0.1439912796884557	-0.13233493784864855	-0.10778900412246067
		JP19_08477@JP19_26189	na	na	na	na
		JP19_08476@JP19_26188	0.9184842598132912	na	2.9928202237823456	3.200516548128313
'''

infile = sys.argv[1]
outfig = sys.argv[2]
log2_fold = 1
xlim = [-10, 10]
ylim = [0, 600]
ncols = 4
nrows = 4

data = defaultdict(list)
with open(infile, 'r' ) as FI:
	header = next(FI).rstrip().split('\t')
	for line in FI:
		entry = dict(zip(header[1:], line.rstrip().split('\t')[1:]))
		for t in entry:
			if entry[t] != 'na' and float(entry[t]) != 0:
				data[t].append(float(entry[t]))
			

fig = plt.figure(figsize=(10, 10))
gs = fig.add_gridspec(nrows=nrows, ncols=ncols)
fig.subplots_adjust(wspace=0.2, hspace=0.4)

axs = []
yticks_list = []
yup = []
print('tissue\tleft\tmiddle\tright')
for i, t in enumerate(header[1:]):
	r = i // nrows
	c = i % nrows
	ax = fig.add_subplot(gs[r, c])
	axs.append(ax)

	plot_data = data[t]
	#colors = ['lightblue']*(xlim[0]-1) + ['lightgray', 'lightgray'] + ['lightred']*(xlim[0]-1) 
	count, _, bar_container = ax.hist(plot_data, bins=list(range(xlim[0], xlim[1]+1)), color='lightgray')

	ax.set_title(t, fontsize=12)
	ax.tick_params(labelsize=8)
	ax.set_ylim(ylim)

	for i in range(xlim[1]-log2_fold):

		bar_container.patches[i].set(facecolor=(55/255,126/255,184/255,1), linewidth=None)

	for i in range(xlim[1]+log2_fold, xlim[1]*2):
		bar_container.patches[i].set(facecolor=(228/255,26/255,28/255,1), linewidth=None)

	for i in bar_container.patches:
		i.set(width=0.8)

	ltext = 'N=' + str(len([i for i in plot_data if i<=-log2_fold]))
	ax.text(xlim[0]/2, ylim[1]*0.75, ltext, ha='center')
	rtext = 'N=' + str(len([i for i in plot_data if i>=log2_fold]))
	ax.text(xlim[1]/2, ylim[1]*0.75, rtext, ha='center')

	print(t + '\t' + str(len([i for i in plot_data if i<=-log2_fold])) 
			+ '\t' + str(len([i for i in plot_data if -log2_fold<i<log2_fold]))
			+ '\t' + str(len([i for i in plot_data if i>=log2_fold]))
			)


	yup.append(max(count))
	yticks_list.append(ax.get_yticks())

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	if c != 0:
		ax.spines['left'].set_visible(False)


for idx in range(0, len(data), nrows):
	mi = yup.index(max(yup[idx:idx+nrows]))
	yticks = yticks_list[mi]

	axs[idx].set_yticks(yticks)
	axs[idx].set_yticklabels(map(int, yticks))

	for a in range(1, nrows):
		ax_i = idx  + a
		if ax_i < len(data):
			axs[ax_i].set_yticks([])


fig.savefig(outfig)
