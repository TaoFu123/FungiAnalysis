import argparse
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('agg')

'''
	输入文件格式:
	x, y, color

	example:
	FJChr01	0.8182	#48A462
	FJChr02	0.0195	#999999
	FJChr03	0.0341	#999999
	FJChr04	0.4498	#48A462
	FJChr05	0.0237	#999999
	FJChr06	0.0324	#999999
'''

parser = argparse.ArgumentParser()
parser.add_argument('data')
parser.add_argument('-o', '--outfile', default='bar.svg')
parser.add_argument('-d', '--width', default='0.8', type=float)
parser.add_argument('--title', default='')
parser.add_argument('--figsize', default='6.8,4.8', help='width,height, default:6.8,4.8')
args = parser.parse_args()


bars = []
with open(args.data, 'r') as DT:
	for line in DT:
		if line.startswith('#'):
			continue

		entry = line.rstrip().split()
		x, y, color = entry
		y = float(y)
		bars.append([x, y, color])

figsize = list(map(float, args.figsize.split(',')))
fig, ax = plt.subplots(figsize=figsize)
for i, b in enumerate(bars, 1):
	ax.bar(i, b[1], width=args.width, color=b[2], linewidth=0, alpha=0.8)

ax.set_xticks(range(1, len(bars)+1))
ax.set_xticklabels([i[0] for i in bars], fontsize=8, rotation=90)
ax.set_ylabel('Coverage')
ax.set_ylim([0, 1])
ax.set_title(args.title)
#fig.
fig.savefig(args.outfile)
