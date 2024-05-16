import numpy as np
from collections import defaultdict
from sklearn.decomposition import PCA

import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('agg')

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

n_components = 5

expr_dict = defaultdict(dict)
gids = set()
with open('expr.data.form.xls', 'r') as FI:
	header = next(FI).rstrip().split('\t')
	for line in FI:
		entry = dict(zip(header, line.rstrip().split('\t')))
		sid = entry['sid']
		gid = entry['gid']
		#if 'poplar' not in sid:
		expr_dict[sid][gid] = np.log2(float(entry['rTPM'])+1)
		gids.add(gid)

gids = sorted(gids)

filt_expr = defaultdict(dict)
filt_gids = set()
for g in gids:
	ex = []
	for t in expr_dict:
		ex.append(expr_dict[t][g])

	if any(ex):
		filt_gids.add(g)
		for t in expr_dict:
			filt_expr[t][g] = expr_dict[t][g]

expr_dict = filt_expr
gids = sorted(filt_gids)

expr_data = [[expr_dict[si].get(gi, 0) for gi in gids] for si in sorted(expr_dict)]

	

with open('form.dataFrame.xls', 'w') as DF:
	DF.write('spec' + '\t' + '\t'.join(gids) + '\n')
	for si, ex in zip(sorted(expr_dict), expr_data):
		DF.write(si + '\t' + '\t'.join(map(str, ex)) + '\n')

#归一化
norm_data1 = []
for i in expr_data:
	min_i = np.min(i)
	max_i = np.max(i)
	dist = max_i - min_i
	norm_data1.append([(j-min_i)/dist for j in i])


#标准化
norm_data2 = []
for i in norm_data1:
	m = np.mean(i)
	std = np.std(i)
	norm_data2.append([(j-m)/std for j in i])

norm_data2 = np.array(norm_data2)
print(norm_data2.shape)
pca = PCA(n_components=n_components)
pca.fit(norm_data2)
print(pca.explained_variance_ratio_)
x_reduction = pca.transform(norm_data2)

markers = ['o', '^', 's', 'x', '*', 'p', 'd']
colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628"]

t = set()
d = set()
for sid in expr_dict:
	t.add(sid.split('-')[0])
	d.add(sid.split('-')[1])

marker_dict = {i:j for i, j in zip(sorted(t), markers)}
color_dict = {i:j for i, j in zip(sorted(d), colors)}

plot_marker = [marker_dict[sid.split('-')[0]] for sid in sorted(expr_dict)]
plot_color = [color_dict[sid.split('-')[1]] for sid in sorted(expr_dict)]

if n_components == 2:
	fig, ax = plt.subplots()
	ax.scatter(x_reduction[:,0], x_reduction[:,1], c=plot_color, marker=plot_marker)
	fig.savefig('pca.pdf')

else:
	fig = plt.figure(figsize=(12, 12))
	gs = fig.add_gridspec(nrows=n_components-1, ncols=n_components-1)
	fig.subplots_adjust(wspace=0.05, hspace=0.05)
	
	axs = []
	for i in range(0, n_components-1):
		for j in range(1, n_components):
			r = j-1
			c = i
			ax = fig.add_subplot(gs[r, c])
			axs.append(ax)

			ax.set_xticks([])
			ax.set_yticks([])

			ax.xaxis.set_label_position('top') 
			if r == 0:
				ax.set_xlabel('PC{}({}%)'.format(i+1,
												round(pca.explained_variance_ratio_[i]*100, 2)))
			if c == 0:
				ax.set_ylabel('PC{}({}%)'.format(j+1,
												round(pca.explained_variance_ratio_[j]*100, 2)))
			
			if i < j:
				plot_data = defaultdict(list)
				for x, y, c, m in zip(x_reduction[:,i], x_reduction[:,j], plot_color, plot_marker):
					plot_data[m].append([x, y, c])

				for m in plot_data:
					ax.scatter([i[0] for i in plot_data[m]], [i[1] for i in plot_data[m]], c=[i[2] for i in plot_data[m]], marker=m, alpha=0.8)

			if r == 0 and c == n_components-2:
				handles = []
				for t in marker_dict:
					m = marker_dict[t]
					handles.append(ax.scatter([], [], c='grey', marker=m, label=t))
				
				for d in color_dict:
					c = color_dict[d]
					handles.append(ax.scatter([], [], c=c, marker='o', label=d))

				ax.legend(handles=handles, ncol=2, loc='upper right', frameon=False)

	
	#axs[-1].legend()

fig.savefig('pca.pdf')


