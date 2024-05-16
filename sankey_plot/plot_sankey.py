import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from pysankey2.datasets import load_fruits
from pysankey2.utils import setColorConf
from pysankey2 import Sankey

matplotlib.use('agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

infile = sys.argv[1]
outfig = sys.argv[2]

data = pd.read_csv(infile, sep='\t', index_col=False, header=0)
cls = list(reversed(['AHP', 'HP', 'AP', 'BLP', 'LP', 'BP', 'MP', 'PL', 'UC']))
colors = setColorConf(len(cls), colors='tab10')
cls_map = dict(zip(cls, colors))


layer_labels = {k:cls for k in data.columns}

print(layer_labels)
sky = Sankey(data, colorDict=cls_map, colorMode="global", stripColor='left', layerLabels=layer_labels)
#sky = Sankey(data, colorDict=cls_map, colorMode="global", stripColor='left')
#sky.layerLabels = {k:cls  for k in sky.layerLabels}
print(sky.layerLabels)
print(sky.labels)
print(sky.colnameMaps)
#print(sky.layerPos)

fig, ax = sky.plot()
fig.savefig(outfig)
