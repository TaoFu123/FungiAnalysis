import argparse
import numpy as np
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument('mummer_coords')
parser.add_argument('-q', '--stats-qry', dest='stats_qry', action='store_true')
parser.add_argument('-s', '--stats-sbj', dest='stats_sbj', action='store_true')
parser.add_argument('-l', '--lengths')
args = parser.parse_args()

lengths = {}
with open(args.lengths, 'r') as LN:
	for line in LN:
		if line.startswith('#'):
			continue
		chrom, leng = line.rstrip().split()[:2]
		lengths[chrom] = int(leng)


match_map = defaultdict(list)
is_content = False
with open(args.mummer_coords, 'r') as FI:
	for line in FI:
		if line.startswith('==='):
			is_content = True
			continue
		
		if is_content:
			entry = line.rstrip().split('|')
			if args.stats_sbj:
				start, end = sorted(map(int, entry[0].strip().split()))
				chrom = entry[-1].strip().split()[0]
				match_map[chrom].append([start, end])
			elif args.stats_qry:
				start, end = sorted(map(int, entry[1].strip().split()))
				chrom = entry[-1].strip().split()[1]
				match_map[chrom].append([start, end])

covs = {}
tcov = 0
tcov_len = 0
for c in lengths:
	arr = np.zeros(lengths[c], dtype=bool)
	for s, e in match_map[c]:
		arr[s-1:e] = True
	
	cov_len = np.sum(arr)
	covs[c] = round(cov_len/lengths[c], 4)
	
	tcov += cov_len
	tcov_len += lengths[c]

	print(c + '\t' + str(covs[c]))

print('total_t\t' + str(round(tcov/tcov_len, 4)))
print('total\t' + str(round(tcov/sum(lengths.values()), 4)))

