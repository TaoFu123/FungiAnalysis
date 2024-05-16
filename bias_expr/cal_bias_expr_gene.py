import sys
import numpy as np
from collections import defaultdict

expr_file = sys.argv[1]
homo_file = sys.argv[2]

fold = 2

expr_data = defaultdict(dict)
with open(expr_file, 'r') as EP:
	next(EP)
	for line in EP:
		entry = line.rstrip().split('\t')
		tis = entry[0].split('-')[1]
		g = entry[3]
		expr_data[g].setdefault(tis, []).append(float(entry[-1]))
		#expr_data[tis].setdefault(g, []).append(float(entry[-1]))

with open(homo_file, 'r') as HF:
	print('pair' + '\t' + '\t'.join(sorted(expr_data[g])))
	for line in HF:
		if line.startswith('#'):
			continue

		g1, g2 = line.rstrip().split('\t')[:2]
		out = [g1 + '@' + g2]
		for tis in sorted(expr_data[g1].keys()):
			tpm1 = np.median(expr_data[g1][tis])
			tpm2 = np.median(expr_data[g2][tis])

			if tpm1 < 1 and tpm2 < 1:
				out.append('b')
			elif tpm1 < 1:
				if tpm2 >= fold:
					out.append('2')
				else:
					out.append('b')
			elif tpm2 < 1:
				if tpm1 >= fold:
					out.append('1')
				else:
					out.append('b')
			else:
				if tpm1 >= tpm2:
					if tpm1/tpm2 >= fold:
						out.append('1')
					else:
						out.append('b')
				else:
					if tpm2/tpm1 >= fold:
						out.append('2')
					else:
						out.append('b')

		print('\t'.join(out))

				
