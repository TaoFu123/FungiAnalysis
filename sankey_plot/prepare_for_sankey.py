import sys
import pandas as pd

'''
	>> cut -f 1,2,3 25.dom.ase.tsv|sort -u
	pnl	cond	H
	R19 NaCl.0.3M	NaCl.0.3M	R19
	R19 NaCl.0.8M	NaCl.0.8M	R19
	R19 NaCl.0M	NaCl.0M	R19
	R19 poplar	poplar	R19
	R8 NaCl.0.3M	NaCl.0.3M	R8
	R8 NaCl.0.8M	NaCl.0.8M	R8
	R8 NaCl.0M	NaCl.0M	R8
	R8 poplar	poplar	R8

'''

infile = sys.argv[1]
spec = sys.argv[2] # R19 or R8
outfile = sys.argv[3]

sts_to_idx = {'NaCl.0.3M': 0,
			'NaCl.0.8M': 1,
			'NaCl.0M' : 2,
			'poplar' : 3}

result = {}
with open(infile, 'r') as FI:
	header = next(FI).rstrip().split('\t')
	for line in FI:
		entry = dict(zip(header, line.rstrip().split('\t')))
		if entry['H'] == spec:
			assert entry['Dom'].strip() or entry['Dom2'].strip(), line
			dm = entry['Dom'].strip() if entry['Dom'].strip() else entry['Dom2'].strip()
			#print(entry['gid'] + '\t' + dm)
			
			if dm == 'PD_L':
				dm = 'LP'
			elif dm == 'PD_H':
				dm = 'HP'
			
			result.setdefault(entry['gid'], ['UC']*4)[sts_to_idx[entry['cond']]] = dm
			# UC: Unclassified

			
data = pd.DataFrame(list(result.values()), columns=['NaCl.0.3M', 'NaCl.0.8M', 'NaCl.0M', 'poplar'])
data.to_csv(outfile, sep='\t', index=False)
