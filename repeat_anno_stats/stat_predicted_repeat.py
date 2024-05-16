import re
import sys
import argparse
import numpy as np
import logging
from collections import defaultdict, OrderedDict

def stat_repeat_category(ingffs, class_database, seq_len, level, chrom):
	if args.chrom:
		chrom_set = set(chrom.split(','))
	else:
		chrom_set = set(seq_len.keys())

	logging.info('Read class database: ' + class_database)
	class_to_categ = {}
	with open(class_database) as RC:
		for line in RC:
			if line.startswith('#'):
				continue

			k, v = line.rstrip().split()
			#logging.debug(k, v)
			class_to_categ[k] = v


	logging.info('Read gff file')
	repeat_region = {}
	unknown_class = set()
	class_regrex = re.compile(r'Class=([^;]+)')
	logging.debug('Inputs:' + str(ingffs))
	for gff in ingffs:
		with open(gff, 'r') as RG:
			for line in RG:
				if line.startswith('#'):
					continue

				entry = line.rstrip().split('\t')
				cls = class_regrex.search(entry[-1]).group(1)
				
				categ = class_to_categ.get(cls)
				if categ is None:
					class_to_categ[cls] = 'Other'
					unknown_class.add(cls)

				start, end = sorted([int(entry[3]), int(entry[4])])
				seqid = entry[0]
				repeat_region.setdefault(seqid, defaultdict(list))[categ].append([start, end])
	
	logging.warning('Repeat class: {} not in category file, will be added to "Other"'.format(' '.join(unknown_class)))


	ouput_order = ['Retro/LTR/Copia', 'Retro/LTR/Gypsy', 'Retro/LTR/Other', 'Retro/LTR/Unknown', 'Retro/SINE', 'Retro/LINE', 'Retro/Other', \
				    'DNA/Academ', 'DNA/CACTA', 'DNA/Crypton', 'DNA/Dada', 'DNA/Ginger', 'DNA/hAT', 'DNA/Helitron', \
					'DNA/Kolobok', 'DNA/Mutator-like', 'DNA/P_Element', 'DNA/PIF-Harbinger', 'DNA/PiggyBac',\
					'DNA/Sola', 'DNA/TcMar', 'DNA/Other', \
					'Other', 'Unknown']
	
	categ_stat = {i:[] for i in ouput_order + ['Total']}

	for seqid in sorted(repeat_region):
		tot_arr = np.zeros(seq_len[seqid], dtype=bool)
		for cat in ouput_order:
			cat_arr = np.zeros(seq_len[seqid], dtype=bool)

			for start, end in repeat_region[seqid][cat]:
				cat_arr[start-1 : end] = True
				
			categ_stat[cat].append(np.sum(cat_arr))
			tot_arr += cat_arr

		categ_stat['Total'].append(np.sum(tot_arr))

	

	if level == 'whole':
		genome_size = sum(seq_len.values())
		print('Type\tLength(bp)\t% of genome')
		for cat in ouput_order + ['Total']:
			output = [cat, str(sum(categ_stat[cat])), str(round(sum(categ_stat[cat])/genome_size *100, 5))]
			print('\t'.join(output))
	
	elif level == 'chr':
		header = ['Type']
		for chrm in sorted(repeat_region):
			if chrm in chrom_set:
				for j in ['Length(bp)', '% of chromosome']:
					header.append(chrm + ':' + j)
		print('\t'.join(header))

		for cat in ouput_order + ['Total']:
			output = [cat]
			for l, seqid in zip(categ_stat[cat], sorted(repeat_region)):
				if seqid in chrom_set:
					output.append(str(l))
					output.append(str(round(l/seq_len[seqid]*100, 5)))
			print('\t'.join(output))


def stat_repeat_by_class(ingffs, class_database, seq_len, outprefix):
	class_to_categ = {}
	categ_to_class = defaultdict(list)
	with open(class_database) as RC:
		for line in RC:
			k, v = line.rstrip().split()
			class_to_categ[k] = v
			categ_to_class[v].append(k)


	unknown_class = set()
	repeat_region = {}
	class_regrex = re.compile(r'Class=([^;]+)')
	logging.debug('Inputs:' + str(ingffs))
	for gff in ingffs:
		with open(gff, 'r') as RG:
			for line in RG:
				if line.startswith('#'):
					continue
		
				entry = line.rstrip().split('\t')
				k = class_regrex.search(entry[-1]).group(1)
			
				v = class_to_categ.get(k)
				if v is None:
					categ_to_class['Other'].append(k)
					class_to_categ[k] = 'Other'
					unknown_class.add(k)


				start, end = sorted([int(entry[3]), int(entry[4])])

				seqid = entry[0]

				repeat_region.setdefault(seqid, defaultdict(list))[k].append([start, end])

	logging.warning('Repeat class: {} not in category file, will be added to "Other"'.format(' '.join(unknown_class)))


	ouput_order = ['Retro/LTR/Copia', 'Retro/LTR/Gypsy', 'Retro/LTR/Other', 'Retro/LTR/Unknown', 'Retro/SINE', 'Retro/LINE', 'Retro/Other', \
				    'DNA/Academ', 'DNA/CACTA', 'DNA/Crypton', 'DNA/Dada', 'DNA/Ginger', 'DNA/hAT', 'DNA/Helitron', \
	                'DNA/Kolobok', 'DNA/Mutator-like', 'DNA/P_Element', 'DNA/PIF-Harbinger', 'DNA/PiggyBac',\
	                'DNA/Sola', 'DNA/TcMar', 'DNA/Other', \
	                'Other', 'Unknown']

	if level == 'whole':
		main_cat_stat = {i:0 for i in ouput_order + ['Total']}
		exha_cls_stat = defaultdict(int)
		genome_size = sum(seq_len.values())

		for seqid in repeat_region:
			total_arr = np.zeros(seq_len[seqid], dtype=bool)
			for cat in ouput_order:
				cat_arr = np.zeros(seq_len[seqid], dtype=bool)

				for cls in categ_to_class[cat]:

					if cls in repeat_region[seqid]:
						cls_arr = np.zeros(seq_len[seqid], dtype=bool)

						for start, end in repeat_region[seqid][cls]:
							cls_arr[start-1 : end] = True

						exha_cls_stat[cls] += np.sum(cls_arr)
						cat_arr += cls_arr

				main_cat_stat[cat] += np.sum(cat_arr)
				total_arr += cat_arr

			main_cat_stat['Total'] += np.sum(total_arr)


		with open(outprefix + '.main_class.stat.xls', 'w') as MC:
			MC.write('Type\tLength(bp)\t% of genome\n')
			for cat in ouput_order + ['Total']:
				output = [cat, str(main_cat_stat[cat]), str(round(main_cat_stat[cat]/genome_size *100, 5))]
				MC.write('\t'.join(output) + '\n')

		with open(outprefix + '.exhaustive_class.stat.xls', 'w') as EC:
			EC.write('Class\tCategory\tLength(bp)\t% of genome\n')
			for cat in ouput_order:
				for cls in categ_to_class[cat]:
					if cls in exha_cls_stat:
						output = [cls, cat, str(exha_cls_stat[cls]), str(round(exha_cls_stat[cls]/genome_size *100, 5))]
						EC.write('\t'.join(output) + '\n')
	
	elif level == 'chr':
		with open(outprefix + '.main_class.stat.xls', 'w') as MC, open(outprefix + '.exhaustive_class.stat.xls', 'w') as EC:
			output = ['Type'] + [chrm + ':' + j for chrm in sorted(repeat_region) for j in ['Length(bp)', '% of chromosome']]
			MC.write('\t'.join(output) + '\n')

			output = ['Class', 'Category'] + [chrm + ':' + j for chrm in sorted(repeat_region) for j in ['Length(bp)', '% of chromosome']]
			EC.write('\t'.join(output) + '\n')

		for cat in ouput_order:
			for seqid in sorted(repeat_region):
				arr = np.zeros(seq_len[seqid], dtype=bool)

				for cls in categ_to_class[cat]:
					rs = repeat_region[seqid].get(cls, None)
					if rs:
						for start, end in rs:
							arr[start-1 : end] = True
							
			
def parse_genome_size(infile):
	logging.info('Read genome')
	import gzip

	if infile.endswith(('gz', 'GZ')):
		FA = gzip.open(infile, 'rt')
	else:
		FA = open(infile)

	size = 0
	for line in FA:
		if line.startswith('>'):
			continue

		sequence = line.rstrip().upper()
		size += len(sequence) - sequence.count('N')

	return size


def parse_seq_len_from_genome(infile):
	logging.info('Read sequence length.')

	import gzip

	if infile.endswith(('gz', 'GZ')):
		FA = gzip.open(infile, 'rt')
	else:
		FA = open(infile)

	seq_len = {}
	for line in FA:
		if line.startswith('>'):
			seqid = line.rstrip().split()[0][1:]
			seq_len[seqid] = 0
		else:
			sequence = line.rstrip().upper()
			seq_len[seqid] += len(sequence) - sequence.count('N')

	return seq_len


def stat_repeat_by_method(types, ingffs, genome_size):
	types = list(types)

	repeats = {}
	seq_lens = defaultdict(int)
	for t, gff in zip(types, ingffs):
		with open(gff, 'r') as GF:
			for line in GF:
				if line.startswith(('#', '\n')):
					continue

				entry = line.rstrip().split('\t')
				seqid = entry[0]
				start, end = sorted(map(int, entry[3:5]))
				
				repeats.setdefault(seqid, defaultdict(list))[t].append([start, end])
				
				if end > seq_lens[seqid]:
					seq_lens[seqid] = end


	stats = OrderedDict((t, 0) for t in types + ['Total'])

	for seqid in repeats:
		tarr = np.zeros(seq_lens[seqid], dtype=bool)

		for t in types:
			arr = np.zeros(seq_lens[seqid], dtype=bool)

			for s, e in repeats[seqid][t]:
				arr[s-1: e] = True
				tarr[s-1: e] = True

			stats[t] += np.sum(arr)
		
		stats['Total'] += np.sum(tarr)
		
	
	print('Type\tLength\t% of genome')
	for t in stats:
		print('{}\t{}\t{}'.format(t, stats[t], round(stats[t]/genome_size*100, 4)))



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	subparser = parser.add_subparsers()

	mth_parser = subparser.add_parser('method')
	mth_parser.add_argument('-t', '--type', action='append', help=' TRF RepeatProteinMask RepeatMasker denovo LTRharvest  ')
	mth_parser.add_argument('-i', '--gff', action='append')
	mth_parser.add_argument('-g', '--genome', help='input genome file for genome size')
	#mth_parser.add_argument('--genomeSize', dest='genome_size', type=int)
	mth_parser.set_defaults(action='method')


	cls_parser = subparser.add_parser('class')
	cls_parser.add_argument('gff', nargs='+')
	cls_parser.add_argument('-g', '--genome', help='input genome file for genome size')
	cls_parser.add_argument('-c', '--chrom', help='When --level is set chr, you can specify chromosome to stats, example: chr01,chr02')
	#cls_parser.add_argument('--genomeSize', dest='genome_size', type=int)
	cls_parser.add_argument('-l', '--level', default='whole', choices={'whole', 'chr'}, help='Statistic level, whole genome or single chromosome, when it is set chr, option "--genome" is required.')
	cls_parser.add_argument('-e', '--exhaustive', action='store_true', help='Output exhaustive statistic of repeat class')
	cls_parser.add_argument('--classDatabase', default='/home/futao/bin/pipe/ann/repeat/scripts/Transposon.class.map')
	cls_parser.add_argument('--outprefix', default='transopon')
	cls_parser.set_defaults(action='class')

	args = parser.parse_args(sys.argv[1:] if sys.argv[1:] else ['-h'])

	logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.DEBUG, datefmt='%Y/%m/%d %H:%M:%S')


	if args.action == 'method':
		if not args.type:
			args.type = 'TRF RepeatProteinMask RepeatMasker denovo LTRharvest'.split()

		genome_size = parse_genome_size(args.genome) if args.genome else args.genome_size
		stat_repeat_by_method(args.type, args.gff, genome_size)

	elif args.action == 'class':
		#genome_size = parse_genome_size(args.genome) if args.genome else args.genome_size
		seq_len = parse_seq_len_from_genome(args.genome)
		if args.exhaustive:
			stat_repeat_by_class(args.gff, args.classDatabase, genome_size, args.exhaustive)
		else:
			stat_repeat_category(args.gff, args.classDatabase, seq_len, args.level, args.chrom)	

