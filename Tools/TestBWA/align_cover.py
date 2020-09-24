import matplotlib.pyplot as plt
import numpy as np
import os
import sys



GENO_SIZE   = 500
HEIGHT_STEP = 0.1



READS_BASE   = 'SampleReads'
SAM_PATH     = os.path.join(READS_BASE, 'Aligned.out.sam')
SAM_POS_PATH = os.path.join(READS_BASE, 'pos-aligned.txt')



def SAM2Pos(sam_path, out_path):
	min_pos = 9999999999
	max_pos = -1
	pos = []
	pos_str = ''
	sam_f = open(sam_path, 'r')
	pos_f = open(out_path, 'w')
	for l in sam_f:
		# Skip comments.
		if l.startswith('@'):
			continue

		# Parse start and end position.
		sp = l.split('\t')
		seq_id = sp[0]
		st = sp[3]
		seq_len = len(sp[9].strip())
		ed = str(int(st) + int(seq_len))

		# Format pos string.
		pos_str = seq_id + '\t' + st + '\t' + ed + '\n'
		pos_f.write(pos_str)

		# Append postion to the list.
		rec = (int(st) - 1, int(ed) - 1)
		pos.append(rec)

		# Update genome range.
		min_pos = min(rec[0], min_pos)
		max_pos = max(rec[1], max_pos)

	# Close input and output files.
	sam_f.close()
	pos_f.close()

	# Make genome range and return outputs.
	genome_range = [max(0, min_pos - 10), max_pos + 10]
	return genome_range, pos



def HasIntersection(p1, p2):
	return p1[0] <= p2[1] and p1[1] >= p2[0]



def HasDrawnIntersection(p, drawn_pos):
	for dp in drawn_pos:
		if HasIntersection(p, dp):
			return True
	return False



def GetPoints(pos):
	num_drawn = 0
	NUM_POS = len(pos)
	h = HEIGHT_STEP
	ys = [ 0 for i in range(NUM_POS) ]
	while num_drawn != NUM_POS:
		for i in range(0, NUM_POS):
			drawn_pos = []
			if ys[i] != 0 or HasDrawnIntersection(pos[i], drawn_pos):
				continue

			num_drawn += 1
			ys[i] = h
			drawn_pos.append(pos[i])
			for j in range(i + 1, NUM_POS):
				if ys[j] != 0 or HasDrawnIntersection(pos[j], drawn_pos):
					continue

				num_drawn += 1
				ys[j] = h
				drawn_pos.append(pos[j])
			h += HEIGHT_STEP
	minx = [ p[0] for p in pos ]
	maxx = [ p[1] for p in pos ]
	return ys, minx, maxx



def main():
	# Read input SAM file and write postions in pos_file.
	print('reading SAM file . . .')
	genome_range, pos = SAM2Pos(SAM_PATH, SAM_POS_PATH)
	print('genome range is : ' + str(genome_range))

	# Arrange cover height and x positions.
	print('find Ys . . .')
	ys, minx, maxx = GetPoints(pos)

	# Plot read cover.
	print('plot read cover . . .')
	plt.plot(genome_range, [0, 0])
	plt.hlines(ys, minx, maxx)
	plt.show()



if __name__ == '__main__':
	main()

