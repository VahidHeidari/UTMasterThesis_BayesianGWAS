import matplotlib.pyplot as plt
import numpy as np
import sys



GENO_SIZE   = 500
HEIGHT_STEP = 0.1

DEF_PATH    = 'SampleReads\\pos.txt'



def PosPred(a, b):
	if a[0] < b[0]:
		return -1
	if a[0] > b[0]:
		return 1
	if a[1] < b[1]:
		return -1
	if a[1] > b[1]:
		return 1
	return 0



def ReadPos(pos_path):
	pos = []
	with open(pos_path, 'r') as f:
		for l in f:
			sp = l.split('\t')
			st = int(sp[1]) - 1
			ed = int(sp[2]) - 1
			rec = (st, ed)
			pos.append(rec)
	return pos



def HasIntersection(p1, p2):
	return p1[0] <= p2[1] and p1[1] >= p2[0]



def GetPointsSorted(pos):
	pos = sorted(pos, cmp=PosPred)
	h = HEIGHT_STEP
	NUM_POS = len(pos)
	ys = [ 0 for i in range(NUM_POS) ]
	num_drawn = 0
	while num_drawn != NUM_POS:
		for i in range(NUM_POS):
			if ys[i] != 0:
				continue

			ys[i] = h
			num_drawn += 1
			p = pos[i]
			for j in range(i + 1, NUM_POS):
				if ys[j] != 0 or HasIntersection(p, pos[j]):
					continue

				ys[j] = h
				num_drawn += 1
				p = pos[j]
			h += HEIGHT_STEP
	minx = [ p[0] for p in pos ]
	maxx = [ p[1] for p in pos ]
	return ys, minx, maxx



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
			for j in range(0, NUM_POS):
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
	sys.argv.append(DEF_PATH)
	pos_path = sys.argv[1]
	pos = ReadPos(pos_path)
	ys, minx, maxx = GetPoints(pos)
	#ys, minx, maxx = GetPointsSorted(pos)
	plt.plot([0, GENO_SIZE], [0, 0])
	plt.hlines(ys, minx, maxx)
	plt.show()


if __name__ == '__main__':
	main()

