import sys



sys.argv.append('SRR1982742-sorted.sam-1.gtp')
cur_rname = ''
pos = 0
line_num = 0
with open(sys.argv[1], 'r') as f:
	for l in f:
		line_num += 1
		if l.startswith('>>>'):
			cur_rname = l.split()[1].strip()
			pos = 0
			print('cur_rname:' + cur_rname)
			continue

		cur_pos = int(l.split('\t')[0][1:])
		if pos > cur_pos:
			print('Line:', line_num, 'Ref:', cur_rname, 'Pos:', pos, 'CurPos:', cur_pos)
			print(l.strip())
			break

		pos = cur_pos

print('Done!')

