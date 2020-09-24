import intervals



SAM_FIELDS = {
	'QNAME' : 0,
	'FLAG'  : 1,
	'RNAME' : 2,
	'POS'   : 3,
	'MAPQ'  : 4,
	'CIGAR' : 5,
	'RNEXT' : 6,
	'PNEXT' : 7,
	'TLEN'  : 8,
	'SEQ'   : 9,
	'QUAL'  : 10
}

SAM_FLAGS = {
	'FLAGS_MULTI_SEGMENTS'              : 0x0001,
	'FLAGS_ALIGNED'                     : 0x0002,
	'FLAGS_UNMAPPED'                    : 0x0004,
	'FLAGS_NEXT_SEGMENT_UNMAPPED'       : 0x0008,
	'FLAGS_REVERSED_COMPLEMENTED'       : 0x0010,
	'FLAGS_NEXT_SEGMENT_REVERSED'       : 0x0020,
	'FLAGS_FIRST_SEGMENT'               : 0x0040,
	'FLAGS_LAST_SEGMENT'                : 0x0080,
	'FLAGS_SECONDARY_ALIGNMENT'         : 0x0100,
	'FLAGS_NOT_PASSING_QUALITY_CONTROL' : 0x0200,
	'FLAGS_PCR_OR_OPTICAL_DUPLICATED'   : 0x0400
}


def GetQName(rec):
	return rec[SAM_FIELDS['QNAME']]


def GetFlag(rec):
	return int(rec[SAM_FIELDS['FLAG']])


def HasFlag(rec, flag_name):
	flags = GetFlag(rec)
	flag = SAM_FLAGS[flag_name]
	return (flags & flag) == flag


def FlagToString(rec):
	flags_str = ''
	for flag_name in SAM_FLAGS:
		if HasFlag(rec, flag_name):
			if len(flags_str) != 0:
				flags_str += ' | '
			flags_str += flag_name
	return flags_str


def GetCIGAR(rec):
	return rec[SAM_FIELDS['CIGAR']]


def SplitCIGAR(rec):
	cigar = GetCIGAR(rec)
	ops = []
	ops_len = []
	i = 0
	while i < len(cigar):
		j = i
		while j < len(cigar) and '0' <= cigar[j] and cigar[j] <= '9':
			j += 1

		op_len = int(cigar[i : j])
		ops_len.append(op_len)
		ops.append(cigar[j])
		i = j + 1
	return ops, ops_len


def GetPos(rec):
	return int(rec[SAM_FIELDS['POS']])


def GetSeq(rec):
	return rec[SAM_FIELDS['SEQ']]


def GetSeqLen(rec):
	if GetCIGAR(rec) == '*':
		return len(GetSeq(rec))

	seq_len = len(GetSeq(rec))
	ops, ops_len = SplitCIGAR(rec)
	for i in range(len(ops)):
		op = ops[i]
		op_len = ops_len[i]
		if op == 'N':
			seq_len += op_len
		elif op == 'I':
			seq_len -= op_len
	return seq_len


def GetSAMPosWithSoftClip(rec):
	if GetCIGAR(rec) == '*':
		return len(GetSeq(rec))

	pos = GetPos(rec)
	ops, ops_len = SplitCIGAR(rec)
	pos -= ops_len[0] if ops[0] == 'S' else 0
	return pos


def GetRName(rec):
	return rec[SAM_FIELDS['RNAME']]



class SAMFile:
	def __init__(self, path):
		self.header = []
		self.reads = []
		with open(path, 'r') as f:
			for l in f:
				sp = l.split('\t')
				if len(sp) > 0:
					lst_idx = len(sp) - 1
					sp[lst_idx] = sp[lst_idx].strip()

				if sp[0].startswith('@'):
					self.header.append(sp)
				else:
					self.reads.append(sp)


	def IsSorted(self, is_print_stats=False):
		cur_pos = 0
		line_num = 0
		out_of_order_pos = 0
		cur_rname = ''
		out_of_order_rname = 0
		for rec in self.reads:
			rname = rec[SAM_FIELDS['RNAME']]
			if rname != cur_rname:
				out_of_order_rname += 1
				cur_pos = 0
			cur_rname = rname

			pos = int(rec[SAM_FIELDS['POS']])
			if pos < cur_pos:
				out_of_order_pos += 1
			cur_pos = pos

		if is_print_stats:
			print('Num out of order pos   : ' + str(out_of_order_pos))
			print('Num out of order rname : ' + str(out_of_order_rname))

		return out_of_order_pos == 0



class TruncatedSAMFile:
	def __init__(self, path, ref_name, ref_start, ref_end):
		self.header = []
		self.reads = []
		self.intr = intervals.Intervals()
		self.min_sam_pos = 999999999999
		with open(path, 'r') as f:
			for l in f:
				sp = l.split('\t')
				if len(sp) > 0:
					lst_idx = len(sp) - 1
					sp[lst_idx] = sp[lst_idx].strip()

				if sp[0].startswith('@'):
					self.header.append(sp)
				else:
					if ref_name != GetRName(sp):
						continue

					if self.HasIntersection(sp, ref_start, ref_end):
						# Add to reads list.
						self.reads.append(sp)

						# Update minimum position.
						rec_pos = GetSAMPosWithSoftClip(sp)
						self.min_sam_pos = min(self.min_sam_pos, rec_pos)

						# Update multisegment intervals.
						if HasFlag(sp, 'FLAGS_MULTI_SEGMENTS'):
							ops, ops_len = SplitCIGAR(sp)
							for j in range(len(ops)):
								if ops[j] == 'H':
									continue

								if ops[j] in ['I', 'P']:
									self.intr.AddInterval(rec_pos, rec_pos + ops_len[j])
								rec_pos += ops_len[j]


	def HasIntersection(self, rec, ref_start, ref_end):
		seq_start = GetSAMPosWithSoftClip(rec)
		seq_len = GetSeqLen(rec)
		seq_end = seq_start + seq_len
		return ref_start <= seq_end and seq_start <= ref_end


	def Dump(self, path):
		with open(path, 'w') as f:
			for h in self.header:
				for i in range(len(h)):
					f.write(h[i])
					if i + 1 < len(h):
						f.write('\t')
				f.write('\n')
			for r in self.reads:
				for i in range(len(r)):
					f.write(r[i])
					if i + 1 < len(r):
						f.write('\t')
				f.write('\n')



if __name__ == '__main__':
	print('This is a module and should be imported! Don\'t run directly!')
