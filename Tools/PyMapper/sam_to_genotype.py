import os
import sys

import logger
import SAM_file



NUCLEOTIDE_FLAGS = {
	'A' : 0x01,
	'C' : 0x02,
	'G' : 0x04,
	'T' : 0x08,
	'N' : 0x10
}

APPROX_REC_SIZE  = (24 + 24)
MAX_REF_POS      = 999999999999
MAX_MEMORY_LIMIT = 1024 * 1024 * 10



class GenotypeRec:
	def __init__(self, pos, nucleotides):
		self.pos = pos
		self.nucleotides = nucleotides


	def CheckEq(self, other):
		return self.pos == other.pos and	\
			self.nucleotides == other.nucleotides


	def ToString(self, op):
		nucleotides_str = ''
		for k in NUCLEOTIDE_FLAGS:
			flg = NUCLEOTIDE_FLAGS[k]
			if (flg & self.nucleotides) == flg:
				nucleotides_str += k
		return 'op:' + op + '  pos:' + str(self.pos) + '  nuc:' + nucleotides_str


class GenotypeTestRec(GenotypeRec):
	def __init__(self, op, pos, nucleotides):
		GenotypeRec.__init__(self, pos, nucleotides)
		self.op = op


class PrintFunctor:
	def Process(self, r, op):
		print(r.ToString(op))



class WriteFunctor:
	def __init__(self, out_file):
		self.out_file = out_file


	def Process(self, r, op):
		if op == '*':
			pos_str = ''
		else:
			pos_str = str(r.pos)
		self.out_file.write(op + str(pos_str) + '\t')
		for k in NUCLEOTIDE_FLAGS:
			flg = NUCLEOTIDE_FLAGS[k]
			if (flg & r.nucleotides) == flg:
				self.out_file.write(k)
		self.out_file.write('\n')



class GenotypeList:
	def __init__(self):
		self.genotypes = []
		self.inserts = []


	def Reset(self):
		self.genotypes = []
		self.inserts = []


	def Add(self, sam_pos, nucleotide):
		idx = self.GetGenotypeIdxBinSearch(sam_pos)
		if idx == len(self.genotypes):
			new_rec = GenotypeRec(sam_pos, NUCLEOTIDE_FLAGS[nucleotide])
			self.genotypes.append(new_rec)
			return

		if self.genotypes[idx].pos == sam_pos:
			flg = NUCLEOTIDE_FLAGS[nucleotide]
			if (self.genotypes[idx].nucleotides & flg) != flg:
				new_nucleotides = self.genotypes[idx].nucleotides | flg
				new_rec = GenotypeRec(sam_pos, new_nucleotides)
				self.genotypes[idx] = new_rec
		else:
			new_rec = GenotypeRec(sam_pos, NUCLEOTIDE_FLAGS[nucleotide])
			self.genotypes.insert(idx, new_rec)


	def AddInsert(self, sam_pos, insert_idx, nucleotide):
		sam_pos -= 1
		idx = self.GetGenotypeIdx(sam_pos)
		if idx == len(self.inserts):
			new_rec = GenotypeRec(sam_pos, NUCLEOTIDE_FLAGS[nucleotide])
			self.inserts.append(new_rec)
			return

		# Move to insert_idx, extend if needed.
		for i in range(insert_idx):
			# Extend records.
			if len(self.inserts) == idx or self.inserts[idx].pos != sam_pos:
				new_rec = GenotypeRec(sam_pos, 0)
				self.inserts.insert(idx, new_rec)
			idx += 1
		# Check last index, and extend if needed.
		if len(self.inserts) == idx or self.inserts[idx].pos != sam_pos:
			new_rec = GenotypeRec(sam_pos, 0)
			self.inserts.insert(idx, new_rec)

		# Append nucleotide at current index.
		flg = NUCLEOTIDE_FLAGS[nucleotide]
		if (flg & self.inserts[idx].nucleotides) != flg:
			new_nucleotides = self.inserts[idx].nucleotides | flg
			new_rec = GenotypeRec(sam_pos, new_nucleotides)
			self.inserts[idx] = new_rec


	def GetGenotypeIdx(self, sam_pos):
		if len(self.inserts) == 0:
			return 0

		for i in range(len(self.inserts)):
			if self.inserts[i].pos >= sam_pos:
				return i

		return len(self.inserts)


	def GetGenotypeIdxBinSearch(self, sam_pos):
		if len(self.genotypes) == 0:
			return 0

		left = 0
		right = len(self.genotypes) - 1
		while left <= right:
			if self.genotypes[left].pos == sam_pos:
				return left

			if self.genotypes[right].pos == sam_pos:
				return right

			mid = left + (right - left) // 2
			if self.genotypes[mid].pos == sam_pos:
				return mid

			if self.genotypes[mid].pos < sam_pos:
				left = mid + 1
			else:
				right = mid - 1

		return len(self.genotypes)



	def CheckEq(self, genotype_list):
		if len(self.genotypes) + len(self.inserts) != len(genotype_list):
			print('Total lengths are not equal -> ' +	\
				str(len(self.genotypes)) + ' + ' +		\
				str(len(self.inserts)) + ' != ' +		\
				str(len(genotype_list)))
			return False

		num_normal_genotype = 0
		num_insert_genotype = 0
		for rec in genotype_list:
			if rec.op == ' ':
				num_normal_genotype += 1
			else:
				num_insert_genotype += 1
		if len(self.genotypes) != num_normal_genotype:
			print('Number of noraml genotypes are not equal -> ' +	\
				str(len(self.gentoypes)) + ' != ' + str(num_normal_genotype))
			return False

		if len(self.inserts) != num_insert_genotype:
			print('Number of insert gnotypes are not equal -> ' +	\
				str(len(self.inserts)) + ' != ' + str(num_insert_genotype))
			return False

		g_idx = i_idx = 0
		for i in range(len(genotype_list)):
			rec = None
			if genotype_list[i].op == ' ':
				rec = self.genotypes[g_idx]
				g_idx += 1
			else:
				rec = self.inserts[i_idx]
				i_idx += 1
			if not rec.CheckEq(genotype_list[i]):
				rec_idx = g_idx if rec.op == ' ' else i_idx
				print('rec[' + str(rec_idx) +				\
				'] not eq list[' + str(i) + '] -> \n' +		\
					'     ' + rec.ToString() + '\n'			\
					'     ' + genotype_list[i].ToString())
				return False

		return True


	def PrintLists(self):
		for r in self.genotypes:
			print(r.ToString(' '))
		print('-----')
		for r in self.inserts:
			print(r.ToString('*'))


	def PrintListsCombined(self, functor=PrintFunctor(), max_ref_pos=MAX_REF_POS):
		g_idx = i_idx = 0
		num_combined_recs = len(self.genotypes) + len(self.inserts)
		for i in range(num_combined_recs):
			if (g_idx < len(self.genotypes) and i_idx < len(self.inserts) and		\
					self.genotypes[g_idx].pos <= self.inserts[i_idx].pos) or			\
					i_idx >= len(self.inserts):
				if self.genotypes[g_idx].pos >= max_ref_pos:
					return

				functor.Process(self.genotypes[g_idx], ' ')
				g_idx += 1
			elif (i_idx < len(self.inserts) and g_idx < len(self.genotypes) and		\
					self.inserts[i_idx].pos <= self.genotypes[g_idx].pos) or			\
					g_idx >= len(self.genotypes):
				if self.inserts[i_idx].pos >= max_ref_pos:
					return

				functor.Process(self.inserts[i_idx], '*')
				i_idx += 1


	def Dump(self, out_file):
		wf = WriteFunctor(out_file)
		self.PrintListsCombined(wf)


	def Flush(self, max_ref_pos, out_file):
		wf = WriteFunctor(out_file)
		self.PrintListsCombined(wf, max_ref_pos)
		while len(self.genotypes) and self.genotypes[0].pos < max_ref_pos:
			self.genotypes.pop(0)
		while len(self.inserts) and self.inserts[0].pos < max_ref_pos:
			self.inserts.pop(0)


def GetHumanSize(sz):
	sz_str = [ 'Bytes', 'KB', 'MB', 'GB', 'TB' ]
	i = 0
	while sz // 1000 ** (i + 1):
		i += 1
	sz = float(sz) / (1024 ** i)
	sz_out = str(round(sz, 1)) + ' ' + sz_str[i]
	return sz_out



def WriteGenotype(sam_path, out_path):
	logger.Log('\n\n\n---------- START ----------')
	logger.Log('SAM FILE -> ' + sam_path)
	logger.Log('OUT FILE -> ' + out_path)

	cur_rname = ''
	is_rname_dumped = False
	genotype = GenotypeList()
	sam_file = open(sam_path, 'r')
	genotype_file = open(out_path, 'w')
	num_sam_recs = 0
	for l in sam_file:
		# Skip headers.
		if l.startswith('@'):
			continue

		genotype_mem_size = genotype.genotypes.__sizeof__() * APPROX_REC_SIZE
		num_loci = len(genotype.genotypes) + len(genotype.inserts)
		mem_size = num_loci * APPROX_REC_SIZE
		num_sam_recs += 1
		if num_sam_recs % 1000 == 0:
			logger.Log('Number of SAM records: ' + str(num_sam_recs) + \
				'   cached size: ' + GetHumanSize(genotype_mem_size) +	\
				'   real size: ' + GetHumanSize(mem_size))

		# Split and trim record.
		rec = l.split('\t')
		if len(rec) > 0:
			lst_idx = len(rec) - 1
			rec[lst_idx] = rec[lst_idx].strip()

		if SAM_file.HasFlag(rec, 'FLAGS_UNMAPPED'):
			continue

		rname = SAM_file.GetRName(rec)
		if rname != cur_rname:
			if cur_rname != '':
				# Dump to file and reset for the next chromosome.
				logger.Log('Dump reference `' + cur_rname + '\' to output')
				if not is_rname_dumped:
					genotype_file.write('>>> ' + cur_rname + '\n')
				genotype.Dump(genotype_file)
				genotype.Reset()
				is_rname_dumped = False
			cur_rname = rname

		if SAM_file.GetCIGAR(rec) == '*':
			continue

		# Find genotype list position.
		seq = SAM_file.GetSeq(rec)
		seq_idx = 0
		start_sam_pos = sam_pos = SAM_file.GetPos(rec)
		ops, ops_len = SAM_file.SplitCIGAR(rec)
		for i in range(len(ops)):
			if ops[i] in [ 'H', 'P' ]:
				continue

			if ops[i] == 'S':
				seq_idx += ops_len[i]
				continue

			if ops[i] in [ 'D', 'N' ]:
				sam_pos += ops_len[i]
				continue

			sub_seq = seq[seq_idx : seq_idx + ops_len[i]]
			if ops[i] == 'I':
				j_start = ops_len[i - 1] if i - 1 >= 0 and ops[i - 1] == 'P' else 0
				for j in range(j_start, len(sub_seq)):
					s = sub_seq[j]
					genotype.AddInsert(sam_pos, j, s)
				seq_idx += ops_len[i]
				continue

			for s in sub_seq:
				genotype.Add(sam_pos, s)
				sam_pos += 1
			seq_idx += ops_len[i]

		# Flush lists if it exceeds from memory limits.
		if mem_size > MAX_MEMORY_LIMIT:
			logger.Log('Flush genotypes (' + GetHumanSize(mem_size) + ')   ' +	\
				'num loci: ' + str(len(genotype.genotypes) + len(genotype.inserts)))
			if not is_rname_dumped:
				genotype_file.write('>>> ' + cur_rname + '\n')
				is_rname_dumped = True
			genotype.Flush(start_sam_pos, genotype_file)
			genotype_file.flush()
			num_loci = len(genotype.genotypes) + len(genotype.inserts)
			logger.Log('After flush   real size: ' + GetHumanSize(num_loci * APPROX_REC_SIZE) +	\
				'   num loci: ' + str(num_loci))

	# Dump last chromosome.
	logger.Log('Number of SAM records: ' + str(num_sam_recs))
	logger.Log('Dump last refernece `' + cur_rname + '\' to output')
	genotype_file.write('>>> ' + cur_rname + '\n')
	genotype.Dump(genotype_file)

	# Close files.
	sam_file.close()
	genotype_file.close()
	logger.Log('Done!')



def main():
	if len(sys.argv) < 2:
		print('SAM file name is required!')
		exit(1)

	# Process input command line options.
	sam_file = sys.argv[1]
	out_file = os.path.split(sam_file)[1] + '.gtp'
	if len(sys.argv) > 2:
		out_file = os.path.join(sys.argv[2], out_file)

	WriteGenotype(sam_file, out_file)


if __name__ == '__main__':
	main()
