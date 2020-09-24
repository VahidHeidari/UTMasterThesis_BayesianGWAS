
#DEBUG_NEW_LINE = '\n'
DEBUG_NEW_LINE = ''



class FASTAFile:
	def __init__(self, path):
		with open(path, 'r') as f:
			self.chromosomes = []
			self.sequences = []
			for l in f:
				if l.startswith(';'):
					continue

				if l.startswith('>'):
					chromosome_name = l[1:].strip().split()[0]
					self.chromosomes.append(chromosome_name)
					self.sequences.append('')
					continue

				self.sequences[len(self.sequences) - 1] += l.strip()


	def IndexOfChromosome(self, chromosome_name):
		return self.chromosomes.index(chromosome_name)



class TruncatedFASTAFile:
	def __init__(self, path, ref_name, ref_start, ref_end, is_print_chromosomes=False):
		if ref_end == -1:
			ref_end = 999999999999
		if ref_end < ref_start:
			tmp = ref_start; ref_start = ref_end; ref_end = tmp

		self.ref_start = ref_start
		self.ref_end = ref_end
		self.sequence = ''
		with open(path, 'r') as f:
			for l in f:
				if l.startswith(';'):
					continue

				if is_print_chromosomes:
					if l.startswith('>'):
						print(l.strip())

				if l.startswith('>' + ref_name):
					cur_pos = 1
					for l in f:
						if l.startswith('>'):
							return

						line_len = len(l.strip())
						if cur_pos + line_len < ref_start:
							cur_pos += line_len
							continue

						st = max(ref_start - cur_pos, 0)
						ed = min(line_len, ref_end - cur_pos + 1)
						self.sequence += l.strip()[st : ed] + DEBUG_NEW_LINE
						cur_pos += line_len
						if ref_end < cur_pos:
							return


	def GetSequence(self, start_pos, end_pos):
		seq_start = start_pos - self.ref_start
		seq_end = seq_start + (end_pos - start_pos) + 1
		return self.sequence[seq_start : seq_end]
