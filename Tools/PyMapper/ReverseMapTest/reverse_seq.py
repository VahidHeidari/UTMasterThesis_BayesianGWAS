import sys



def ReverseComplementSeq(seq):
	NUCLEOTIDES = [ 'A', 'C', 'G', 'T' ]
	COMP_NUCLEOTIDES = [ 'T', 'G', 'C', 'A' ]
	rev_seq = ''
	seq_len = len(seq)
	for i in range(seq_len):
		n_idx = seq_len - i - 1
		n = seq[n_idx]
		nuc_idx = NUCLEOTIDES.index(n)
		c = COMP_NUCLEOTIDES[nuc_idx]
		rev_seq += c
	return rev_seq


def main():
	seq = sys.argv[1]
	rev_seq = ReverseComplementSeq(seq)
	print(rev_seq)


if __name__ == '__main__':
	main()
