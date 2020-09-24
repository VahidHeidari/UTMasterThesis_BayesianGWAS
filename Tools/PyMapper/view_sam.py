import os.path
import sys

import FASTA_file
import intervals
import logger
import SAM_file



SKIP_SOFT_CLIP     = False
DUMP_EXTRACTED_SAM = False



def NumDigits(num):
	if num == 0:
		return 1

	num = -num if num < 0 else num
	i = 0
	while num > 0:
		i += 1
		num = num // 10
	return i


def GetDigit(dig_idx, num):
	num_digs = NumDigits(num)
	if dig_idx >= num_digs or dig_idx < 0:
		return ' '

	for i in range(dig_idx):
		num = num // 10
	num = num % 10
	return str(num)



def WriteSAMToFile(ref_name, ref_start, ref_end, fasta_file, sam_file, output_file):
	# Find all alignment withint ref [start, end] boundary.
	logger.Log('Reading SAM file . . .')
	sam = SAM_file.TruncatedSAMFile(sam_file, ref_name, ref_start, ref_end)
	max_seq_len = max([ len(SAM_file.GetQName(rec)) for rec in sam.reads ]) + 1
	left_margin = ''.join([' ' for i in range(max(0, ref_start - sam.min_sam_pos) + max_seq_len)])
	if DUMP_EXTRACTED_SAM:
		sam.Dump('extr-sam.sam')

	# Find reference assembly sequence within [start, end] boundary.
	logger.Log('Reading FASTA file . . .')
	fasta = FASTA_file.TruncatedFASTAFile(fasta_file, ref_name, ref_start, ref_end)

	# Adjust ref_end if it is -1.
	if len(fasta.sequence) == 0:
		logger.LogError('Reference chromosome does not exists in FASTA file!')
		exit(2)

	# Adjust ref_end if it is equal to -1 or greater than end of reference.
	ref_seq_len = ref_start + len(fasta.sequence) - 1
	if ref_end == -1:
		ref_end = ref_seq_len
	ref_end = min(ref_end, ref_seq_len)

	# Open output file.
	logger.Log('Writing to output . . .')
	f = open(output_file, 'w')

	# Write position coordinates.
	num_coor_rows = NumDigits(ref_end)
	for r in range(num_coor_rows):
		intr_idx = 0 if len(sam.intr.intervals) else -1
		f.write('    ' if r < num_coor_rows - 1 else 'Coor')
		f.write(left_margin)
		for i in range(ref_start, ref_end + 1):
			dig_idx = num_coor_rows - r - 1
			dig = GetDigit(dig_idx, i)
			f.write(dig)
			if intr_idx >= 0 and sam.intr.intervals[intr_idx][0] - 1 == i:
				intr_len = sam.intr.intervals[intr_idx][1] - sam.intr.intervals[intr_idx][0]
				f.write(''.join([' ' for j in range(intr_len)]))
				intr_idx = intr_idx + 1 if intr_idx + 1 < len(sam.intr.intervals) else -1
		f.write('\n')

	# Write reference sequence.
	f.write('ref ' + left_margin)
	ref_seq = ''
	ref_pos = ref_start
	for itr in sam.intr.intervals:
		ref_seq += fasta.GetSequence(ref_pos, itr[0])
		ref_seq += ''.join(['*' for i in range(itr[1] - itr[0])])
		ref_pos = itr[0]
	ref_seq += fasta.GetSequence(ref_pos, ref_end)
	f.write(ref_seq)
	f.write('\n\n')

	# Write sequences alignment.
	for i in range(len(sam.reads)):
		rec = sam.reads[i]

		# Prepare QNAME.
		qname = SAM_file.GetQName(rec)
		if SAM_file.HasFlag(rec, 'FLAGS_REVERSED_COMPLEMENTED'):
			qname = '-' + qname
		fmt = '{0:>' + str(max_seq_len) + 's}    '
		f.write(fmt.format(qname))

		# Write SEQ padding white space.
		pos = SAM_file.GetSAMPosWithSoftClip(rec) - min(ref_start, sam.min_sam_pos)
		f.write(''.join([' ' for i in range(pos)]))		# Write white spaces to start of sequence position.

		# Write white spaces for padded reference.
		for itr in sam.intr.intervals:
			if itr[0] - min(ref_start, sam.min_sam_pos) >= pos:
				break

			f.write(''.join([' ' for j in range(itr[1] - itr[0])]))

		# Prepare SEQ.
		ops, ops_len = SAM_file.SplitCIGAR(rec)
		seq = SAM_file.GetSeq(rec)
		idx = 0
		for i in range(len(ops)):
			SKIP_OPS = [ 'D', 'P', 'N', 'H' ]
			if ops[i] in SKIP_OPS:
				SKIP_CHAR = ['*', '-', '.', ''][SKIP_OPS.index(ops[i])]
				f.write(''.join([ SKIP_CHAR for j in range(ops_len[i]) ]))
				continue

			sub_seq = seq[idx : idx + ops_len[i]]
			if ops[i] == 'S':
				if SKIP_SOFT_CLIP:
					sub_seq = ''.join([' ' for p in range(ops_len[i])])
				else:
					sub_seq = sub_seq.lower()
			f.write(sub_seq)
			idx += ops_len[i]

		f.write('\n')
	f.close()



def main():
	if len(sys.argv) < 3:
		exe_name = os.path.split(sys.argv[0])[1]
		print('Usage : ' + exe_name + '    ref_start_end   fasta_file   sam_file   [out_path]')
		print('')
		print('            ref_start_end : Fasta reference assembly 1-based start to end position')
		print('                            format is CHR:ST,ED   CHR is chromosome name, ST start position,')
		print('                            and ED is end position, and -1 for Ed means to end of reference')
		print('            fasta_file    : Fasta file path of reference assembly')
		print('            sam_file      : SAM file path for alignments')
		print('            out_path      : Optional output path')
		exit(1)

	# Process command line inputs.
	ref_name = sys.argv[1][0 : sys.argv[1].index(':')]
	ref_start = int(sys.argv[1][sys.argv[1].index(':') + 1 : sys.argv[1].index(',')])
	ref_end = int(sys.argv[1][sys.argv[1].index(',') + 1 : ])
	fasta_file = sys.argv[2]
	sam_file = sys.argv[3]

	# Prepare output file.
	output_file = ref_name + '_' + str(ref_start) + '_' + str(ref_end) + \
			'_' + os.path.split(fasta_file)[1] + '_' + os.path.split(sam_file)[1] + '.txt'
	if len(sys.argv) > 4:
		output_file = os.path.join(sys.argv[4], output_file)
	logger.Log('Output file is -> ' + output_file)

	# Write SAM to output file.
	WriteSAMToFile(ref_name, ref_start, ref_end, fasta_file, sam_file, output_file)


if __name__ == '__main__':
	main()
