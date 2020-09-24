import os
import random



GEN_SIZE    = 500
NUM_READS   = 200
READ_LEN    = 50

GENOME_BASE = 'SampleGenome'
READS_BASE  = 'SampleReads'
GENOME_PATH = os.path.join(GENOME_BASE, 'genome.fasta')
GTF_PATH    = os.path.join(GENOME_BASE, 'genome.gtf')
READ_PATH   = os.path.join(READS_BASE, 'reads.fasta')
POS_PATH    = os.path.join(READS_BASE, 'pos.txt')

ANNOTATION = [
	[ (10, 490), [ (190, 290) ] ]
]



def CreateReferenceGenome(genome_size, out_path):
	genome = ''.join([ "ACGT"[int(random.uniform(0, 4))] for i in range(genome_size) ])
	fasta = '>chr1\n' + genome
	f = open(out_path, 'w')
	f.write(fasta)
	f.close()
	print(fasta)
	return genome



def HasIntersection(st1, ed1, st2, ed2):
	return st1 <= ed2 and ed1 >= st2



def IsInside(pt_st, pt_ed, reg_st, reg_ed):
	return pt_st >= reg_st and pt_ed <= reg_ed



def IsIntron(start, end, annotation):
	for gene in annotation:
		st = gene[0][0]
		ed = gene[0][1]
		if not IsInside(start, end, st, ed):
			continue

		for exon in gene[1]:
			est = exon[0]
			eed = exon[1]
			if HasIntersection(start, end, est, eed):
				return False
		return True

	return False



def CreateReads(genome, annotation, read_len, out_path):
	intervals = []
	reads_fasta = ''
	pos = ''
	GEN_SIZE = len(genome)

	for i in range(NUM_READS):
		st = int(random.uniform(0, GEN_SIZE - read_len))        # Generate random start positoin.
		ed = st + READ_LEN                                      # Calculate end position.
		while not IsIntron(st, ed, annotation):
			st = int(random.uniform(0, GEN_SIZE - read_len))    # Generate random start positoin.
			ed = st + READ_LEN                                  # Calculate end position.
		read_seq = genome[st : ed]                              # Create read.
		seq_id = '>Sequence_' + str(i + 1)                      # Make sequence ID.
		pos += seq_id + '\t' + str(st + 1) + '\t' + str(ed + 1) + '\n'  # Make position record.
		reads_fasta += seq_id + '\n' + read_seq + '\n'          # Make sequence itself.
		read_int = (st, ed)                                     # Save read position.
		intervals.append(read_int)
	f = open(out_path, 'w')
	f.write(reads_fasta)
	f.close()
	print(intervals)
	print(reads_fasta)
	return pos, intervals



def WritePos(pos, out_path):
	f = open(out_path, 'w')
	f.write(pos)
	f.close()
	print(pos)


def CreateGTF(annotation, out_path):
	gtf = ''
	gene_id = 1
	for gene in annotation:
		st = gene[0][0]
		exons = gene[1]
		gene_id_str = 'gene_id ' + str(gene_id) + '; transcript_id 1;'
		# Write introns
		for exon in exons:
			ed = exon[0] - 1
			#      seqname source feature  start            end              score strand frame   attrib
			gtf += 'chr1\tsmpl\tCDS\t' +   str(st) + '\t' + str(ed) + '\t' + '.\t+\t0' +          gene_id_str + '\n'
			st = exon[1] + 1
		ed = gene[0][1]
		#      seqname source feature  start            end              score strand frame   attrib
		gtf += 'chr1\tsmpl\tCDS\t' +   str(st) + '\t' + str(ed) + '\t' + '.\t+\t0' +          gene_id_str + '\n'

		# Write exons.
		for exon in exons:
			est = exon[0]
			eed = exon[1]
			#      seqname source feature  start             end               score strand frame  attrib
			gtf += 'chr1\tsmpl\texon\t' +  str(est) + '\t' + str(eed) + '\t' + '.\t+\t0' +         gene_id_str + '\n'

		# Write start and end codons.
		st = gene[0][0]
		ed = gene[0][1]
		#      seqname source feature  start            end              score strand frame             attrib
		gtf += 'chr1\tsmpl\tstart_codone\t' + str(st) + '\t' + str(st + 2) + '\t' + '.\t+\t0' +         gene_id_str + '\n'
		#      seqname source feature         start                end              score strand frame  attrib
		gtf += 'chr1\tsmpl\tend_codone\t' +   str(ed - 2) + '\t' + str(ed) + '\t' + '.\t+\t0' +         gene_id_str + '\n'

		gene_id += 1
	f = open(out_path, 'w')
	f.write(gtf)
	f.close()



def CHECK_TRUE(val):
	if val == True:
		print('OK')
	else:
		print('Failed!')

def CHECK_FALSE(val):
	if val == False:
		print('OK')
	else:
		print('Failed!')

def TestIsIntron():
	TEST_ANNOTATION = [
		[ (10, 490), [ (190, 290) ] ]
	]
	CHECK_TRUE(IsIntron(10, 20, TEST_ANNOTATION))
	CHECK_TRUE(IsIntron(11, 22, TEST_ANNOTATION))
	CHECK_TRUE(IsIntron(40, 60, TEST_ANNOTATION))
	CHECK_TRUE(IsIntron(45, 55, TEST_ANNOTATION))
	CHECK_TRUE(IsIntron(55, 65, TEST_ANNOTATION))
	CHECK_TRUE(IsIntron(61, 89, TEST_ANNOTATION))
	CHECK_FALSE(IsIntron(  1,   9, TEST_ANNOTATION))
	CHECK_FALSE(IsIntron(  5,  15, TEST_ANNOTATION))
	CHECK_FALSE(IsIntron(180, 200, TEST_ANNOTATION))
	CHECK_FALSE(IsIntron(180, 200, TEST_ANNOTATION))
	CHECK_FALSE(IsIntron(450, 500, TEST_ANNOTATION))
	exit(0)



def AreFilesCreated():
	return os.path.isfile(GENOME_PATH) and		\
			os.path.isfile(GTF_PATH) and		\
			os.path.isfile(READ_PATH) and		\
			os.path.isfile(POS_PATH)



def main():
	if not os.path.isdir(GENOME_BASE):
		os.makedirs(GENOME_BASE)
	if not os.path.isdir(READS_BASE):
		os.makedirs(READS_BASE)

	if AreFilesCreated():
		print('Sample files already created!')
		exit(0)

	genome = CreateReferenceGenome(GEN_SIZE, GENOME_PATH)
	pos, intervals = CreateReads(genome, ANNOTATION, READ_LEN, READ_PATH)
	WritePos(pos, POS_PATH)
	CreateGTF(ANNOTATION, GTF_PATH)



if __name__ == '__main__':
	#TestIsIntron()
	main()

