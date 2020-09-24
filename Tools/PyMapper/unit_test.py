import FASTA_file
import intervals
import SAM_file
import sam_to_genotype



IS_DUMP_GENOTYPE_LISTS = False



is_failed = False
def CheckEq(val, exp):
	if val == exp:
		print('OK!')
	else:
		global is_failed
		is_failed = True
		print('Failed!')
		print('val -> ' + str(val))
		print('exp -> ' + str(exp))


def CheckTrue(val):
	CheckEq(val, True)


def CheckFalse(val):
	CheckEq(val, False)



def TestSAMFile(path, num_hdrs, num_reads, is_sorted):
	print('Reading SAM file from `' + path + '\' . . .')
	sam = SAM_file.SAMFile(path)
	print('HDR       : ' + str(len(sam.header)))
	CheckEq(len(sam.header), num_hdrs)
	print('RDS       : ' + str(len(sam.reads)))
	CheckEq(len(sam.reads), num_reads)
	print('HDR + RDS : ' + str(len(sam.header) + len(sam.reads)))
	CheckEq(len(sam.header) + len(sam.reads), num_hdrs + num_reads)

	print('Check sorted by coordinate . . .')
	is_sam_sorted = sam.IsSorted(True)
	print('Is Sorted : ' + str(is_sam_sorted))
	CheckEq(is_sam_sorted, is_sorted)
	return sam


def TestSortedSAMFile():
	print('--- TestSortedSAMFile() ---')
	sam = TestSAMFile('D:\\Datasets\\Programs\\Mapper\\SAMSpecSamples\\SAM_BAM_spec.sam', 2, 6, True)
	ops, ops_len = SAM_file.SplitCIGAR(sam.reads[0])
	CheckEq(ops, ['M', 'I', 'M', 'D', 'M'])
	CheckEq(ops_len, [8, 2, 4, 1, 3])

	sam = TestSAMFile('D:\\Datasets\\GSE68086\\Dataset\\FastQ\\HC\\SRR2095004-sorted.sam', 197, 3320807, True)
	print('Checking flags . . .')
	for rec in sam.reads:
		if SAM_file.HasFlag(rec, 'FLAGS_REVERSED_COMPLEMENTED'):
			print(rec)
			print('FLAG  : ' + SAM_file.FlagToString(rec))
			print('CIGAR : ' + SAM_file.GetCIGAR(rec))
			print('POS   : ' + str(SAM_file.GetPos(rec)))
			print('SEQ   : ' + SAM_file.GetSeq(rec))
			print(len(SAM_file.GetSeq(rec)))
			break



def TestIntersectionFunctions():
	CheckTrue(intervals.HasIntersection((10, 100), (0, 10)))
	CheckFalse(intervals.IsInside((0, 10), (10, 100)))

	CheckTrue(intervals.HasIntersection((10, 100), (20, 30)))
	CheckTrue(intervals.IsInside((20, 30), (10, 100)))
	CheckFalse(intervals.IsInside((10, 100), (20, 30)))

	CheckFalse(intervals.HasIntersection((10, 20), (40, 50)))
	CheckFalse(intervals.IsInside((10, 20), (40, 50)))


def TestAddNonIntersecting():
	intr = intervals.Intervals()
	intr.AddInterval(20, 30)
	CheckEq(intr.intervals, [(20, 30)])
	intr.AddInterval(220, 330)
	CheckEq(intr.intervals, [(20, 30), (220, 330)])
	intr.AddInterval(50, 100)
	CheckEq(intr.intervals, [(20, 30), (50, 100), (220, 330)])
	intr.AddInterval(150, 200)
	CheckEq(intr.intervals, [(20, 30), (50, 100), (150, 200), (220, 330)])
	intr.AddInterval(5, 10)
	CheckEq(intr.intervals, [(5, 10), (20, 30), (50, 100), (150, 200), (220, 330)])
	intr.AddInterval(500, 600)
	CheckEq(intr.intervals, [(5, 10), (20, 30), (50, 100), (150, 200), (220, 330), (500, 600)])


def TestAddIntersecting():
	intr = intervals.Intervals()
	intr.AddInterval(20, 30)
	CheckEq(intr.intervals, [(20, 30)])
	intr.AddInterval(15, 25)
	CheckEq(intr.intervals, [(15, 30)])
	intr.AddInterval(25, 35)
	CheckEq(intr.intervals, [(15, 35)])


def TestAddInside():
	intr = intervals.Intervals()
	intr.AddInterval(20, 30)
	CheckEq(intr.intervals, [(20, 30)])
	intr.AddInterval(22, 28)
	CheckEq(intr.intervals, [(20, 30)])


def TestAddMerge():
	intr = intervals.Intervals()
	intr.AddInterval(20, 30)
	CheckEq(intr.intervals, [(20, 30)])
	intr.AddInterval(40, 50)
	CheckEq(intr.intervals, [(20, 30), (40, 50)])
	intr.AddInterval(60, 70)
	CheckEq(intr.intervals, [(20, 30), (40, 50), (60, 70)])
	intr.AddInterval(25, 45)
	CheckEq(intr.intervals, [(20, 50), (60, 70)])


def RunIntervalsTests():
	print('--- RunIntervalsTests() ---')
	TestIntersectionFunctions()
	TestAddNonIntersecting()
	TestAddIntersecting()
	TestAddInside()
	TestAddMerge()



def TestSingleChromosomeFASTAFile():
	HG38_REF_FILE = 'TruncatedFASTAFileTest\\hg38.fasta'

	fasta = FASTA_file.TruncatedFASTAFile(HG38_REF_FILE, '1', 166 * 60 + 40, 167 * 60)
	CheckEq(fasta.sequence, 'NTAACCCTAACCCTAACCCTA')

	fasta = FASTA_file.TruncatedFASTAFile(HG38_REF_FILE, '1', 167 * 60 + 1, 168 * 60)
	CheckEq(fasta.sequence, 'aCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTb')

	fasta = FASTA_file.TruncatedFASTAFile(HG38_REF_FILE, '1', 166 * 60 + 38, 166 * 60 + 43)
	CheckEq(fasta.sequence, 'NNNTAA')

	fasta = FASTA_file.TruncatedFASTAFile(HG38_REF_FILE, '1', 9961, 10019)
	CheckEq(fasta.sequence, 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCT')

	fasta = FASTA_file.TruncatedFASTAFile(HG38_REF_FILE, '1', 166 * 60 + 1, 167 * 60)
	CheckEq(fasta.sequence, 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTA')

	fasta = FASTA_file.TruncatedFASTAFile(HG38_REF_FILE, '1', 166 * 60 + 40, 168 * 60)
	CheckEq(fasta.sequence,                'NTAACCCTAACCCTAACCCTA' +	\
	'aCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTb')

	fasta = FASTA_file.TruncatedFASTAFile(HG38_REF_FILE, '1', 166 * 60 + 40, 168 * 60 + 1)
	CheckEq(fasta.sequence,                'NTAACCCTAACCCTAACCCTA' +	\
	'aCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTb' +	\
	'c')

	fasta = FASTA_file.TruncatedFASTAFile(HG38_REF_FILE, '1', 166 * 60 + 40, 168 * 60 + 41)
	CheckEq(fasta.sequence,                'NTAACCCTAACCCTAACCCTA' +	\
	'aCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTb' +	\
	'cCCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTA')

	fasta = FASTA_file.TruncatedFASTAFile(HG38_REF_FILE, '1', 166 * 60 + 40, 169 * 60 + 41)
	CheckEq(fasta.sequence,                'NTAACCCTAACCCTAACCCTA' +	\
	'aCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTb' +	\
	'cCCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAd' +	\
	'eCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTA')

	fasta = FASTA_file.TruncatedFASTAFile(HG38_REF_FILE, '1', 166 * 60 + 40, 170 * 60)
	CheckEq(fasta.sequence,                'NTAACCCTAACCCTAACCCTA' +	\
	'aCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTb' +	\
	'cCCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAd' +	\
	'eCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAf')

	REF_FILE = 'TruncatedFASTAFileTest\\reference.fasta'

	fasta = FASTA_file.TruncatedFASTAFile(REF_FILE, 'chr2', 1, 10)
	CheckEq(fasta.sequence, 'QWERTYUIOP')

	fasta = FASTA_file.TruncatedFASTAFile(REF_FILE, 'chr2', 1, 20)
	CheckEq(fasta.sequence, 'QWERTYUIOP' +	\
	'ASDFGHJKL:')

	fasta = FASTA_file.TruncatedFASTAFile(REF_FILE, 'chr2', 1, 30)
	CheckEq(fasta.sequence, 'QWERTYUIOP' +	\
	'ASDFGHJKL:' +	\
	'ZXCVBNM<>?')

	fasta = FASTA_file.TruncatedFASTAFile(REF_FILE, 'chr2', 1, 40)
	CheckEq(fasta.sequence, 'QWERTYUIOP' +	\
	'ASDFGHJKL:' +	\
	'ZXCVBNM<>?' +	\
	'[];\',./=-\\')

	fasta = FASTA_file.TruncatedFASTAFile(REF_FILE, 'chr2', 7, 10)
	CheckEq(fasta.sequence, 'UIOP')

	fasta = FASTA_file.TruncatedFASTAFile(REF_FILE, 'chr2', 7, 20)
	CheckEq(fasta.sequence, 'UIOP' +	\
	'ASDFGHJKL:')

	fasta = FASTA_file.TruncatedFASTAFile(REF_FILE, 'chr2', 7, 30)
	CheckEq(fasta.sequence, 'UIOP' +	\
	'ASDFGHJKL:' +	\
	'ZXCVBNM<>?')

	fasta = FASTA_file.TruncatedFASTAFile(REF_FILE, 'chr2', 7, 40)
	CheckEq(fasta.sequence, 'UIOP' +	\
	'ASDFGHJKL:' +	\
	'ZXCVBNM<>?' +	\
	'[];\',./=-\\')

	fasta = FASTA_file.TruncatedFASTAFile(REF_FILE, 'chr2', 7, 34)
	CheckEq(fasta.sequence, 'UIOP' +	\
	'ASDFGHJKL:' +	\
	'ZXCVBNM<>?' +	\
	'[];\'')



def PrintAndDump(genotype, dump_path):
	print('--- Print Lists ---')
	genotype.PrintLists()
	print('--- Print Lists Combined ----')
	genotype.PrintListsCombined()
	if IS_DUMP_GENOTYPE_LISTS:
		with open(dump_path, 'w') as f:
			genotype.Dump(f)

def TestGenotypeListAddInsert():
	genotype = sam_to_genotype.GenotypeList()

	# Test normal sequence adding.
	genotype.Add(7, 'A')
	genotype.Add(8, 'C')
	genotype.Add(9, 'C')
	genotype.Add(10, 'G')
	genotype.Add(11, 'T')
	test_list = [
		sam_to_genotype.GenotypeTestRec(' ', 7, sam_to_genotype.NUCLEOTIDE_FLAGS['A']),
		sam_to_genotype.GenotypeTestRec(' ', 8, sam_to_genotype.NUCLEOTIDE_FLAGS['C']),
		sam_to_genotype.GenotypeTestRec(' ', 9, sam_to_genotype.NUCLEOTIDE_FLAGS['C']),
		sam_to_genotype.GenotypeTestRec(' ', 10, sam_to_genotype.NUCLEOTIDE_FLAGS['G']),
		sam_to_genotype.GenotypeTestRec(' ', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['T'])
	]
	CheckTrue(genotype.CheckEq(test_list))

	# Test insert sequence adding at the end.
	genotype.AddInsert(12, 0, 'C')
	genotype.AddInsert(12, 1, 'G')
	genotype.AddInsert(12, 2, 'A')
	test_list.append(sam_to_genotype.GenotypeTestRec('*', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['C']))
	test_list.append(sam_to_genotype.GenotypeTestRec('*', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['G']))
	test_list.append(sam_to_genotype.GenotypeTestRec('*', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['A']))
	CheckTrue(genotype.CheckEq(test_list))

	# Test adding at the middel.
	genotype.AddInsert(10, 1, 'G')
	genotype.AddInsert(10, 2, 'G')
	test_list.insert(3, sam_to_genotype.GenotypeTestRec('*', 9, 0))
	test_list.insert(4, sam_to_genotype.GenotypeTestRec('*', 9, sam_to_genotype.NUCLEOTIDE_FLAGS['G']))
	test_list.insert(5, sam_to_genotype.GenotypeTestRec('*', 9, sam_to_genotype.NUCLEOTIDE_FLAGS['G']))
	CheckTrue(genotype.CheckEq(test_list))

	print('Add/Insert resutl:')
	PrintAndDump(genotype, 'test_genotype_list_add_insert.gtp')

def TestGenotypeListEmptyInsert():
	genotype = sam_to_genotype.GenotypeList()

	# Test normal sequence adding.
	genotype.Add(7, 'A')
	genotype.Add(8, 'G')
	genotype.Add(9, 'C')
	genotype.Add(10, 'T')
	genotype.Add(11, 'T')
	test_list = [
		sam_to_genotype.GenotypeTestRec(' ', 7, sam_to_genotype.NUCLEOTIDE_FLAGS['A']),
		sam_to_genotype.GenotypeTestRec(' ', 8, sam_to_genotype.NUCLEOTIDE_FLAGS['G']),
		sam_to_genotype.GenotypeTestRec(' ', 9, sam_to_genotype.NUCLEOTIDE_FLAGS['C']),
		sam_to_genotype.GenotypeTestRec(' ', 10, sam_to_genotype.NUCLEOTIDE_FLAGS['T']),
		sam_to_genotype.GenotypeTestRec(' ', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['T'])
	]
	CheckTrue(genotype.CheckEq(test_list))

	print('Emtpy insert list result:')
	PrintAndDump(genotype, 'test_genotype_list_empty_insert.gtp')

def TestGenotypeListEmptyGenotype():
	genotype = sam_to_genotype.GenotypeList()

	# Test insert sequence adding at the end.
	genotype.AddInsert(12, 0, 'C')
	genotype.AddInsert(12, 1, 'G')
	genotype.AddInsert(12, 2, 'A')
	test_list = [
		sam_to_genotype.GenotypeTestRec('*', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['C']),
		sam_to_genotype.GenotypeTestRec('*', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['G']),
		sam_to_genotype.GenotypeTestRec('*', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['A'])
	]
	CheckTrue(genotype.CheckEq(test_list))

	print('Empty genotype list result:')
	PrintAndDump(genotype, 'test_genotype_list_empty_genotype.gtp')

def TestGenotypeListMultiAllelicGenotype():
	genotype = sam_to_genotype.GenotypeList()

	# Test normal sequence adding.
	genotype.Add(7, 'A')
	genotype.Add(8, 'G')
	genotype.Add(9, 'C')
	genotype.Add(10, 'T')
	genotype.Add(11, 'T')
	test_list = [
		sam_to_genotype.GenotypeTestRec(' ', 7, sam_to_genotype.NUCLEOTIDE_FLAGS['A']),
		sam_to_genotype.GenotypeTestRec(' ', 8, sam_to_genotype.NUCLEOTIDE_FLAGS['G']),
		sam_to_genotype.GenotypeTestRec(' ', 9, sam_to_genotype.NUCLEOTIDE_FLAGS['C']),
		sam_to_genotype.GenotypeTestRec(' ', 10, sam_to_genotype.NUCLEOTIDE_FLAGS['T']),
		sam_to_genotype.GenotypeTestRec(' ', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['T'])
	]
	CheckTrue(genotype.CheckEq(test_list))

	genotype.Add(8, 'T')
	test_list[1] = sam_to_genotype.GenotypeTestRec(' ', 8, sam_to_genotype.NUCLEOTIDE_FLAGS['G'] | sam_to_genotype.NUCLEOTIDE_FLAGS['T'])
	CheckTrue(genotype.CheckEq(test_list))

	genotype.Add(8, 'T')
	CheckTrue(genotype.CheckEq(test_list))

	print('Multi-Allelic list result:')
	PrintAndDump(genotype, 'test_genotype_list_multi_allelic.gtp')


def TestGenotypeListFlush():
	genotype = sam_to_genotype.GenotypeList()
	genotype.Add(7, 'A')
	genotype.Add(8, 'G')
	genotype.Add(9, 'C')
	genotype.Add(10, 'T')
	genotype.Add(11, 'T')
	genotype.Add(8, 'T')
	genotype.Add(8, 'A')
	test_list = [
		sam_to_genotype.GenotypeTestRec(' ', 7, sam_to_genotype.NUCLEOTIDE_FLAGS['A']),
		sam_to_genotype.GenotypeTestRec(' ', 8, sam_to_genotype.NUCLEOTIDE_FLAGS['A'] | sam_to_genotype.NUCLEOTIDE_FLAGS['G'] | sam_to_genotype.NUCLEOTIDE_FLAGS['T']),
		sam_to_genotype.GenotypeTestRec(' ', 9, sam_to_genotype.NUCLEOTIDE_FLAGS['C']),
		sam_to_genotype.GenotypeTestRec(' ', 10, sam_to_genotype.NUCLEOTIDE_FLAGS['T']),
		sam_to_genotype.GenotypeTestRec(' ', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['T'])
	]
	CheckTrue(genotype.CheckEq(test_list))
	fl = open('dmp.gtp', 'w')
	genotype.Dump(fl)
	fl.close()
	return

	with open('dmp.gtp', 'a') as f:
		genotype.Flush(8, f)
	test_list = [
		sam_to_genotype.GenotypeTestRec(' ', 8, sam_to_genotype.NUCLEOTIDE_FLAGS['A'] | sam_to_genotype.NUCLEOTIDE_FLAGS['G'] | sam_to_genotype.NUCLEOTIDE_FLAGS['T']),
		sam_to_genotype.GenotypeTestRec(' ', 9, sam_to_genotype.NUCLEOTIDE_FLAGS['C']),
		sam_to_genotype.GenotypeTestRec(' ', 10, sam_to_genotype.NUCLEOTIDE_FLAGS['T']),
		sam_to_genotype.GenotypeTestRec(' ', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['T'])
	]
	CheckTrue(genotype.CheckEq(test_list))

	with open('dmp.gtp', 'a') as f:
		genotype.Flush(11, f)
	test_list = [
		sam_to_genotype.GenotypeTestRec(' ', 11, sam_to_genotype.NUCLEOTIDE_FLAGS['T'])
	]
	CheckTrue(genotype.CheckEq(test_list))

	with open('dmp.gtp', 'a') as f:
		genotype.Dump(f)
	CheckTrue(genotype.CheckEq(test_list))


def TestSAMToGenotype():
	TestGenotypeListAddInsert()
	TestGenotypeListEmptyInsert()
	TestGenotypeListEmptyGenotype()
	TestGenotypeListMultiAllelicGenotype()
	TestGenotypeListFlush()



def main():
	TestSortedSAMFile()
	RunIntervalsTests()
	TestSingleChromosomeFASTAFile()
	TestSAMToGenotype()

	global is_finised
	if is_failed:
		print('  Some test FAILED!')
	else:
		print('---- All tests PASSED! ----')


if __name__ == '__main__':
	main()
