import struct
import sys
import time



def ReadRName(f):
	rname = ''
	ch = f.read(1)
	while ch != '\0':
		rname += ch
		ch = f.read(1)
	return rname


def ReadUint32(f):
	res = ord(f.read(1))
	res |= ord(f.read(1)) << 8
	res |= ord(f.read(1)) << 16
	res |= ord(f.read(1)) << 24
	return res


def GetNucleotides(nucleotides):
	nucleotides_str = ''
	if nucleotides & 0x01:
		nucleotides_str += 'A'
	if nucleotides & 0x02:
		nucleotides_str += 'C'
	if nucleotides & 0x04:
		nucleotides_str += 'G'
	if nucleotides & 0x08:
		nucleotides_str += 'T'
	if nucleotides & 0x10:
		nucleotides_str += 'N'
	return nucleotides_str


def TestFile(path):
	f = open(sys.argv[1], 'rb')
	f.seek(0, 2)
	f_sz = f.tell()
	f.seek(0, 0)
	print('file size : ' + str(f_sz))
	while f.tell() < f_sz:
		rname = ReadRName(f)
		num_genotypes = ReadUint32(f)
		num_inserts = ReadUint32(f)
		num_recs = num_genotypes + num_inserts
		print('>>> ' + rname + '    ' + str(num_genotypes) + '  ' + str(num_inserts))
		cur_pos = 0
		for i in range(num_recs):
			op = f.read(1)
			pos = ReadUint32(f)
			nucleotides = ReadUint32(f)
			if pos < cur_pos:
				print('Error at record #' + str(i) + '!')
				print('          >>> ' + rname)
				print('     cur_pos: ' + str(cur_pos))
				print('         op :`' + op + '\'')
				print('         pos:' + str(pos))
				print(' nucleotides:' + GetNucleotides(nucleotides))
				f.close()
				return
			cur_pos = pos
	f.close()


def main():
	print(time.ctime(time.time()))
	sys.argv.append('SRR1982742-sorted.sam-4.gtp')
	TestFile(sys.argv[1])
	print(time.ctime(time.time()))
	print('Done!')


if __name__ == '__main__':
	main()

