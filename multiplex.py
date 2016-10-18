import gzip
import shutil

file = "Undetermined_S0_R1_001.fastq.gz"

indices = {"TTGACT", "GGAACT", "TGACAT", "GGACGG", "CTCTAC", "GCGGAC"}

indexed_file = open("indexed.fastq", "w")


def hamming_dist(s1, s2): 
	assert len(s1) == len(s2)
	return sum(c1 != c2 for c1, c2 in zip(s1, s2))

with gzip.open(file, "rb") as f: #need gzip.open for fastq file
	count = -1
	for line in f:
		# print line
		count += 1
		if (count % 4 == 0):
			close = False
			for ind in indices:
				close = close or (hamming_dist(line[45:51], ind) < 2)
			if close:
				indexed_file.write(line)
				indexed_file.write(f.next())
				indexed_file.write(f.next())
				indexed_file.write(f.next())
				count += 3

with open('indexed.fastq', 'rb') as f_in, gzip.open('indexed.fastq.gz', 'wb') as f_out:
    f_out.write(f_in.read())