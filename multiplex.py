import gzip
import shutil

file = "ex.fastq.gz"

indices = {"TTGACT", "GGAACT", "TGACAT", "GGACGG", "CTCTAC", "GCGGAC"}

indexed_file = open("indexed.fastq", "w")

with gzip.open(file, "rb") as f: #need gzip.open for fastq file
	count = -1
	for line in f:
		print line
		count += 1
		if (count % 4 == 0) and (line[45:51] in indices) == True:
			indexed_file.write(line)
			indexed_file.write(f.next())
			indexed_file.write(f.next())
			indexed_file.write(f.next())
			count += 3

with open('indexed.fastq', 'rb') as f_in, gzip.open('indexed.fastq.gz', 'wb') as f_out:
    f_out.write(f_in.read())