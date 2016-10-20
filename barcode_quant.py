# barcode_quant.py

# To read csv, use command df = pd.read_csv(fn, index_col=['pos','aa'/'cod'])

# e0 is day 1, e1 is day 2, and c0 is the control. These are the *expts*
# *cqs* are the calculated quantities: slope (relative fitness), r (correlation coefficient),
# p (not clear on its meaning in this context), stde (standard error)
#
# To generate a heatmap for a given quantity, use the following code:
# df_*expt* = df['*expt*_	*cq*'].reset_index().pivot('pos', 'aa')
# df_*expt*.columns = df_*expt*.columns.droplevel()
# df_e0 = df_e0.T
# sns.heatmap(df_*expt*, vmin=-0.6, vmax=0.6) # vmin/vmax are the min and max of the colorbar
# plt.xticks(rotation=90)
# plt.show()
import gzip
import itertools
import sys
# from scipy.stats import linregress
import cPickle as pic
import pandas as pd
from Bio import Seq
# import statsmodels.api as sm
from scipy.stats import linregress
import numpy as np
import seaborn as sns

# create a dictionary associating nucleotide index sequences with indices
INDEX_DAY = {'TTGACT': 'e0t0', 'GGAACT': 'e0t1', 'TGACAT': 'e0t2', 
				'GGACGG': 'e1t0', 'CTCTAC': 'e1t1', 'GCGGAC': 'e1t2', 
				'ATTCCG': 'c0t0', 'AGCTAG': 'c0t1', 'GTATAG': 'c0t2'}
DAY_INDEX = {0: 'TTGACT', 1: 'GGAACT', 2: 'TGACAT', 3: 'GGACGG', 4: 'CTCTAC', 5: 'GCGGAC'}
# create a dictionary associating indices with timepoints
INDEX_TIMES = {0: 0.0, 1: 2.02, 2: 3.69, 3: 0.0, 4: 1.44, 5: 3.67}
X1 = [0.0, 2.02, 3.69]
X2 = [0.0, 1.44, 3.67]
X3 = [0.0, 2.02, 3.50]

# filename is specified in the command line
fn = sys.argv[1]

f1 = open("allele_dic_with_WT.pkl", "rU")
f2 = open("aminotonumber.pkl")
f3 = open("translate.pkl")
allele_dic = pic.load(f1)
amino_to_num = pic.load(f2)
translate = pic.load(f3)


def hamming_dist_2(s1, s2):
	assert len(s1) == len(s2)
	dist = 0
	c1 = list(s1)
	c2 = list(s2)
	for i in range(0, len(s2)):
		if dist > 2:
			break
		else:
			if c1[i] != c2[i]:
				dist += 1
	return dist


def count_barcodes(filename):
	# the dictionary of dictionaries
	# first level keys are timepoints
	# second level of keys are barcodes
	# second level of values are number of reads
	samples = {}
	# read in the compressed file line by line
	with gzip.open(filename, 'rb') as f:
		# for each pair of lines in the file
		for l1, l2 in itertools.izip_longest(*[f]*2):
			# strip the first line; it should be the index
			index = l1.strip()
			# strip the second line; it contains the barcode and quality score
			v = l2.strip()
			# the barcode is all but the last two characters of the second line
			barcode = v[:18]
			seq = str(Seq.Seq(barcode, Seq.Alphabet.SingleLetterAlphabet()).reverse_complement())
			# the score is the last two characters of the second line
			# score = v[-2:]
			# first, check to see if the barcode has been added
			# if not, create a k,v pair with k=timepoint and v=barcode
			if (seq, index) not in samples:
				samples[(seq, index)] = 1.0
			# otherwise, the barcode has been added, 
			# so check to see if barcode has been added
			else:
				samples[(seq, index)] += 1.0
		return samples

# this method counts barcodes for each timepoints
def count_barcodes_hamming(filename):
	# the dictionary of dictionaries
	# first level keys are timepoints
	# second level of keys are barcodes
	# second level of values are number of reads
	samples = {}
	# read in the compressed file line by line
	with gzip.open(filename, 'rb') as f:
		# for each pair of lines in the file
		for l1, l2 in itertools.izip_longest(*[f]*2):
			# strip the first line; it should be the index
			index = l1.strip()
			# strip the second line; it contains the barcode and quality score
			v = l2.strip()
			# the barcode is all but the last two characters of the second line
			barcode = v[:18]
			if barcode not in allele_dic:
				for k in allele_dic:
					if hamming_dist_2(barcode, k) < 3:
						barcode = k
						break
			seq = str(Seq.Seq(barcode, Seq.Alphabet.SingleLetterAlphabet()).reverse_complement())
			# the score is the last two characters of the second line
			# score = v[-2:]
			# first, check to see if the barcode has been added
			# if not, create a k,v pair with k=timepoint and v=barcode
			if (seq, index) not in samples:
				samples[(seq, index)] = 1.0
			# otherwise, the barcode has been added, 
			# so check to see if barcode has been added
			else:
				samples[(seq, index)] += 1.0
		return samples
# now, convert this into codon reads and amino acid reads
# the fitness score is defined as k_mut / k_wt
# k_mut is the slope of the line formed by three timepoints of the mut growth
# this method calculates the fitness score for each barcode
# outline:


def fitness_scores(input):
	df = pd.DataFrame([{'barcode': k[0], 'time': k[1], 'count': input[k]} for k in input])
	pt = df.pivot(index='barcode', columns='time', values='count')
	pt['pos'] = [allele_dic[k][0] if k in allele_dic else -1 
					for k in pt.index]
	pt['codon'] = [allele_dic[k][1].replace('T', 'U') if k in allele_dic 
					else "none" for k in pt.index]
	pt['aa'] = [translate[cod] if cod in translate 
					else "none" for cod in pt['codon']]
	pt.fillna(value=0.1, inplace=True)
	df_cod = pt.groupby(['pos', 'codon']).sum()
	df_cod = df_cod.div(df_cod.iloc[1], axis='columns')
	df_cod = np.log2(df_cod)
	df_cod.columns = [INDEX_DAY[x] for x in df_cod.columns]
	df_cod = pd.merge(
		df_cod,
		pd.DataFrame([linregress(X1, np.asarray(r[1])) for r in 
		df_cod[['e0t0', 'e0t1', 'e0t2']].iterrows()], 
		columns=['e0_slope', 'e0_itcp', 'e0_r', 'e0_p', 'e0_stde'], 
		index=df_cod.index), left_index=True, right_index=True)
	df_cod = pd.merge(
		df_cod,
		pd.DataFrame([linregress(X2, np.asarray(r[1])) for r in 
		df_cod[['e1t0', 'e1t1', 'e1t2']].iterrows()], 
		columns=['e1_slope', 'e1_itcp', 'e1_r', 'e1_p', 'e1_stde'], 
		index=df_cod.index), left_index=True, right_index=True)
	df_cod = pd.merge(
		df_cod,
		pd.DataFrame([linregress(X3, np.asarray(r[1])) for r in 
		df_cod[['c0t0', 'c0t1', 'c0t2']].iterrows()], 
		columns=['c0_slope', 'c0_itcp', 'c0_r', 'c0_p', 'c0_stde'], 
		index=df_cod.index), left_index=True, right_index=True)
	df_aa = pt.groupby(['pos', 'aa']).sum()
	df_aa = df_aa.div(df_aa.iloc[1], axis='columns')
	df_aa = np.log2(df_aa)
	df_aa.columns = [INDEX_DAY[x] for x in df_aa.columns]
	df_aa = pd.merge(
		df_aa,
		pd.DataFrame([linregress(X1, np.asarray(r[1])) for r in 
		df_aa[['e0t0', 'e0t1', 'e0t2']].iterrows()], 
		columns=['e0_slope', 'e0_itcp', 'e0_r', 'e0_p', 'e0_stde'], 
		index=df_aa.index), left_index=True, right_index=True)
	df_aa = pd.merge(
		df_aa,
		pd.DataFrame([linregress(X2, np.asarray(r[1])) for r in 
		df_aa[['e1t0', 'e1t1', 'e1t2']].iterrows()], 
		columns=['e1_slope', 'e1_itcp', 'e1_r', 'e1_p', 'e1_stde'], 
		index=df_aa.index), left_index=True, right_index=True)
	df_aa = pd.merge(
		df_aa,
		pd.DataFrame([linregress(X3, np.asarray(r[1])) for r in 
		df_aa[['c0t0', 'c0t1', 'c0t2']].iterrows()], 
		columns=['c0_slope', 'c0_itcp', 'c0_r', 'c0_p', 'c0_stde'], 
		index=df_aa.index), left_index=True, right_index=True)
	df_cod.to_csv('codons.csv')
	df_aa.to_csv('aas.csv')
	return df_cod, df_aa
cod, aa = fitness_scores(count_barcodes(fn))
