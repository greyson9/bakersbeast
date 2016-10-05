import cPickle as pickle
from matplotlib import pyplot as plt
import numpy as np

# key: barcode sequence   value: (codon number, codon)
allele_dic = pickle.load(open('allele_dic_with_WT.pkl','rb'))

# key: amino acid   value: number
amino_to_number = pickle.load(open('aminotonumber.pkl', 'rb'))

# key: codon   value: amino acid
translate = pickle.load(open('translate.pkl', 'rb'))

# list to store dictionary of amino_to_number
ub = [{k:0 for k in amino_to_number} for k in range(0,77)]

for barcode in allele_dic:
	codon_num, codon = allele_dic[barcode]
	if codon_num != 0:
		codon = codon.replace('T', 'U')
		amino = translate[codon]

		count_dic = ub[codon_num-1]
		count_dic[amino] += 1


a = np.zeros((21,77))

for i in range(len(ub)):
	for key in ub[i]:
		a[amino_to_number[key],i] = ub[i][key]

print a

plt.figure(figsize=(12,6))
plt.imshow(a, interpolation='nearest', aspect='auto')
plt.colorbar()
rows = ["" for k in range(21)]
for amino in amino_to_number:
	rows[amino_to_number[amino]] = amino

plt.yticks(range(21), rows)
plt.xlabel('residue #')
plt.ylabel('residue identity')
plt.savefig('doug_figure.png', bbox_inches='tight')

pickle.dump(a, open('doug_array.pkl', 'w'))

###################################################

frequencies = [0 for k in range(77)]

b = a.sum(axis=0)

plt.clf()
plt.plot(b)
plt.xlabel('residue #')
plt.ylabel('# of mutants in library')
plt.xlim(0,76)
plt.savefig('doug_figure2.png', bbox_inches='tight')

