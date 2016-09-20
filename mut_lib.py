import cPickle as pic
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
f1 = open("allele_dic.pkl", "rU")
f2 = open("aminotonumber.pkl")
f3 = open("translate.pkl")
allele_dic = pic.load(f1)
amino_to_num = pic.load(f2)
translate = pic.load(f3)
print(amino_to_num)
vals = allele_dic.values()
vals = [string for val in vals for string in val]
#print(vals)
codon_coverage = {}
for val in vals:
    t = val.split("_")
    if int(t[0]) not in codon_coverage:
        codon_coverage[int(t[0])] = [t[1]]
    else:
        codon_coverage[int(t[0])].append(t[1])
#print(coverage)
#for k in codon_coverage:
    #print("Total codon coverage in residue " + str(k) + " is " + str(100 * len(set(codon_coverage[k])) / 64) + "%")

aa_redundancy = {}
codon_redundancy = {}
for k in codon_coverage:
    aas = [translate[cod.replace('T','U')] for cod in codon_coverage[k]]
    aaredun = {}
    for aa in aas:
        if aa in aaredun:
            aaredun[aa] += 1
        else:
            aaredun[aa] = 1
    codredun = {}
    for cod in codon_coverage[k]:
        if cod in codredun:
            codredun[cod] += 1
        else:
            codredun[cod] = 1
    aa_redundancy[k] = aaredun
    codon_redundancy[k] = codredun
#print aa_redundancy
#print codon_redundancy
aa_df = pd.DataFrame.from_dict(aa_redundancy)
aa_df.fillna(0, inplace=True)
aa_df.sort_index()
#print(aa_df)
#sb.heatmap(aa_df)
#c_df = pd.DataFrame.from_dict(codon_redundancy)
c_df.fillna(0,inplace=True)
c_df.sort_index()
sb.heatmap(c_df)
plt.show()
#print(aa_redundancy['12'])