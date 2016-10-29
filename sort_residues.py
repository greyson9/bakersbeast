import csv
import collections
import numpy as np
from operator import itemgetter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
mydict1 = collections.defaultdict(list)
mydict2 = collections.defaultdict(list)
mydictc = collections.defaultdict(list)

with open('aas.csv', mode='r') as infile:
    reader = csv.reader(infile)

    for row in reader:
        if row[0] == '-1' or row[0] == '0' or row[0] =='pos' or row[0]=='77':
            continue
        elif row[1] == 'STOP':
            continue
        else:
            mydict1[float(row[0])].append([row[1],float(row[11])])
            mydict2[float(row[0])].append([row[1],float(row[16])])
            mydictc[float(row[0])].append([row[1],float(row[21])])

polar=['S','T','Q','C','N']
nonpolar=['G','A','V','I','L','P','M']
poscharge=['K','H','R']
negcharge=['E','D']
aromatic=['F','Y','W']

def aa_sort(dict,dict2):
    output=[]
    for key in dict:
        p=[]
        nonp=[]
        pos=[]
        neg=[]
        a=[]
        p1=[]
        nonp1=[]
        pos1=[]
        neg1=[]
        a1=[]
        for item in dict[key]:
            if item[0] in polar:
                p.append(item[1])
            elif item[0] in nonpolar:
                nonp.append(item[1])
            elif item[0] in poscharge:
                pos.append(item[1])
            elif item[0] in negcharge:
                neg.append(item[1])
            elif item[0] in aromatic:
                a.append(item[1])
        for item in dict2[key]:
            if item[0] in polar:
                p1.append(item[1])
            elif item[0] in nonpolar:
                nonp1.append(item[1])
            elif item[0] in poscharge:
                pos1.append(item[1])
            elif item[0] in negcharge:
                neg1.append(item[1])
            elif item[0] in aromatic:
                a1.append(item[1])
        output.append([key,round(np.median(np.array(p)-np.array(p1)),6),round(np.median(np.array(nonp)-np.array(nonp1)),6),round(np.median(np.array(pos)-np.array(pos1)),6),round(np.median(np.array(neg)-np.array(neg1)),6),round(np.median(np.array(a)-np.array(a1)),6)])
    output=sorted(output)
    return output
sort_day1= aa_sort(mydict1,mydictc)
sort_day2= aa_sort(mydict2,mydictc)
# difference_day1={}
# difference_day2={}



# 1 is polar
# 2 is nonpolar
# 3 is pos charge
# 4 is neg charge
# 5 is aromatic

# for i in range(0,75):
#     for j in ['0','1','2','3','4']:
#         print j
#         if j in difference_day1:
#             difference_day1[j].append(round(sort_day1[i][int(j)]-sort_control[i][int(j)],6))
#             difference_day2[j].append(round(sort_day2[i][int(j)]-sort_control[i][int(j)],6))
#         else:
#             difference_day1[j]=round(sort_day1[i][int(j)]-sort_control[i][int(j)],6)
#             difference_day2[j]=round(sort_day2[i][int(j)]-sort_control[i][int(j)],6)
#


for item in sort_day1:
     del item[0]
for item in sort_day2:
     del item[0]
# for item in sort_control:
#     del item[0]
#
dataframe_day1=pd.DataFrame(sort_day1,index=range(2,77),columns=['polar','nonpolar','pos charge','neg charge','aromatic'])
heatmap1=sns.heatmap(pd.DataFrame.transpose(dataframe_day1))
plt.show()
dataframe_day2=pd.DataFrame(sort_day2,index=range(2,77),columns=['polar','nonpolar','pos charge','neg charge','aromatic'])
# dataframe_control=pd.DataFrame(sort_control,index=range(2,77),columns=['polar','nonpolar','pos charge','neg charge','aromatic'])
#
dataframe_day1.to_csv('residue_type_day1.csv')
dataframe_day2.to_csv('residue_type_day2.csv')
# dataframe_control.to_csv('residue_type_control.csv')
