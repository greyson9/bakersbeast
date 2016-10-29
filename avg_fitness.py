import csv
import collections
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from operator import itemgetter
mydict1 = collections.defaultdict(list)
mydictc = collections.defaultdict(list)

with open('20161028_2248_e1_c0_aa_hd1.csv', mode='r') as infile:
    reader = csv.reader(infile)

    for row in reader:
        if row[0] == '-1' or row[0] == '0' or row[0] =='pos':
            continue
        elif row[1] == 'STOP':
            continue
        else:
            mydict1[float(row[0])].append(float(row[5]))
            mydictc[float(row[0])].append(float(row[10]))
average_day1=[]
average_control=[]

for key in mydict1:
    average_day1.append([key,np.median(mydict1[key])])
    average_control.append([key,np.median(mydictc[key])])
average_day1=sorted(average_day1)
average_control=sorted(average_control)
List1=[]
for i in range(0,75):
    List1.append(round(average_day1[i][1]-average_control[i][1],6))
print List1
