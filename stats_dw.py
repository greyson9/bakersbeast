import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
from scipy.stats import mannwhitneyu as mwu 
import seaborn as sns
from scipy.stats import ttest_ind as ttest 
sns.set()
from random import sample


# IMAGING STATISTICS


# data1 = pd.read_hdf('dumps/Plate000_WellC01_Seq0000.nd2.hf5')
# data2 = pd.read_hdf('dumps/Plate000_WellD01_Seq0023.nd2.hf5')

# bins1 = np.linspace(data1['Distance'].min(), data1['Distance'].max(), 11)
# binned1 = data1.groupby(pd.cut(data1['Distance'], bins1))

# bins2 = np.linspace(data2['Distance'].min(), data2['Distance'].max(), 11)
# binned2 = data2.groupby(pd.cut(data2['Distance'], bins2))

# first = []
# second = []

# for group, df in binned1:
# 	first.append(df['Signal'])

# for group, df in binned2:
# 	second.append(df['Signal'])

# for i in range(len(first)):
# 	t, p = ttest(first[i],second[i])
# 	print i, t, p


# SEQUENCING STATISTICS

data1 = pd.read_csv('aas.csv')


