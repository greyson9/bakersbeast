import numpy as np 
import seaborn as sns
import pandas as pd 
from matplotlib import pyplot as plt
from matplotlib import mlab
sns.set()

# import data
data = pd.read_csv('aas.csv')

dataset1 = 'e1'
dataset2 = 'c0'




a = data[dataset1 + '_slope']
b = data[dataset2 + '_slope']



# scatter plot of data along diagonal
plt.figure(figsize=(10,10), dpi=600)

lims = [-5,5]

plt.plot(a, b, 'ro', lims, lims, '--k')
plt.xlabel(dataset1)
plt.ylabel(dataset2)
plt.axis([-2,2,-2,2])
plt.savefig('seq_figures/scatter_' + dataset1 + '_' + dataset2 + '.png')

plt.close()


# histogram and distribution of ratios
cutoff = 3

plt.figure(figsize=(10,10), dpi=600)

ratios = a - b
mean = np.mean(ratios)
sigma = np.std(ratios)
x = np.linspace(min(ratios),max(ratios),100)
axes = plt.gca()
plt.hist(ratios, normed=True, bins=20)
plt.plot(x, mlab.normpdf(x, mean, sigma), '-k')
plt.plot((-sigma*cutoff + mean, -sigma*cutoff + mean), (axes.get_ylim()), '--r')
plt.plot((sigma*cutoff + mean, sigma*cutoff + mean), (axes.get_ylim()), '--r')
plt.xlim(min(ratios), max(ratios))
plt.xlabel(dataset1 + ' - ' + dataset2)
plt.ylabel('')

plt.savefig('seq_figures/hist_' + dataset1 + '_' + dataset2 + '.png')

plt.close()



# show stdev of the fitness ratio of each mutant
sd = []
for mut in ratios:
	sd.append((mut-mean)/sigma)

plt.figure(figsize=(10,10), dpi=600)
plt.plot(sd, 'ro')
axes = plt.gca()
plt.xlabel('mutant')
plt.ylabel('standard deviations from mean')
plt.plot(axes.get_xlim(),(cutoff,cutoff), '--k')
plt.plot(axes.get_xlim(),(-cutoff,-cutoff), '--k')

for index, mut in enumerate(sd):
	if mut >= cutoff or mut <= -cutoff:
		label = str(data['pos'][index]) + str(data['aa'][index])
		plt.annotate(label, xy=(index, mut))

plt.savefig('seq_figures/sigma_' + dataset1 + '_' + dataset2 + '.png')

plt.clf()

