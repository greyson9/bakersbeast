import numpy as np 
import seaborn as sns
import pandas as pd 
from matplotlib import pyplot as plt
from matplotlib import mlab


# generate random input data
a = np.random.rand(1,1000)*2
df = pd.DataFrame(a)

df.to_hdf('test.hf5','w')

a = np.random.rand(1,1000)*2
df = pd.DataFrame(a)

df.to_hdf('test1.hf5','w')

# import data
data = pd.read_hdf('test.hf5')
control = pd.read_hdf('test1.hf5')

a = data.iloc[0]
b = control.iloc[0]



# scatter plot of data along diagonal
plt.figure(figsize=(6,6))

lims = [-5,5]

plt.plot(a, b, 'ro', lims, lims, '--k')
plt.axis([0,2,0,2])
plt.xlabel('no drug')
plt.ylabel('drug')
plt.show()

plt.clf()




# histogram and distribution of ratios
cutoff = 3

ratios = np.log(a / b)
mean = np.mean(ratios)
sigma = np.std(ratios)
x = np.linspace(min(ratios),max(ratios),100)
axes = plt.gca()
plt.hist(ratios, normed=True, bins=20)
plt.plot(x, mlab.normpdf(x, mean, sigma), '-k')
plt.plot((-sigma*cutoff, -sigma*cutoff), (axes.get_ylim()), '--r')
plt.plot((sigma*cutoff, sigma*cutoff), (axes.get_ylim()), '--r')
plt.xlim(min(ratios), max(ratios))
plt.xlabel('ln( drug / no drug )')
plt.ylabel('frequency (%)')
plt.show()

plt.clf()



# show stdev of the fitness ratio of each mutant
sd = []
for mut in ratios:
	sd.append(abs(mut-mean)/sigma)

plt.plot(sd, 'ro')
axes = plt.gca()
plt.plot(axes.get_xlim(),(cutoff,cutoff), '--k')

for index, mut in enumerate(sd):
	if mut >= cutoff:
		plt.annotate(index, xy=(index, mut))

plt.show()

plt.clf()

