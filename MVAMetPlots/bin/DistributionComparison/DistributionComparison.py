import numpy as np

import matplotlib.pyplot as plt
minUni = 150
maxUni = 160
mu = (maxUni+minUni)/2
sigma = 29.5
samples = 100000


uniformDist = np.random.uniform(minUni,maxUni,samples)
normalDist = np.random.normal(mu,sigma,samples)
plt.ylabel('Frequency distribution', fontsize=17)
plt.xlabel('X', fontsize=17)
count, bins, ignored = plt.hist([uniformDist,normalDist], 50, normed=True, range=[80,130], label=['Uniform','Normal'])
plt.text(0.7*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.3*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$Min_U = %i$''\n'r'$Max_U = %i$''\n'r'$\mu_N = %i$''\n'r'$\sigma_N = %i$'%(minUni,maxUni,mu,sigma),color = 'k',fontsize=16)
plt.legend(loc='best')
plt.savefig('DistributionInitial.png')
plt.clf()
count, bins, ignored = plt.hist([(uniformDist/normalDist),(normalDist/uniformDist)], 50, normed=True, label = ['Uniform/Normal','Normal/Uniform'], range=[0.9,1.1])
plt.clf()
meanUN = (uniformDist/normalDist).mean()
meanNU = (normalDist/uniformDist).mean()
stdUN = (uniformDist/normalDist).std()
stdNU = (normalDist/uniformDist).std()
print((uniformDist/normalDist).mean())
print((normalDist/uniformDist).mean())
print((uniformDist/normalDist).std())
print((normalDist/uniformDist).std())
bincenters = bins[:-1] + (bins[1:]-bins[:-1])/2.
plt.plot(bincenters, count[0], 'r--', label = 'Uniform/Normal')
plt.plot(bincenters, count[1], 'b--', label = 'Normal/Uniform')
plt.legend(loc='best')
plt.ylabel('Frequency distribution', fontsize=17)
plt.xlabel('Ratio', fontsize=17)
plt.text(0.7*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.2*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$Min_U = %i$''\n'r'$Max_U = %i$''\n'r'$\mu_N = %i$''\n'r'$\sigma_N = %i$''\n'r'$\mu_{U/N} = %.3f$''\n'r'$\sigma_{U/N} = %.3f$''\n'r'$\mu_{N/U} = %.3f$''\n'r'$\sigma_{N/U} = %.3f$''\n'%(minUni,maxUni,mu,sigma,meanUN,stdUN,meanNU,stdNU),color = 'k',fontsize=16)
plt.savefig('DistributionFraction.png')
