from matplotlib import use
use('WXAgg')
import pylab as plt
import numpy as np

plt.figure(figsize=(8,5))

ax = plt.subplot(111)
fil = './spectrum_sciamachy.OUT'
data = np.loadtxt(fil)
y       = data[:,1]
x       = data[:,0]

pl_list = [] 
pl, = ax.plot(x,y,'r')
pl_list.append(pl) 

#plt.xlim([425,450])
#plt.ylim([0,2000])

plt.ylabel(r"Nadir reflectivity", fontsize = 12)

plt.xlabel(r"Wavelength (nm)", fontsize = 12)

#plt.show()

plt.savefig('spectrum_sciamachy.png')
