from matplotlib import use
use('WXAgg')
import pylab as plt
import numpy as np

plt.figure(figsize=(8,5))

ax = plt.subplot(111)
fil = './spectrum_solar_sbdart.OUT'
data = np.loadtxt(fil)
y       = data[:,1]
x       = data[:,0]

pl_list = [] 
pl, = ax.plot(x,y,'b')
pl_list.append(pl) 

y       = data[:,2]
pl, = ax.plot(x,y,'g')
pl_list.append(pl) 

y       = data[:,3]
pl, = ax.plot(x,y,'r')
pl_list.append(pl) 

plt.xlim([280,4000])
plt.ylim([0,2000])

plt.ylabel(r"Irradiance (mW/(m$^2$ nm))", fontsize = 12)

from matplotlib.legend import Legend 
l0 = Legend(ax, pl_list[0:1], ('Global (direct+diffuse downward)',), loc=(0.3,0.85)) 
#ltext  = l0.get_texts()  # all the text.Text instance in the legend
#plt.setp(ltext, fontsize='small', linespacing=0)    # the legend text fontsize
l0.draw_frame(False)           # don't draw the legend frame
ax.add_artist(l0) 
l0 = Legend(ax, pl_list[1:2], ('Direct',), loc=(0.3,0.75)) 
#ltext  = l0.get_texts()  # all the text.Text instance in the legend
#plt.setp(ltext, fontsize='small', linespacing=0)    # the legend text fontsize
l0.draw_frame(False)           # don't draw the legend frame
ax.add_artist(l0) 
l0 = Legend(ax, pl_list[2:3], ('Upward',), loc=(0.3,0.65)) 
#ltext  = l0.get_texts()  # all the text.Text instance in the legend
#plt.setp(ltext, fontsize='small', linespacing=0)    # the legend text fontsize
l0.draw_frame(False)           # don't draw the legend frame
ax.add_artist(l0) 


plt.xlabel(r"Wavelength (nm)", fontsize = 12)

plt.show()
