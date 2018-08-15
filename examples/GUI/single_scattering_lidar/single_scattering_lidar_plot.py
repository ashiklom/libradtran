from pylab import *
#from matplotlib import *
#from scipy import *

fontsize=15

def prepol(fontsize):
    rc("axes", linewidth=2.0)
    rc("lines", linewidth=3.0)
    rc("lines", markeredgewidth=2.0)
    ax = gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

data=loadtxt('./single_scattering_lidar.OUT')

semilogx((data[:,1]+data[:,3]),data[:,0],"k",label='I+Q')
semilogx((data[:,1]-data[:,3]),data[:,0],"r",label='I-Q')
semilogx(data[:,2],data[:,0],"k:",label='lidar ratio')

xlabel('Signal strength [photons per range bin]',fontsize=fontsize)
ylabel('range/height [km]',fontsize=fontsize)
title("Lidar simulation with polarisation",fontsize=fontsize)
legend(loc="upper right")
prepol(fontsize)
show()

xlim((1e-1,1e5))
ylim((0,8.2))
