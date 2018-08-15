from matplotlib import use
use('TkAgg')
#use('pdf')
from matplotlib.patches import Patch
import numpy as np
from matplotlib.colors import colorConverter
import math
import matplotlib.pyplot as plt
from matplotlib.legend import Legend 
from matplotlib.legend import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def my_legend(ax,pl_list, text,x,y,size,which, plt):
    l = Legend(ax, pl_list[which-1:which], (text,), loc=(x,y)) 
    ltext  = l.get_texts()  # all the text.Text instance in the legend
    #plt.setp(ltext, fontsize='medium', linespacing=1)    # the legend text fontsize
    plt.setp(ltext, fontsize=size, linespacing=1)    # the legend text fontsize
    l.draw_frame(False)           # don't draw the legend frame
    ax.add_artist(l) 

def read_two_column(fn):

    f  = open(fn,'r')

    x  = []
    y = []
    for line in f:
        if line.find('#') != -1:
            next
        else:
            ary = map(float,line.split())
            x.append(ary[0])
            y.append(ary[1])

    f.close()
    x=np.array(x)
    y=np.array(y)
    return x,y

def read_output(fn):

    f  = open(fn,'r')

    wvl  = []
    edir = []
    uu   = []
    for line in f:
        ary = map(float,line.split())
        wvl.append(ary[0])
        edir.append(ary[1])
        uu.append(ary[3])

    f.close()
    wvl=np.array(wvl)
    edir=np.array(edir)
    uu=np.array(uu)
    return wvl,edir,uu

def conv(wvl,data,slit_function):
    out = 'tmp_conv'
    sav = 'tmp_sav'
    np.savetxt(sav,np.transpose((wvl,data)))
    cmd = home+'/libRadtran/bin/conv '+sav+' '+slit_function+ ' > ' + out
    p   = call(cmd,shell=True,stdin=PIPE,stdout=PIPE)
    data = np.loadtxt(out)
    data = data[:,1]
    return data


fsize=20
fig = plt.figure(figsize=(11,7))
pl_list = [] 


host = fig.add_subplot(111)

fn1 = 'spectrum_flex_noraman_noflu.out'
(wvl,edir,uuo) = read_output(fn1)
col = 'b'
x  = wvl
y  = uuo
col = 'b'
ratnoflu=uuo
pl, = host.plot(x, y, color=col)
pl_list.append(pl) 

fn1 = 'spectrum_flex_noraman_flu.out'
(wvl,edir,uuof) = read_output(fn1)
col = 'r'
x  = wvl
y  = uuof
ratwithflu=uuof
pl, = host.plot(x, y, color=col)
pl_list.append(pl) 

rat =   100*(ratwithflu/ratnoflu-1)
par  = host.twinx()
par.set_ylabel('Fluorescence signal (%)', fontsize = fsize)
x  = wvl
y  = rat
col = 'k'
pl, =par.plot(x, y, color=col)
pl_list.append(pl) 
par.set_ylim(0,10)

ylabel='L [photons/(s nm cm$^2$ sr)]'
ymin= 0.2e+13
ymax= 1.0e+13
ymaL= 0.2e+13
ymiL= 0.5e+12
xmin=677
xmax=697
xmaL= 5
xmiL= 1
txt = '0.3 nm resolution'
my_legend(host, pl_list, txt, 0.02,0.87, fsize, 1, plt )
txt = '0.3 nm resolution, fluorescence included'
my_legend(host, pl_list, txt, 0.02,0.92, fsize, 2,plt )
txt = 'Fluorescence signal'
my_legend(host, pl_list, txt, 0.02,0.8, fsize, 3,plt )
ymajorLocator   = MultipleLocator(2)
par.yaxis.set_major_locator(ymajorLocator)
yminorLocator   = MultipleLocator(0.2)
par.yaxis.set_minor_locator(yminorLocator)

rcParams.update({'font.size':fsize})

host.set_xlim(xmin,xmax)
host.set_ylim(ymin,ymax)
host.set_xlabel(r"Wavelength [nm]", fontsize = fsize)
host.set_ylabel(ylabel, fontsize = fsize)
xminorLocator   = MultipleLocator(xmiL)
xmajorLocator   = MultipleLocator(xmaL)
host.xaxis.set_minor_locator(xminorLocator)
host.xaxis.set_major_locator(xmajorLocator)

ymajorLocator   = MultipleLocator(ymaL)
host.yaxis.set_major_locator(ymajorLocator)
yminorLocator   = MultipleLocator(ymiL)
host.yaxis.set_minor_locator(yminorLocator)

plt.show()
exit()



