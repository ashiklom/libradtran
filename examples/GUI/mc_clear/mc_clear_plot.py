from pylab import *
from matplotlib.legend import Legend 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def my_legend(ax,pl_list, text,x,y,size,which):
    l = Legend(ax, pl_list[which-1:which], (text,), loc=(x,y)) 
    l.draw_frame(False)           # don't draw the legend frame
    ax.add_artist(l) 

rad_clear=loadtxt('mc_clear.rad')
rad_opac=loadtxt('mc_opac.rad')


figure(1, figsize=(7,9))

labels=['I', 'Q', 'U']

pl_list = [] 

for i in range(3):

    ind_0=range(i,shape(rad_clear)[0],8)
    ind_180=range(i+4,shape(rad_clear)[0],8)
    
    ax = subplot(3,1,i+1)
    ax.plot(270-rad_clear[ind_0,2]-90, rad_clear[ind_0,7], 'b')
    pl, = ax.plot(rad_clear[ind_180,2]-90-90, rad_clear[ind_180,7], 'b')
    if i==0:
        pl_list.append(pl) 


    ax.plot(270-rad_opac[ind_0,2]-90, rad_opac[ind_0,7], 'r')
    pl, = ax.plot(rad_opac[ind_180,2]-90-90, rad_opac[ind_180,7], 'r')
    if i==0:
        pl_list.append(pl) 

    xlim([-90,90])
    ylabel(labels[i])
    if i==0:
        title('Normalized Reflectance')
        my_legend(ax, pl_list, r'Rayleigh only', 0.4,0.77, 'large', 1 )
        my_legend(ax, pl_list, r'Rayleigh and aerosol', 0.4,0.67, 'large', 2 )

    elif i==2: 
        xlabel(r"Viewing zenith angle $\theta$ (degree)", fontsize = 12)




subplots_adjust(left=0.15)



show()
