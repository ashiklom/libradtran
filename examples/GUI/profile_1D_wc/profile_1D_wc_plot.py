# 
# This file is part of libRadtran.
# Copyright (c) 2010 by Arve Kylling.
# 
# libRadtran is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# libRadtran is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with libRadtran.  If not, see <http://www.gnu.org/licenses/>.

from matplotlib import use
use('WXAgg')
import pylab as plt
import numpy as np
from matplotlib.legend import Legend 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def my_legend(ax,pl_list, text,x,y,size,which):
    l = Legend(ax, pl_list[which-1:which], (text,), loc=(x,y)) 
    #ltext  = la.get_texts()  # all the text.Text instance in the legend
    #plt.setp(ltext, fontsize='medium', linespacing=1)    # the legend text fontsize
    l.draw_frame(False)           # don't draw the legend frame
    ax.add_artist(l) 

plt.figure(figsize=(8,5))

ax = plt.subplot(111)
fil = './profile_1D_wc.OUT'
data = np.loadtxt(fil)
y       = data[:,1]


pl_list = [] 
color   = ['k', 'r', 'g', 'b', 'b--', 'g--' ]
i       = 0
for col in [2, 3, 4, 5, 6, 7]:
    x       = data[:,col]
    pl, = ax.plot(x, y, color[i])
    pl_list.append(pl) 
    i = i + 1


my_legend(ax, pl_list, r'Global', 0.23,0.92, 'large', 1 )
my_legend(ax, pl_list, r'Direct', 0.23,0.87, 'large', 2 )
my_legend(ax, pl_list, r'Diffuse up', 0.23,0.82, 'large', 3 )
my_legend(ax, pl_list, r'Diffuse down', 0.23,0.77, 'large', 4 )
my_legend(ax, pl_list, r'Nadir radiance', 0.23,0.72, 'large', 5 )
my_legend(ax, pl_list, r'Zenith radiance', 0.23,0.67, 'large', 6 )


plt.xlim([0,700])
plt.ylim([0,10])

plt.xlabel(r"Irradiance and radiance (W/(m$^2$ nm))", fontsize = 12)


plt.ylabel(r"Altitude (km)", fontsize = 12)

plt.show()
