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

pl_list = [] 

fil = './profile_1D_wc_heating.OUT'
data = np.loadtxt(fil)
y       = data[:,1]
heating = data[:,3]
fil = './profile_1D_wc_cooling.OUT'
data = np.loadtxt(fil)
y       = data[:,1]
cooling = data[:,3]

net     = heating + cooling

pl, = ax.plot(heating, y, 'r')
pl_list.append(pl) 
pl, = ax.plot(cooling, y, 'b')
pl_list.append(pl) 
pl, = ax.plot(net, y, 'k')
pl_list.append(pl) 

my_legend(ax, pl_list, r'Heating', 0.23,0.72, 'large', 1 )
my_legend(ax, pl_list, r'Cooling', 0.23,0.67, 'large', 2 )
my_legend(ax, pl_list, r'Net', 0.23,0.62, 'large', 3 )


plt.xlim([-50,10])
plt.ylim([0,10])

plt.xlabel(r"Heating/cooling rate (K/day))", fontsize = 12)


plt.ylabel(r"Altitude (km)", fontsize = 12)

plt.show()
