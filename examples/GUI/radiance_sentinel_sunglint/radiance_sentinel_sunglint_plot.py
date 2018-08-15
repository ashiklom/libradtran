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
import re
import sys
import math
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def read_inp(fname):
	f = open(fname, "r")
	txt = f.read()
	f.close()
	
	pattern = "((umu)|(phi)) (([\s|\d|[-]|[.])*)"
	out = {}

	search = re.search(pattern, txt)
	while search:
		out[search.group(1)] = [float(i) for i in search.group(4).split()]
		txt = txt[search.end():]
		search = re.search(pattern, txt)

	return out


inp_fname = "./radiance_sentinel_sunglint.INP"
umu, phi = read_inp(inp_fname).items()
umu_txt, umus = umu

out_fname = "./radiance_sentinel_sunglint.OUT"
data = list(np.loadtxt(out_fname))

# The reading of radiances is a bit complicated. They all come on one line for each
# wavelength and altitude when "output_user uu" is specified. In addition, the first
# number is for the phi[0], the second for phi[1] and so on. Since we have many umu
# and two phis (0 and 180 to get the fill principal plane for umu between 0 and 90)
# in the input file the following will unpack, sort and get the right viewing angles.
#

# Two phis, each uu for one phi is in every second element of data.
uu1       =  data[1::2]
uu2       =  data[::2]

# Reverse uu1 since it is for phi=0 and we want viewing angle later on.

uu1.reverse()

# Get the complete list of all uu

uu        =  uu2 + uu1

# A solar zenith angle of 30 corresponds to a umu of -30 (sign convention).
# For convenience we change the sign of the umus to be the same as of the solar
# zenith angle.

angles = [(180-math.acos(float(-umu))*180/3.14159) for umu in umus]

# Get angles for phi=0

bangles = list(angles)

# In reverse as for uu1

bangles.reverse()

bangles = [-angle for angle in bangles]


# And get the complete list of all angles

all_angles = angles + bangles

# Next plot

x = all_angles
y = uu


plt.figure(figsize=(8,5))

ax = plt.subplot(111)

pl_list = [] 
pl, = ax.plot(x,y,'b')
pl_list.append(pl) 

plt.xlim([-90,90])
plt.ylim([0,110])

plt.ylabel(r"Radiance", fontsize = 12)
plt.xlabel(r"Viewing nadir angle $\theta$ (degree)", fontsize = 12)

xmajorLocator   = MultipleLocator(20)
ymajorLocator   = MultipleLocator(20)
majorFormatter  = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(5)
yminorLocator   = MultipleLocator(5)

ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
#for the minor ticks, use no labels; default NullFormatter
ax.xaxis.set_minor_locator(xminorLocator)

ax.yaxis.set_major_locator(ymajorLocator)
ax.yaxis.set_minor_locator(yminorLocator)


plt.show()



