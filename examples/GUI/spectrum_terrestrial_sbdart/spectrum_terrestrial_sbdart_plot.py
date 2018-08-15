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

plt.figure(figsize=(8,5))

ax = plt.subplot(111)
fil = './spectrum_terrestrial_sbdart.OUT'
data = np.loadtxt(fil)
y       = data[:,2]
x       = data[:,1]

pl_list = [] 
pl, = ax.plot(x,y,'b')
pl_list.append(pl) 

plt.xlim([3000,100])
plt.ylim([0,0.5])

plt.ylabel(r"Downward irradiance (W/(m$^2$ nm))", fontsize = 12)


plt.xlabel(r"Wavenumber (cm$^{-1}$)", fontsize = 12)

plt.show()
