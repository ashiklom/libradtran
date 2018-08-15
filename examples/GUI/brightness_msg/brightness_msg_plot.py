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
fil = './brightness_msg.OUT'
data = np.loadtxt(fil)
y       = data[1:,1]
x       = data[1:,0]

pl_list = [] 
pl, = ax.plot(x,y,'b')
pl_list.append(pl) 

plt.xlim([0,1])
plt.ylim([235,265])

plt.ylabel(r"Brightness temperature", fontsize = 12)


plt.xlabel(r"Cosine of viewing nadir angle", fontsize = 12)

plt.show()
