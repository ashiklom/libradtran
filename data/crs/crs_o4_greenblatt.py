"""
Calculate the O2-O2 collision pair absorption cross section at T=296K
based on Greenblatt et al (1990).

@ARTICLE{Greenblatt1990,
    AUTHOR  = "Gary D. Greenblatt and John J. Orlando and James
                  B. Burkholder and A. R. Ravishankara",
    TITLE   = "{Absorption measurements of oxygen between 330 and 1140
                  nm}",
    JOURNAL = JGR,
    VOLUME  = "95",
    YEAR    = "1990",
    PAGES   = "18,577-18,582"}

"""

import os
from subprocess import Popen,PIPE, STDOUT

# Wavelength (nm), cross section (1.e-46 cm5 molec-2), fwhm
crs = [[ 343.4,  1.20,  4.2],
       [ 360.5,  4.10,  4.8],
       [ 380.2,  2.40,  4.4],
       [ 446.7,  0.57,  5.6],
       [ 477.3,  6.30,  6.2],
       [ 532.2,  1.00, 10.2],
       [ 577.2, 11.00, 11.6],
       [ 630.0,  7.20, 13.8],
       [1065.2, 12.00, 27.3]]

resolution = 0.1  # Resolution of resultant cross section in nm

d = {}
for band in crs:
#    print band
    wvl  = float(band[0])
    crs  = float(band[1])
    fwhm = float(band[2])
    cmd = '../../bin/make_slitfunction -t 3 -n 20 -f '+str(fwhm)+' -r '+str(resolution) 
    p      = Popen(cmd,shell=True,stdin=PIPE,stdout=PIPE)
    data   = p.stdout.read()
    data   = data.split('\n')
    data   = data[0:len(data)-1]  # rm blank line at end of output

    for line in data:
        (rwvl,rcrs) = line.split()
        rwvl = ('{0:12.6f}'.format(wvl+float(rwvl)))
        if rwvl in d:
            d[rwvl] = d[rwvl] + crs*float(rcrs)
        else:
            d[rwvl] = crs*float(rcrs)
#        print ('{0:12.6f} {1:12.6e}'.format(wvl+float(rwvl),crs*float(rcrs)))
        

keys = d.keys()
keys.sort()
for key in keys:
    print key, d[key]


