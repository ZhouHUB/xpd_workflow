"""
Copyright (c) 2014 Brookhaven National Laboratory All rights reserved.
Use is subject to license terms and conditions.

@author: Christopher J. Wright"""
__author__ = 'Christopher J. Wright'

import numpy as np
import matplotlib.pyplot as plt

import IO
import Calculate


poni_file = '/media/Seagate_Backup_Plus_Drive/Mar2014/Ni/17_21_24_NiSTD_300K-00002.poni'
geometry_object = IO.loadgeometry(poni_file)

# load image
xray_image = IO.loadimage(
    '/media/Seagate_Backup_Plus_Drive/Mar2014/Ni/17_21_24_NiSTD_300K-00002.tif')

# generate radial array based on image shape
# TODO: use a more Q/angle based system
Radi = geometry_object.rArray(xray_image.shape)
angles = geometry_object.chiArray(xray_image.shape)

# create integer array of radi by division by the pixel resolution
roundR = np.around(Radi / geometry_object.pixel1).astype(int)

# define the maximum radius for creation of numpy arrays
# make maxr into an integer so it can be used as a counter
maxr = int(np.ceil(np.amax(roundR)))

pols = []
print 'start opts'
for r in range(0, maxr, 100):
    pols.append(
        Calculate.optomize_polarization(r, xray_image, geometry_object, roundR, angles))
pols = np.array(pols)
plt.plot(pols)
plt.show()
