# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 17:31:43 2014

@author: christopher
"""

from __future__ import division
import tkFileDialog  # maybe depends on wether a wrapper is going to be usied
import math
import os
import warnings

from pyFAI import geometry as geo
import numpy as np
import fabio
import matplotlib.pyplot as plt
import scipy
import pyFAI
import scipy.ndimage.filters as snf


warnings.simplefilter("error")
# load image and geometry

def loadgeometry(ponifilepath=None):
    """Loads the detector geometry information from a poni file or from hard coded information
    Parameters
    ----------
    ponifilepath: str
        File path to the .poni file
    Returns
    -------
    a: an object which contains the geometry information for the detector
    """
    if ponifilepath is None:
        distance = 0.204666698799
        poni1 = 0.205911555168
        poni2 = 0.207264947729
        rot1 = 0.0105322341791
        rot2 = 0.0104587423844
        rot3 = -3.27539683464e-08
        spline_file = None
        wavelength = 1.839e-11
        detector = "Perkin"
        # This is the PyFAI geometry stuff
        a = geo.Geometry(dist=distance, poni1=poni1, poni2=poni2, rot1=rot1,
                         rot2=rot2, rot3=rot3, pixel1=200e-6, pixel2=200e-6,
                         splineFile=None, detector=detector,
                         wavelength=wavelength)
    else:
        a = pyFAI.load(ponifilepath)
    return a


def loadimage(filepath=None):
    """
    May be a placeholder for a different fileIO, may need summing capacity
    """
    if filepath is None:
        # simple file dialog if the filepath is not found somewhere in the intervening code
        filepath = tkFileDialog.askopenfilename()
    xray_image = fabio.open(filepath).data
    return xray_image


def raw1d(img, max=None, step=None, position_matrix=None):
    bins = np.arange(0, max + step, step)
    return np.histogram(position_matrix, bins, weights=img)[0]


def mean1d(img, max=None, step=None, position_matrix=None):
    bins = np.arange(0, max + step, step)
    pixels_per_bin = np.array(np.histogram(position_matrix, bins)[0],
                              dtype=float)
    # plt.plot(np.histogram(position_matrix, bins, weights=img)[0]/pixels_per_bin)
    # plt.show()
    return np.histogram(position_matrix, bins, weights=img)[0] / pixels_per_bin


def median1d(img, max=None, step=None, position_matrix=None):
    bins = np.arange(0, max + step, step)
    pixels_per_bin = np.array(np.histogram(position_matrix, bins)[0],
                              dtype=float)
    # plt.plot(np.histogram(position_matrix, bins, weights=img)[0]/pixels_per_bin)
    # plt.show()
    # TODO: may need bleeding edge scipy
    return \
        scipy.stats.binned_statistic(position_matrix, img, statistic='median',
                                     bins=bins)[0]


def ringstd(img, max=None, step=None, position_matrix=None):
    # <x^2>-<x>^2
    average_squared = mean1d(img, max=max, step=step,
                             position_matrix=position_matrix) ** 2
    square_averaged = mean1d(img ** 2, max=max, step=step,
                             position_matrix=position_matrix)
    return np.sqrt(average_squared - square_averaged)


def variance1d(pic, max=None, step=None, position_matrix=None):
    picavg = snf.uniform_filter(pic, 5, mode='wrap')
    pics2 = (pic - picavg) ** 2
    pvar = snf.uniform_filter(pics2, 5, mode='wrap')

    gain = np.divide(pvar, pic)

    gain[np.isnan(gain)] = 0
    gain[np.isinf(gain)] = 0
    gainmedian = np.median(gain)
    var = pic * gainmedian
    bins = np.arange(0, max + step, step)
    # print(np.shape(bins))
    pixels_per_bin = np.array(np.histogram(position_matrix, bins)[0],
                              dtype=float)
    variance = np.histogram(position_matrix, bins, weights=var)[0]
    # print(np.shape(variance))
    # print variance
    # plt.plot(np.sqrt(variance))
    # plt.show()
    return variance / pixels_per_bin


def maskw(m, edge, xray_image, imax, jmax, maxr, roundR, listcounter):
    """
    Takes in a bunch of parameters for masking the image

    Parameters
    ----------
    Stuff ;

    Returns
    -------
    highper, lowper, totalper, old_files, iMask, StatsA


    """
    # old_files=np.zeros((imax,jmax))
    # highper=np.zeros(maxr+1)
    # lowper=np.zeros(maxr+1)
    # StatsA=[]
    # r=0
    # TODO: Is this a bottle neck? If so maybe implement this in weave/C/GPU?
    # while r<maxr+1:
    # Im contains all the pixel intensities for the ring with radius r
    # Im=xray_image[roundR==r]
    # mean = np.mean(Im)
    # std = np.std(Im)
    # thresh = m * std
    # Im1 filters out the pixels that are beyond the threshold
    # Im1=Im[abs(Im - mean) <= thresh]
    # note that the means, and medians hold the integrated intensities
    # StatsA.append([r, np.mean(Im1),np.median(Im1),mean + thresh,mean - thresh,std,mean])
    # r+=1
    # StatsA = np.array(StatsA)

    mean = mean1d(xray_image, max=maxr, step=1, position_matrix=roundR)
    # std=sqrt(Variance1D(xray_image, max=maxr, step=1, position_matrix=roundR))
    std = ringstd(xray_image, max=maxr, step=1, position_matrix=roundR)
    plt.plot(std)
    plt.show()
    threshold = m * std
    lower = mean - threshold
    upper = mean + threshold
    print('start masking')

    # start masking based on too high/too low/not enough pixels
    # print np.min(roundR)
    # TODO: there is a problem here with the indexing
    toolow = xray_image < lower[roundR - 1]
    toohi = xray_image > upper[roundR - 1]
    lowcounter = np.bincount(roundR[toolow], minlength=maxr + 1)
    highcounter = np.bincount(roundR[toohi], minlength=maxr + 1)
    # TODO: maybe read magic number for number of pixels from config file or calculate based on detector info?
    # combine toolow, toohigh, and number of pixels too low
    Mask = toolow | toohi | (listcounter[roundR] < 10)
    # mask edges
    Mask[:, :edge] = 1
    Mask[:, -edge:] = 1
    Mask[:edge, :] = 1
    Mask[-edge:, :] = 1
    # percentages for graphing
    highper = (highcounter[:] / listcounter[:]) * 100
    lowper = (lowcounter[:] / listcounter[:]) * 100
    totalper = highper + lowper
    iMask = np.abs(Mask - 1)
    print 'finished masking'
    return highper, lowper, totalper, Mask, iMask,  # StatsA


def plotmask(maskwoutput=None, maxr=None, xray_image=None,
             position_matrix=None):
    """
    Plots the mask metrics allowing for evaluation of mask quality
    """

    # highper,lowper,totalper, old_files, iMask,StatsA=maskwoutput
    highper, lowper, totalper, Mask, iMask = maskwoutput
    # May need matplotlib help from Tom on this
    plt.ioff()
    plt.clf()
    # plot the graph of how the mask depends on radial distance
    # TODO: use Q instead of r in the final product?
    plt.plot(np.arange(0, maxr + 1, 1), highper[:], '-',
             np.arange(0, maxr + 1, 1), lowper[:], 'o',
             np.arange(0, maxr + 1, 1), totalper[:], '^')
    plt.show()
    # Plot the mask itself
    # plot the overlay of the mask and the image
    mask = plt.imshow(Mask)
    plt.show()
    imagemask = plt.imshow(xray_image * iMask)
    plt.show()
    t = np.arange(0, maxr, 1)
    plt.plot(t, mean1d(xray_image, max=maxr, step=1,
                       position_matrix=position_matrix * iMask), 'b',
             # t, Median1D(xray_image, max=maxr, step=1, position_matrix=position_matrix*iMask), 'r',\
             # t,Median1D(xray_image, max=maxr, step=1, position_matrix=position_matrix*iMask)-Mean1D(xray_image, max=maxr, step=1, position_matrix=position_matrix*iMask)
    )
    plt.show()

    oldstd = np.sqrt(variance1d(xray_image, max=maxr, step=1,
                                position_matrix=position_matrix))
    newstd = np.sqrt(variance1d(xray_image, max=maxr, step=1,
                                position_matrix=position_matrix * iMask))
    plt.plot(t, oldstd, 'b', t, newstd, 'r', t, oldstd - newstd, 'g')
    plt.show()


def fit2d_save(Mask, path, imagename):
    """
    Compresses and wraps the mask for Fit2D use
    """
    Mask = np.flipud(Mask)
    os.chdir(path)
    currentdir = os.getcwd()
    Maskname = currentdir + "/" + imagename + ".msk"
    b1 = Mask.shape[0]
    b2 = Mask.shape[1]
    b28 = b2 / 8
    C = np.zeros((b1, b28), dtype=np.uint8)
    print 'start compression'
    # bitmap compression for Fit2D io
    # TODO: add to pyFAI via github fork/pull request maybe?
    for i in range(0, b1):
        j = 0
        while j < b2:
            C[i, math.floor(j / 8)] = Mask[i, j]
            for k in range(1, 8):
                C[i, math.floor(j / 8)] += Mask[i, j + k] * (2 ** k)
            j += 8
    fout = open(Maskname, "wb")
    C.tofile(fout, "", "%x")
    print('mask written')
    fout.close()
    # Giant long mask header, must be there for fit2d to work with the mask
    header = 'M\x00\x00\x00A\x00\x00\x00S\x00\x00\x00K\x00\x00\x00\x00\x08\x00\x00\x00\x08\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'
    mask = open(Maskname, 'rb')
    bodypt1 = mask.read(262144)
    bodypt2 = mask.read()
    mask.close()
    # fin.close()
    body = [bodypt1, bodypt2]
    body2 = "".join(body)
    total = [header, body2]
    total2 = "".join(total)
    fout = open(Maskname, "wb")
    fout.write(total2)
    fout.close()
    junk, Maskname2 = os.path.split(Maskname)
    Maskname2, junk = os.path.splitext(Maskname2)
    print('mask finalized')


# print 'call srxplanar'
# subprocess.call(['srxplanar',
# str(imagelocation),
# '-c /media/christopher/My\ Passport/Mar2014/xplanar_config.cfg',
# '--addmask '+str(Maskname2),
# ])
# unpack=np.loadtxt(imagename+'_qspace.chi',skiprows=32)
# plt.plot(unpack[:,0],unpack[:,1:]),plt.show()
# plt.plot(unpack[:,0],unpack[:,2]/unpack[:,1]*100), plt.show()
#
# subprocess.call(['pdfgetx3',
# imagename+'_qspace.chi',
# '--config=/media/christopher/My\ Passport/Mar2014/NiPdpdf.cfg'
# ])


def automask(ponifilepath=None, imagepath=None, imagearray=None,
             initial_m=None, edge=None,
             Fit2Dout=False):
    old_settings = np.seterr(divide='ignore', invalid='ignore')
    m = initial_m if initial_m is not None else 1.3
    edge = edge if edge is not None else 5
    # load geometry
    a = loadgeometry(ponifilepath)
    # load image
    xray_image = fabio.open(imagepath).data
    # generate radial array based on image shape
    Radi = a.rArray(xray_image.shape)
    # create integer array of radi by division by the pixel resolution
    roundR = np.around(Radi / a.pixel1).astype(int)
    # define the maximum radius for creation of numpy arrays
    maxr = np.amax(roundR)
    # make maxr into an integer so it can be used as a counter
    maxr = int(np.ceil(maxr))
    # initialize the maximums for the pixel counters
    imax, jmax = np.shape(roundR)

    # Max number of pixels in a ring: 2pi*r; 2pi~7
    # ML stores the intensities of each pixel in each ring
    ML = np.zeros((maxr + 1, 1024 * 7))
    # counter for the size of the bins for the pixels
    listcounter = np.bincount(roundR.ravel())
    i = 0
    while True:
        maskout = maskw(m, edge, xray_image, imax, jmax, maxr, roundR,
                        listcounter)
        plotmask(maskwoutput=maskout, maxr=maxr, xray_image=xray_image,
                 position_matrix=roundR)
        if raw_input('good?(y/[n])') == 'y':
            break
        else:
            m = float(raw_input('m='))
    if Fit2Dout is True:
        path, imagename = os.path.split(imagepath)
        imagename, ext = os.path.splitext(imagename)
        fit2d_save(maskout[3], path, imagename)
    np.seterr(**old_settings)

    # automask('/media/Seagate_Backup_Plus_Drive/Mar2014/Ni/17_21_24_NiSTD_300K-00002.poni',
    # '/media/Seagate_Backup_Plus_Drive/Mar2014/7_7_7_NiPd_Ambient/NiPd_LongRun/7_7_7_NiPd_AsIs_300K-00020.tif',
    # initial_m=1)