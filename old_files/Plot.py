"""
Copyright (c) 2014 Brookhaven National Laboratory All rights reserved.
Use is subject to license terms and conditions.

@author: Christopher J. Wright"""
__author__ = 'Christopher J. Wright'

import matplotlib.pyplot as plt
import numpy as np
import Calculate


def plotmask(maskwoutput=None, maxr=None, xray_image=None,
             position_matrix=None):
    """
    Plots the mask metrics allowing for evaluation of mask quality
    """

    # highper,lowper,totalper, mask, imask,StatsA=maskwoutput
    highper, lowper, totalper, mask, imask = maskwoutput
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
    mask = plt.imshow(mask)
    plt.show()
    imagemask = plt.imshow(xray_image * imask)
    plt.show()
    t = np.arange(0, maxr, 1)
    plt.plot(t, Calculate.mean1d(xray_image,
                                 position_matrix=position_matrix * imask,
                                 step=1), 'b',
             # t, Median1D(xray_image, max=maxr, step=1, position_matrix=position_matrix*imask), 'r',\
             # t,Median1D(xray_image, max=maxr, step=1, position_matrix=position_matrix*imask)-Mean1D(xray_image, max=maxr, step=1, position_matrix=position_matrix*imask)
    )
    plt.show()

    oldstd = np.sqrt(
        Calculate.variance1d(xray_image, position_matrix=position_matrix,
                             max=maxr, step=1))
    newstd = np.sqrt(Calculate.variance1d(xray_image,
                                          position_matrix=position_matrix * imask,
                                          max=maxr, step=1))
    plt.plot(t, oldstd, 'b', t, newstd, 'r', t, oldstd - newstd, 'g')
    plt.show()


def plotring(image, position_matrix, angles, r, max, step, m):
    im = image[position_matrix == r]
    imang = angles[position_matrix == r]
    imang, im = zip(*sorted(zip(imang, im)))
    mean = Calculate.statistics1d(image, position_matrix=position_matrix,
                                  max=max, step=step, statistic='mean')
    med = Calculate.statistics1d(image, position_matrix=position_matrix,
                                 max=max, step=step, statistic='median')
    std = Calculate.statistics1d(image, position_matrix=position_matrix,
                                 max=max, step=step, statistic=np.std)
    # mean = np.mean(im)
    # med = np.median(im)
    # std = np.std(im)
    print 'std: ', std
    print 'std/mean%: ', std / mean * 100
    print 'abs(mean-median)/mean%: ', np.abs(mean - med) / mean * 100
    thresh = m * std
    size = np.shape(im)[0]
    print 'size: ', size
    meana = np.ones(size) * mean
    meda = np.ones(size) * med
    upperT = np.ones(size) * (mean + thresh)
    lowerT = np.ones(size) * (mean - thresh)
    t = np.arange(0, 2 * np.pi, 2 * np.pi / size)
    plt.ioff()
    plt.plot(t, im, 'r', t, meda, 'b', t, meana, 'g', t, upperT, 'o', t,
             lowerT, 'o')
    plt.show()