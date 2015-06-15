"""
Copyright (c) 2014 Brookhaven National Laboratory All rights reserved.
Use is subject to license terms and conditions.

@author: Christopher J. Wright"""
__author__ = 'Christopher J. Wright'


import scipy.stats
import scipy.optimize as op
import numpy as np
import scipy.ndimage.filters as snf
import matplotlib.pyplot as plt
import subprocess
import os


def corrections(xray_image, geometry_object, solid=False, pol=True,
                polarization_factor=0.0):
    """
    This function applies the geometry based corrections, both solid angle and
    polarization to the image

    Parameters
    -----------
    xray_image: image array
        This is the x-ray image to be corrected
    solid: bool
        Whether or not to use the solid angle correction
    pol: bool
        Wheter or not to use the polarization correction
    geometry_object: instance of pyFAI geometry
        object that holds the geometry as loaded form the poni file
    polarization_factor: float
        The polaraization to use in the polarization correction

    Returns
    -------
    array:
        The corrected image

    """
    # TODO: EXAMPLE CODE

    corrected_image = xray_image
    if solid is True:
        corrected_image /= geometry_object.solidAngleArray(
            xray_image.shape)
    if pol is True:
        corrected_image /= geometry_object.polarization(
            xray_image.shape,
            polarization_factor)
    return corrected_image


def optomize_polarization(r, xray_image, geometry_object, round_radi, angles):
    """
    This function gets the polarization factor for which the named ring has
    the smallest standard deviation.  This results in the ring being most
    flat uniform in intensity.  This is generally applied to the largest ring,
    allowing for the best statistical sampling.

    Parameters
    ----------
    r: float
        The number of the ring to optimize
    xray_image: 2darray
        The image as a numpy array
TODO: fix this to an actual object instance
    geometry_object: instance of pyFAI geometry
        object that holds the geometry as loaded form the poni file
    round_radi: array
        The radi/Q/two_theta array which maps detector x,y positions to
        reciprical space
    angles: array
        The angle array which maps detector x,y positions to azimuthal angles

    Returns
    -------
    float:
        The optimized polaraization factor
    """

    def distance(c):
        modimage = xray_image / geometry_object.polarization(
            xray_image.shape, c)
        Im = modimage[round_radi == r]
        ImAng = angles[round_radi == r]
        ImAng, Im = zip(*sorted(zip(ImAng, Im)))
        minthis = np.std(Im)
        # print maxthis
        return minthis

    res = op.minimize_scalar(distance, bounds=(0, 2.5), method='bounded')
    print 'Polarization factor should be ', res.x
    return res.x


def optimize_background_level(background, integrated, bg_refine_qmin_index=0,
                              bg_refine_qmax_index=None):
    """
    Optimize the background scale factor by matching the integrated
    background to the integrated frame

    Parameters
    ----------
    :param bg_refine_qmin_index:
    background: ndarray
        Integrated background
    integrated: ndarray
        Integrated frame

    Returns
    -------
    float:
        Background scale factor

    """
    print 'start optimization'
    background[np.isnan(background)] = 0
    integrated[np.isnan(integrated)] = 0
    background = background.astype(np.float64)
    integrated = integrated.astype(np.float64)
    if bg_refine_qmax_index is None:
        bg_refine_qmax_index = integrated.shape[0]
    background = background[bg_refine_qmin_index:bg_refine_qmax_index]
    integrated = integrated[bg_refine_qmin_index:bg_refine_qmax_index]

    def minthis(z):
        corrected_background = background * z
        old_settings = np.seterr(divide='ignore')
        minthis_function = (1 / integrated) * (
                                                  integrated - corrected_background) ** 2
        minthis_function[np.isnan(minthis_function)] = 0
        mean_minthis_function = np.mean(minthis_function)
        np.seterr(**old_settings)
        return mean_minthis_function

    res = op.minimize_scalar(minthis, bounds=(0, 6), method='bounded')
    return res.x


def mean1d(img, position_matrix, step):
    """
    This function bins the rings and generates their mean.

    Parameters
    ----------
    img: np.ndarray
        The image as a numpy array
    max: int, I hope
        The maximum number of the ring measurement, radi,Q, two_theta
    step: int, I hope
        The resolution of the pixels in radi Q, two_theta
    position_matrix: np.ndarray
        Map between reciprical and x,y space

    Returns
    -------
    np.ndarray:
        The mean of each ring
    """
    max = np.max(position_matrix)
    bins = np.arange(0, max + step, step)
    pixels_per_bin = np.array(np.histogram(position_matrix, bins)[0],
                              dtype=float)
    return np.histogram(position_matrix, bins, weights=img)[0] / pixels_per_bin


def statistics1d(img, position_matrix, max=None, step=None, statistic=None):
    """
    Generic function which takes in rings and puts out 1d array of scalars
    generated by the statistic function.

    Parameters
    ----------
    img: np.ndarray
        The image as a numpy array
    max: int, I hope
        The maximum number of the ring measurement, radi,Q, two_theta
    step: int, I hope
        The resolution of the pixels in radi Q, two_theta
    position_matrix: np.ndarray
        Map between reciprical and x,y space
    statistic: str or function
        If str, a string which references the statistical function to be used, see scipy documentation
        If function, a function which takes 1d arrays and generates scalars

    Returns
    -------
    np.ndarray:
        The statistic applied to each ring

    ..note:: the 2*step is because the last element descripes the rightmost
    edge, here we are interested in including the last pixel, which is at
    max+step so we need to add two steps
    """
    max = max if max is not None else np.max(position_matrix)
    bins = np.arange(0, max + 2 * step, step)
    flatposMatrix = np.ravel(position_matrix)
    flatimg = np.ravel(img)
    return \
        scipy.stats.binned_statistic(flatposMatrix, flatimg,
                                     statistic=statistic,
                                     bins=bins)[0]


def variance1d(pic, position_matrix, max=None, step=None):
    """
    From SrXplanar, the variance of the 2d detector

    Parameters
    ----------
    pic: np.ndarray
        The image as a numpy array
    max: int, I hope
        The maximum number of the ring measurement, radi,Q, two_theta
    step: int, I hope
        The resolution of the pixels in radi Q, two_theta
    position_matrix: p
        Map between reciprical and x,y space

    Returns
    -------
    np.ndarray:
        The variance of each ring
    """

    picavg = snf.uniform_filter(pic, 5, mode='wrap')
    pics2 = (pic - picavg) ** 2
    pvar = snf.uniform_filter(pics2, 5, mode='wrap')

    old_settings = np.seterr(divide='ignore')
    gain = np.divide(pvar, pic)

    np.seterr(**old_settings)
    gain[np.isnan(gain)] = 0
    gain[np.isinf(gain)] = 0
    gainmedian = np.median(gain)
    var = pic * gainmedian
    bins = np.arange(0, max + 2 * step, step)
    pixels_per_bin = np.array(np.histogram(position_matrix, bins)[0],
                              dtype=float)
    variance = np.histogram(position_matrix, bins, weights=var)[0]
    return variance / pixels_per_bin


def generate_mask(xray_image, position_matrix=None, max=None, average=None,
                  std=None, m=None, edge=10, pixels_per_bin=None):
    """
    This function generates a mask for the given image based upon the
    statistical variances in the rings

    Parameters
    ----------
    img: array
        The image as a numpy array
    integer_position_matrix: array of ints
        Map between reciprical and x,y space

    max: int, I hope
        The maximum number of the ring measurement, radi,Q, two_theta
    average: array
        The ring averages
    std: array
        The ring standard deviations
    m: float
        User specified number which is multiplied by the standard deviation to
        generate the statistical threshold.
        Larger m means less pixels masked, in general
    step: int, I hope
        The resolution of the pixels in radi Q, two_theta



    statistic: str or function
        If str, a string which references the statistical function to be used,
        see scipy documentation
        If function, a function which takes 1d arrays and generates scalars

    Returns
    -------
    highper, lowper: float
        Percentage of pixels which were too high/low
    totalper: float
        Percentage of pixels masked
    Mask: 2d array
        The mask for the image, in this case 1 is maksed and 0 is not
    iMask: 2d array
        The maks inverted, 0 is masked, 1 is not
    """

    threshold = m * std
    lower = average - threshold
    upper = average + threshold
    # start masking based on too high/too low/not enough pixels

    toolow = xray_image < lower[position_matrix]
    toohi = xray_image > upper[position_matrix]
    lowcounter = np.bincount(position_matrix[toolow],
                             minlength=max + 1)
    highcounter = np.bincount(position_matrix[toohi],
                              minlength=max + 1)
    # TODO: maybe read magic number for number of pixels from config file or
    # calculate based on detector info?
    # combine toolow, toohigh, and number of pixels too low
    Mask = toolow | toohi | (pixels_per_bin[position_matrix] < 10)
    # mask edges
    Mask[:, :edge] = 1
    Mask[:, -edge:] = 1
    Mask[:edge, :] = 1
    Mask[-edge:, :] = 1
    # percentages for graphing
    highper = (highcounter[:] / pixels_per_bin[:]) * 100
    lowper = (lowcounter[:] / pixels_per_bin[:]) * 100
    totalper = highper + lowper
    iMask = np.abs(Mask - 1)
    print 'finished generating mask'
    return highper, lowper, totalper, Mask, iMask


def write_mask(xray_image, geometry_object, mask=None, initial_m=None,
               resolution=0.02):
    """
    Wrapper that handles the pre and post processing for the mask_generator

    Parameters
    ----------
    xray_image: ndarray
        The image to be masked
    geometry_object: AzimuthalIntegrator object
        This contains the geometry of the experiment
    mask: str
        file name of additional mask to load
    initial_m: float
        Threshold scale factor
    resolution: float
        Q resolution in 1/A for the experiment

    Returns
    ------
    ndarray:
        Mask, 0 is masked
    """
    print'start mask write'
    a = geometry_object

    Radi = a.rArray(xray_image.shape)
    q_matrix = a.qArray(xray_image.shape) / 10

    # create integer array of radi by division by the pixel resolution
    roundR = np.around(Radi / a.pixel1).astype(int)

    # define the maximum radius for creation of numpy arrays
    # make maxr into an integer so it can be used as a counter
    maxr = int(np.ceil(np.amax(roundR)))

    # apply any optimization based corrections to produce a modulated image,
    # not used in the integration
    # only used for generating flat rings for masking
    modImage = corrections(xray_image, geometry_object=a, solid=False,
                           pol=True, polarization_factor=.95)

    # Define a bunch of statistics on the unmasked image
    rawmean = statistics1d(modImage, position_matrix=roundR, max=maxr, step=1,
                           statistic='mean')

    rawstd = statistics1d(modImage, position_matrix=roundR, max=maxr, step=1,
                          statistic=np.std)

    pixels_per_bin = statistics1d(modImage, position_matrix=roundR, max=maxr,
                                  step=1, statistic='count')

    plt.imshow(modImage)
    plt.show()
    m = initial_m if initial_m is not None else 1
    while True:
        highper, lowper, totalper, Mask, iMask = generate_mask(xray_image,
                                                               position_matrix=roundR,
                                                               max=maxr,
                                                               average=rawmean,
                                                               std=rawstd, m=m,
                                                               edge=5,
                                                               pixels_per_bin=pixels_per_bin)

        plt.imshow(iMask * modImage, interpolation='nearest')
        plt.colorbar()
        plt.show()

        plt.plot(np.arange(0, highper.shape[0], 1), highper[:], '-',
                 np.arange(0, lowper.shape[0], 1), lowper[:], 'o',
                 np.arange(0, totalper.shape[0], 1), totalper[:], '^')
        plt.show()

        newmean = statistics1d(modImage, position_matrix=q_matrix * iMask,
                               max=np.amax(q_matrix), step=resolution,
                               statistic='mean')
        q_rawmean = statistics1d(modImage, position_matrix=q_matrix,
                                 max=np.amax(q_matrix), step=resolution,
                                 statistic='mean')
        q = np.arange(0, np.max(q_matrix) + resolution, resolution)

        fig, ax = plt.subplots()
        ax.plot(q, q_rawmean, 'r', label='No mask')
        ax.plot(q, newmean, 'g', label='Mask')

        ax.legend(loc='upper right')
        plt.show()

        if raw_input('good?(y/[n])') == 'y':
            return iMask
        else:
            m = float(raw_input('m='))


def pdfgetx3(chi_file, background_file, background_level, composition,
             qmax, qmin, qmaxinst, pdf_format, output):
    if background_file is '':
        subprocess.call(['pdfgetx3',
                         '--format=' + str(pdf_format),
                         '--composition=' + str(composition),
                         '--qmaxinst=' + str(qmaxinst),
                         '--qmax=' + str(qmax),
                         '--qmin=' + str(qmin),
                         '--qmaxinst=' + str(qmaxinst),
                         '-o', str(output),
                         '-t', 'fq, gr',
                         '--force',
                         str(chi_file)])
    else:
        subprocess.call(['pdfgetx3',
                         '-b', str(background_file),
                         '--format=' + str(pdf_format),
                         '--composition=' + str(composition),
                         '--bgscale=' + str(background_level),
                         '--qmaxinst=' + str(qmaxinst),
                         '--qmax=' + str(qmax),
                         '--qmin=' + str(qmin),
                         '--qmaxinst=' + str(qmaxinst),
                         '-o', str(output),
                         '-t', 'fq, gr',
                         '--force',
                         str(chi_file)])
    return os.path.splitext(chi_file)[0] + '.gr', background_level, qmax