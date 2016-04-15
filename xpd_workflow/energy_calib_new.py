__author__ = 'christopher'
import itertools
import pyFAI
import pyFAI.geometry
import fabio
import lmfit
import scipy.stats
import numpy as np
import scipy.signal
import quantities as pq
import matplotlib
import matplotlib.pyplot as plt
from mask_tools import *

plt.ioff()
font = {'family': 'normal',
        # 'weight' : 'bold',
        'size': 18}

matplotlib.rc('font', **font)


def lamda_from_bragg(th, d, n):
    return 2 * d * np.sin(th / 2.) / n


def find_peaks(chi, sides=6, intensity_threshold=0):
    # Find all potential peaks
    preliminary_peaks = scipy.signal.argrelmax(chi, order=20)[0]

    # peaks must have at least sides pixels of data to work with
    preliminary_peaks2 = preliminary_peaks[
        np.where(preliminary_peaks < len(chi) - sides)]

    # make certain that a peak has a drop off which causes the peak height to
    # be more than twice the height at sides pixels away
    criteria = chi[preliminary_peaks2] >= 2 * chi[preliminary_peaks2 + sides]
    criteria *= chi[preliminary_peaks2] >= 2 * chi[preliminary_peaks2 - sides]
    criteria *= chi[preliminary_peaks2] >= intensity_threshold

    peaks = preliminary_peaks[np.where(criteria)]

    left_idxs = peaks - sides
    right_idxs = peaks + sides
    peak_centers = peaks
    return left_idxs, right_idxs, peak_centers


def calibrate_energy(stack, calibs, relative_positions, calibration_file,
                     sides=12, intensity_threshold=6000):
    stack_peaks = []
    # For every image in the routine find the peaks
    for s, calib, pos in zip(stack, calibs, relative_positions):
        if not isinstance(s, np.ndarray):
            # img = np.rot90(s[0], 3)
            img = s[0]
        if isinstance(calib, str):
            geo = pyFAI.load(calib)
        elif isinstance(calib, dict):
            geo = pyFAI.geometry.Geometry(**calib)

        r = geo.rArray(img.shape) * pq.m
        # r[r < .01] = -1. * pq.m
        res = geo.pixel1 * pq.m
        bins = np.arange(0. * pq.m, r.max(), res)
        mask = np.ones(img.shape, dtype=bool)
        mask *= mask_edge(img.shape, 30)
        mask *= ring_blur_mask(img, r.magnitude, geo.pixel1, 2,
                               mask=mask)
        bs_kwargs = {'statistic': (np.median,),
                     'bins': bins,
                     'range': [0, r.max()]}

        a, x = binstats(r[mask], img[mask], **bs_kwargs)
        chi = np.asarray(a[0])
        lidxs, ridxs, peak_centers = find_peaks(
            chi, sides=sides, intensity_threshold=intensity_threshold)

        fitted_peak_centers = []
        for lidx, ridx, peak_center in zip(lidxs, ridxs,
                                           peak_centers):
            mod = lmfit.models.VoigtModel()
            pars = mod.guess(chi[lidx: ridx],
                             x=bins[lidx: ridx])
            out = mod.fit(chi[lidx: ridx], pars,
                          x=bins[lidx: ridx])
            center = pq.UncertainQuantity(out.values['center'], pq.m,
                                          out.params['center'].stderr)
            # center = out.values['center'] * pq.m
            # get peak center from out
            fitted_peak_centers.append(center)
        print(fitted_peak_centers)
        stack_peaks.append(fitted_peak_centers)

        # print fitted_peak_centers
        plt.plot(x, chi)
        plt.plot(x[peak_centers], chi[peak_centers], 'ro')
    plt.show()
    # pair two detector positions
    D0s = np.zeros(len(stack_peaks)) * pq.mm
    n_DOs = np.zeros(len(stack_peaks))
    theta_l = []
    es = []
    lams = []
    for (i, peaksi, rdi), (j, peaksj, rdj) in itertools.combinations(
            zip(range(0, len(stack_peaks)), stack_peaks, relative_positions),
            2):
        # blob0 is the closer image, blob1 the farther
        # for each pair obtain the peak movement and the peak position
        dD = rdj - rdi
        # make certain that we only use the peaks that are common to both,
        # some of the peaks in the farther detector position could have fallen
        # off the map.
        peaksi = peaksi[:len(peaksj)]
        D0_l = []
        thetas = []

        for d0, d1 in zip(peaksi, peaksj):
            # print np.arctan(((d1 - d0) / dD).simplified)
            D0_l.append(-d0 * dD / (d0 - d1))
            thetas.append(np.arctan(((d1 - d0) / dD).simplified))
        theta_l.append(thetas)
        ave_D0 = np.sum(D0_l) * pq.mm / len(D0_l)
        D0s[i] += rdi + ave_D0
        D0s[j] += rdj + ave_D0
        n_DOs[i] += 1
        n_DOs[j] += 1
        th = np.asarray(theta_l[-1]) * pq.radians
        d = np.loadtxt(calibration_file) * pq.angstrom
        d = d[:len(th)]
        lam = lamda_from_bragg(th, d, 1)
        e = pq.constants.h * pq.units.speed_of_light / lam
        es.append(e.rescale(pq.keV))
        lams.append(lam)
        # print len(e)
        # print np.average(e.rescale(pq.keV))
        # print np.average(lam)
    nplams = None
    npe = None
    for l, e_ele in zip(lams, es):
        if nplams is None:
            nplams = l
            npe = e_ele
        else:
            nplams = np.concatenate((nplams, l))
            npe = np.concatenate((npe, e_ele))
    print(nplams, npe)
    print np.mean(nplams), np.std(nplams)
    print np.mean(npe), np.std(npe)
    # need a way to correctly average the positions


if __name__ == '__main__':
    import os
    import subprocess

    from pims.tiff_stack import TiffStack_tifffile as TiffStack
    from pims.image_sequence import ImageSequence
    from pims.tiff_stack import TiffSeries

    calibrant = os.path.join('/mnt/bulk-data/research_data/energy_calib2',
                             'Ni03.cal')
    # dir_name = '/mnt/bulk-data/research_data/energy_calib2'
    dir_name = '/mnt/bulk-data/Dropbox/Temp-Data (2)'
    num = [0, 1, 3, 4]
    f_stem = 'Ni_STD_2rev_asc-000'
    f_names = [os.path.join(dir_name, f_stem + str(i).zfill(2))
               for i in num]
    calibs = [f + '.poni' for f in f_names]

    print f_names
    stack = [TiffStack(f + '.tif') for f in f_names]
    relative_positions = np.asarray([27.6 * 1. * i for i in num]) * pq.mm
    calibrate_energy(stack, calibs, relative_positions,
                     calibrant)
