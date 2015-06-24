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


def lamda_from_bragg(th, d, n):
    return 2 * d * np.sin(th / 2.) / n


def find_peaks(chi, sides=6):
    # Find all potential peaks
    preliminary_peaks = scipy.signal.argrelmax(chi, order=20)[0]

    # peaks must have at least sides pixels of data to work with
    preliminary_peaks2 = preliminary_peaks[
        np.where(preliminary_peaks < len(chi) - sides)]

    # make certain that a peak has a drop off which causes the peak height to
    # be more than twice the height at sides pixels away
    criteria = chi[preliminary_peaks2] >= 2 * chi[preliminary_peaks2 + sides]
    criteria *= chi[preliminary_peaks2] >= 2 * chi[preliminary_peaks2 - sides]

    peaks = preliminary_peaks[np.where(criteria)]

    left_idxs = peaks - sides
    right_idxs = peaks + sides
    peak_centers = peaks
    return left_idxs, right_idxs, peak_centers


def calibrate_energy(stack, calibs, relative_positions, calibration_file,
                     sides=6):
    stack_peaks = []
    # For every image in the routine find the peaks
    for s, calib, pos in zip(stack, calibs, relative_positions):
        if not isinstance(s, np.ndarray):
            img = np.rot90(s[0], 3)
        if isinstance(calib, str):
            geo = pyFAI.load(calib)
        elif isinstance(calib, dict):
            geo = pyFAI.geometry.Geometry(**calib)

        r = geo.rArray(img.shape) * pq.m
        res = geo.pixel1 * pq.m
        bins = np.arange(0. * pq.m, r.max(), res)
        chi = scipy.stats.binned_statistic(r.ravel(), img.ravel(),
                                           statistic='mean',
                                           bins=bins)[0]

        lidxs, ridxs, peak_centers = find_peaks(chi, sides=sides)

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
        stack_peaks.append(fitted_peak_centers)

        # print fitted_peak_centers
        # plt.plot(bins[:-1], chi)
        # plt.plot(bins[:-1][peak_centers], chi[peak_centers], 'ro')
    # plt.show()
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
            print np.arctan(((d1 - d0) / dD).simplified)
            D0_l.append(-d0 * dD / (d0 - d1))
            thetas.append(np.arctan(((d1 - d0) / dD).simplified))
        theta_l.append(thetas)
        ave_D0 = np.sum(D0_l) * pq.mm / len(D0_l)
        D0s[i] += rdi + ave_D0
        D0s[j] += rdj + ave_D0
        n_DOs[i] += 1
        n_DOs[j] += 1
        # ave_D0s = D0s / n_DOs
        # print ave_D0s
        th = np.asarray(theta_l[-1]) * pq.radians
        d = np.loadtxt(calibration_file) * pq.angstrom
        d = d[:len(th)]
        # ns = [3, 4, 8, 11, 12]
        ns = [8, 6, 12, 24, 8]
        n = np.asarray(ns)
        lam = lamda_from_bragg(th, d, 1)
        e = pq.constants.h * pq.units.speed_of_light / lam
        es.append(e.rescale(pq.keV))
        lams.append(lam)
        print len(e)
        print np.average(e.rescale(pq.keV))
        print np.average(lam)
    print np.average(lams), np.average(es)
        # need a way to correctly average the positions


if __name__ == '__main__':
    import os
    import matplotlib.pyplot as plt
    from pims.tiff_stack import TiffStack_tifffile as TiffStack
    from pims.image_sequence import ImageSequence
    from pims.tiff_stack import TiffSeries
    # plt.ion()
    dir_name = '/mnt/bulk-data/research_data/energy_calib2'
    num = range(2)
    f_names = [os.path.join(dir_name, 'Ni-STD_Calib_D1-000' + str(i).zfill(2))
               for i in num]
    calibs = [f + '.poni' for f in f_names]

    # print ImageSequence(['/mnt/work-data/dev/xpd_workflow/energy_calib2/Ni-STD_Calib_D1-00000.tif'])
    # stack = TiffSeries([f + '.tif' for f in f_names])
    # stack = ImageSequence([f + '.tif' for f in f_names])
    stack = [TiffStack(f + '.tif') for f in f_names]
    relative_positions = np.asarray([13.80 * i for i in num]) * pq.mm
    calibrate_energy(stack, calibs, relative_positions, os.path.join(dir_name, 'Ni03.cal'))