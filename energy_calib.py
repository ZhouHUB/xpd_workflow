__author__ = 'christopher'
import itertools
import pyFAI
import fabio
import lmfit
import scipy.stats
import numpy as np
import scipy.signal


def find_peaks(chi):
    # Find all potential peaks
    preliminary_peaks = scipy.signal.argrelmax(chi, order=20)[0]
    # peaks must have at least 6 pixels of data to work with
    preliminary_peaks2 = preliminary_peaks[
        np.where(preliminary_peaks < len(chi) - 6)]
    # make certain that a peak has a drop off which causes the peak height to
    # be more than twice the height at 6 pixels away
    criteria = chi[preliminary_peaks2] >= 2 * chi[preliminary_peaks2 + 6]
    criteria *=chi[preliminary_peaks2] >= 2 * chi[preliminary_peaks2 - 6]
    peaks = preliminary_peaks[np.where(criteria)]
    left_idxs = peaks - 6
    right_idxs = peaks + 6
    peak_centers = peaks
    return left_idxs, right_idxs, peak_centers


def calibrate_energy(imgs, calibs, relative_positions):
    all_x = []
    all_peaks = []
    chis = []
    for img, calib, pos in zip(imgs, calibs, relative_positions):
        geo = pyFAI.load(calib)
        r = geo.rArray(img.shape)
        res = geo.pixel1
        bins = np.arange(0, r.max(), res)
        chi = scipy.stats.binned_statistic(r.ravel(), img.ravel(),
                                           statistic='mean',
                                           bins=bins)[0]
        chis.append(chi)
        lidxs, ridxs, peak_centers = find_peaks(chi)
        fitted_peak_centers = []
        for lidx, ridx, peak_center in zip(lidxs, ridxs,
                                                    peak_centers):
            mod = lmfit.models.VoigtModel()
            pars = mod.guess(chi[lidx: ridx],
                             x=bins[lidx: ridx])
            out = mod.fit(chi[lidx: ridx], pars,
                          x=bins[lidx: ridx])
            center = out.values['center']
            # get peak center from out
            fitted_peak_centers.append(center)
        all_peaks.append(fitted_peak_centers)
        all_x.append(bins[:-1])
    # pair two detector positions
    D0s = np.zeros(len(all_x))
    n_DOs = np.zeros(len(all_x))
    for i,  peaksi, rdi, j, peaksj, rdj in itertools.combinations(
            zip(range(0, len(all_x)), all_peaks, relative_positions),2):
        # blob0 is the closer image, blob1 the farther
        # for each pair obtain the peak movement and the peak position
        dD = rdj-rdi
        # make certain that we only use the peaks that are common to both,
        # some of the peaks in the farther detector position could have fallen
        # off the map.
        peaksi = peaksi[:len(peaksj)]
        D0_l = []
        for d0, d1 in zip(peaksi, peaksj):
            D0_l.append(-d0 * dD / (d0 - d1))
        ave_D0 = sum(D0_l) / len(D0_l)
        D0s[i] += rdi + ave_D0
        D0s[j] += rdj + ave_D0
        n_DOs[i] += 1
        n_DOs[j] += 1
    ave_D0s = D0s/n_DOs
    return ave_D0s
    # need a way to correctly average the positions

if __name__ == '__main__':
    import os
    import matplotlib.pyplot as plt
    plt.ion()

    os.chdir('/mnt/bulk-data/research_data/enery_calib')
    calibs = ['Ni_STD_XRD1200-00000.poni', 'Ni_STD_XRD1500-00000.poni']

    f_names = ['Ni_STD_XRD1200-00000.tif', 'Ni_STD_XRD1500-00000.tif']
    imgs = [fabio.open(f).data for f in f_names]

    relative_positions = [0., 300.]

    all_x = []
    all_peaks = []
    chis = []
    for img, calib, pos in zip(imgs, calibs, relative_positions):
        geo = pyFAI.load(calib)
        r = geo.rArray(img.shape)
        res = geo.pixel1
        bins = np.arange(0, r.max(), res)
        chi = scipy.stats.binned_statistic(r.ravel(), img.ravel(),
                                           statistic='mean',
                                           bins=bins)[0]
        chis.append(chi)
        lidxs, ridxs, peak_centers = find_peaks(chi)
        fitted_peak_centers = []
        for lidx, ridx, peak_center in zip(lidxs, ridxs,
                                                    peak_centers):
            mod = lmfit.models.VoigtModel()
            pars = mod.guess(chi[lidx: ridx],
                             x=bins[lidx: ridx])
            out = mod.fit(chi[lidx: ridx], pars,
                          x=bins[lidx: ridx])
            center = out.values['center']
            # get peak center from out
            fitted_peak_centers.append(center)
        all_peaks.append(fitted_peak_centers)
        all_x.append(bins[:-1])

        print fitted_peak_centers
        plt.plot(bins[:-1], chi)
        plt.plot(bins[:-1][peak_centers], chi[peak_centers], 'ro')
    # plt.show()
    # pair two detector positions
    D0s = np.zeros(len(all_x))
    n_DOs = np.zeros(len(all_x))
    for (i, peaksi, rdi), (j, peaksj, rdj) in itertools.combinations(
            zip(range(0, len(all_x)), all_peaks, relative_positions), 2):
        # blob0 is the closer image, blob1 the farther
        # for each pair obtain the peak movement and the peak position
        dD = rdj-rdi
        # make certain that we only use the peaks that are common to both,
        # some of the peaks in the farther detector position could have fallen
        # off the map.
        peaksi = peaksi[:len(peaksj)]
        D0_l = []
        for d0, d1 in zip(peaksi, peaksj):
            D0_l.append(-d0 * dD / (d0 - d1))
        ave_D0 = sum(D0_l) / len(D0_l)
        D0s[i] += rdi + ave_D0
        D0s[j] += rdj + ave_D0
        n_DOs[i] += 1
        n_DOs[j] += 1
        print rdi, rdj, ave_D0
    ave_D0s = D0s/n_DOs
    print ave_D0s