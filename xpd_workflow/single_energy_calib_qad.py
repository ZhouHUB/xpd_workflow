__author__ = 'christopher'
import itertools
import fabio
import lmfit
import scipy.stats
import numpy as np
import scipy.signal
import quantities as pq
import matplotlib
import matplotlib.pyplot as plt

plt.ioff()
font = {'family': 'normal',
        # 'weight' : 'bold',
        'size': 18}

matplotlib.rc('font', **font)


def lamda_from_bragg(th, d, n):
    return 2 * d * np.sin(th / 2.) / n


def find_peaks(intensity, sides=6):
    # Find all potential peaks
    preliminary_peaks = scipy.signal.argrelmax(intensity, order=20)[0]

    # peaks must have at least sides pixels of data to work with
    preliminary_peaks2 = preliminary_peaks[
        np.where(preliminary_peaks < len(intensity) - sides)]

    # make certain that a peak has a drop off which causes the peak height to
    # be more than twice the height at sides pixels away
    criteria = intensity[preliminary_peaks2] >= 2 * intensity[
        preliminary_peaks2 + sides]
    criteria *= intensity[preliminary_peaks2] >= 2 * intensity[
        preliminary_peaks2 - sides]

    peaks = preliminary_peaks[np.where(criteria)]

    left_idxs = peaks - sides
    right_idxs = peaks + sides
    peak_centers = peaks
    return left_idxs, right_idxs, peak_centers


def calibrate_energy(intensity, ttheta, calibration_file, sides=12):
    # Find all the peaks in the diffraction pattern
    lidxs, ridxs, peak_centers = find_peaks(intensity, sides=sides)

    fitted_peak_centers = []
    for lidx, ridx, peak_center in zip(lidxs, ridxs, peak_centers):
        mod = lmfit.models.VoigtModel()
        pars = mod.guess(intensity[lidx: ridx], x=ttheta[lidx: ridx])
        out = mod.fit(intensity[lidx: ridx], pars, x=ttheta[lidx: ridx])
        center = pq.UncertainQuantity(out.values['center'], pq.m,
                                      out.params['center'].stderr)

        # get peak center from out
        fitted_peak_centers.append(center)

    plt.plot(ttheta[:-1], intensity)
    plt.plot(ttheta[:-1][peak_centers], intensity[peak_centers], 'ro')
    plt.show()

    # This allows the quick and dirty part of the code to work, we assume that
    # there are only two peaks in the pattern
    assert len(fitted_peak_centers) == 2
    es = []
    lams = []
    d_spacings = np.loadtxt(calibration_file) * pq.angstrom
    # Note here we make the assumption that our peaks start at the first
    # d spacing
    d = d_spacings[0]
    half_delta_theta = np.sum(np.abs(fitted_peak_centers)) / 2.
    lam = lamda_from_bragg(half_delta_theta, d, 1)
    e = pq.constants.h * pq.units.speed_of_light / lam
    es.append(e.rescale(pq.keV))
    lams.append(lam)
    nplams = None
    npe = None
    for l, e_ele in zip(lams, es):
        if nplams is None:
            nplams = l
            npe = e_ele
        else:
            nplams = np.concatenate((nplams, l))
            npe = np.concatenate((npe, e_ele))
    print np.mean(nplams * pq.angstrom), np.std(nplams)
    print np.mean(npe * pq.keV), np.std(npe)
    # need a way to correctly average the positions


if __name__ == '__main__':
    calibrant = 'path/to/file.txt'
    # load file here
    a = np.loadtxt('file.chi')
    intensity = a[:, 1]
    ttheta = a[:, 0]
    calibrate_energy(intensity, ttheta, calibrant)
