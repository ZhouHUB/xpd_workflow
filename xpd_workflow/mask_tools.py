from __future__ import division
__author__ = 'christopher'

import numpy as np
import scipy.stats as sts
from copy import deepcopy as dc
import math
import os
import matplotlib.pyplot as plt


def mask_beamstop(pixel_positions, beamstop_radius, mask=None):
    if mask is None:
        mask = np.zeros(pixel_positions.shape)
    z = np.where(pixel_positions < beamstop_radius)
    mask[z] = 1
    return mask.astype(int)


def mask_radial_edge(img, pixel_positions, inner_radius, mask=None):
    if mask is None:
        mask = np.zeros(img.shape)
    z = np.where(pixel_positions > inner_radius)
    mask[z] = 1
    return mask.astype(int)


def mask_edge(img_shape, edge_size):
    """
    Mask the edge of an image

    Parameters
    -----------
    img_shape: tuple
        The shape of the image
    edge_size: int
        Number of pixels to mask from the edge
    Returns
    --------
    1darray:
        The raveled mask array, bad pixels are 0
    """
    mask = np.zeros(img_shape)
    mask[:, :edge_size] = 1
    mask[:, -edge_size:] = 1
    mask[:edge_size, :] = 1
    mask[-edge_size:, :] = 1
    return mask.ravel().astype(bool)


def low_pixel_count_mask(img, geometry, bins, min_pixels, mask=None):
    if mask is None:
        mask = np.zeros(img.shape)
    r = geometry.rArray(img.shape)
    int_r = np.around(r / geometry.pixel1).astype(int)
    mr = dc(r.ravel())
    mr[mask] = -1

    pixels_per_bin = \
        sts.binned_statistic(mr, img.ravel(), bins=bins, range=[0, mr.max()],
                             statistic='count')[0]

    mask = mask | pixels_per_bin[int_r] < min_pixels
    return mask.astype(int)


def ring_blur_mask(fimg, fr, rsize, alpha, bins=None, mask=None):
    """
    Perform a annular mask, which checks the ring statistics and masks any
    pixels which have a value greater or less than alpha * std away from the
    mean
    Parameters
    ----------
    fimg: 1darray
        The flattened image
    fr: 1darray
        The flattened array which maps pixels to radii
    alpha: float or 1darray
        Then number of acceptable standard deviations
    bins: int, optional
        Number of bins used in the integration, if not given then max number of
        pixels +1
    mask: 1darray
        A starting flattened mask
    Returns
    --------
    1darray:
        The flattened mask
    """

    if mask is None:
        mask = np.zeros(img.shape).ravel().astype(bool)
    int_r = np.around(fr / rsize).astype(int)
    if bins is None:
        bins = int_r.max() + 1
    print(bins)
    fmsk_img = fimg[np.invert(mask)]
    fmsk_r = fr[np.invert(mask)]

    # integration
    mean = sts.binned_statistic(fmsk_r, fmsk_img, bins=bins,
                                range=[0, fr.max()], statistic='mean')[0]
    std = sts.binned_statistic(fmsk_r, fmsk_img, bins=bins,
                               range=[0, fr.max()], statistic=np.std)[0]
    print(mean)
    print(mean.shape)
    print(std.shape)
    plt.plot(mean)
    plt.plot(std)
    plt.show()
    threshold = alpha * std
    lower = mean - threshold
    upper = mean + threshold

    # single out the too low and too high pixels
    too_low = fmsk_img < lower[int_r]
    too_hi = fmsk_img > upper[int_r]

    mask = mask | too_low | too_hi
    mask = np.invert(mask.astype(bool))
    return mask.astype(bool)


def invert_mask(mask):
    return np.abs(mask - 1).astype(int)


def fit2d_save(Mask, image_file_name):
    """
    Compresses and wraps the mask for Fit2D use
    """
    Mask = np.flipud(Mask)
    Mask = np.abs(Mask - 1)
    print image_file_name
    path, imagenamext = os.path.split(image_file_name)
    image_name = os.path.splitext(imagenamext)[0]
    if path is not '':
        os.chdir(path)
    currentdir = os.getcwd()
    Maskname = os.path.join(currentdir, image_name + ".msk")
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
    print('mask compressed')
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
    print('mask finalized')
    return Maskname


def bin_rings(img, qa, bins, mask=None):
    if mask is not None:
        mask = mask
    else:
        mask = np.zeros(img.shape).astype(bool)

    fmask = mask.flatten()
    fimg = img.flatten()

    mapping = np.digitize(qa.flatten(), bins)
    master_list = [[] for i in range(len(bins))]

    for pixel, m, tf in zip(fimg, mapping - 1, fmask):
        if not tf:
            master_list[m].append(pixel)

    return master_list


if __name__ == '__main__':
    import fabio
    import pyFAI
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from pims.tiff_stack import TiffStack_tifffile as TiffStack


    # load experiment information
    geo = pyFAI.load(
        '/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/Ni_STD/Ni_PDF_60s-00000.poni'
        # '/mnt/bulk-data/research_data/Low T-MST-summed/Ni_STD_XRD300_60S-00000.poni'
        # '/mnt/bulk-data/Dropbox/BNL_Project/misc/CGO_summed/CGO_summed/Ni_stnd/Ni_STD_60s-00004.poni'
    )

    # start_mask = fabio.open(
        # '/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/xpd_pdf_beamstop.msk'
        # '/mnt/bulk-data/research_data/Low T-MST-summed/low_T_mst.msk'
    # ).data
    # start_mask = np.zeros((2048, 2048)).astype(bool)
    # start_mask = np.flipud(start_mask)
    start_mask = np.load(
        '/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/Au2nm_Cold/xpd_pdf_mask_Au2nm_60s_120K_sum.npy'
        )

    # load images
    # folder, name = ('/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/Au2nm_Temp/Au2nm_300-400', 'Au2nm_300-400')
    # folder = '/mnt/bulk-data/research_data/Low T-MST-summed/'
    # folder = '/mnt/bulk-data/Dropbox/BNL_Project/misc/CGO_summed/CGO_summed/Sample1_350um'
    folder = '/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/Au2nm_Cold/'
    # folder = '/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/SiGe/'

    files = [os.path.join(folder, f) for f in os.listdir(folder) if
             f.endswith('.tif') and 'raw' not in f and 'dark' not in f
             # and 'Sample-2' not in f
             ]
    sort_files = [f.rstrip('.tif') for f in os.listdir(folder) if
                  f.endswith('.tif') and 'raw' not in f and 'dark' not in f
                  # and 'Sample-2' not in f
                  ]
    img_stack = [TiffStack(f) for f in files]
    imgs = [np.rot90(s[0], 3) for s in img_stack]
    # Uncomment the following line if you want to sum all the images in the folder
    img = np.sum(imgs[:-1], 0)
    fname = 'Sample-2_sum'
    # correct image
    # img /= geo.polarization(shape=img.shape, factor=.95)
    # img *= geo.solidAngleArray(img.shape)

    # produce masks
    msk0 = mask_beamstop(geo, .005)
    msk1 = mask_edge(img, 10)
    msk2 = mask_radial_edge(img, geo, 310)

    initial_mask = msk0 | msk1 | msk2 | start_mask
    tmsk = msk0 | msk1 | msk2 | start_mask
    plt.imshow(invert_mask(tmsk) * (img - np.min(img) + .1),
               norm=LogNorm(),
               interpolation='none', aspect='auto'), plt.colorbar()
    plt.show()

    for i in [10,
              # 9, 8, 7, 6,
              # 5, 4.5, 4
              ]:
        print i
        rbmsk = ring_blur_mask(img, geo, i, mask=tmsk)
        print 'total masked pixels', tmsk.sum()
        print 'new masked pixels', rbmsk.sum() - tmsk.sum()
        print 'new masked pixels', (
                                       rbmsk.sum() - tmsk.sum()) / tmsk.sum() * 100, '%'
        print 'pixels masked', (rbmsk.sum() - tmsk.sum()) / img.size * 100, '%'
        tmsk = tmsk | rbmsk

    tmsk = tmsk.astype(np.bool)
    r = geo.qArray(img.shape)
    bins = 2000

    # integration
    median = sts.binned_statistic(r.ravel(), img.ravel(), bins=bins,
                                  range=[0, r.max()], statistic='median')
    mean = sts.binned_statistic(r.ravel(), img.ravel(), bins=bins,
                                range=[0, r.max()], statistic='mean')
    std = sts.binned_statistic(r.ravel(), img.ravel(), bins=bins,
                               range=[0, r.max()], statistic=np.std)

    mr = dc(r)
    mr[tmsk.astype(np.bool)] = -1

    msk_median = sts.binned_statistic(mr.ravel(), img.ravel(), bins=bins,
                                      range=[0, mr.max()], statistic='median')
    msk_mean = sts.binned_statistic(mr.ravel(), img.ravel(), bins=bins,
                                    range=[0, mr.max()], statistic='mean')
    msk_std = sts.binned_statistic(mr.ravel(), img.ravel(), bins=bins,
                                   range=[0, mr.max()], statistic=np.std)
    # Comment out the following line if you want to view the mask and
    # integrated plots, note that these are blocking so you will need to
    # close each window for the results to be saved.
    # '''
    plt.imshow(invert_mask(tmsk) * (img - np.min(img) + .1),
               norm=LogNorm(),
               interpolation='none', aspect='auto'), plt.colorbar()
    plt.show()

    plt.plot(msk_median[1][:-1], msk_median[0], label='mask')
    plt.plot(median[1][:-1], median[0], label='no mask')
    plt.plot(median[1][:-1], msk_median[0] - median[0], label='mask-no_mask')
    plt.title('median')
    plt.legend()
    plt.show()

    plt.plot(msk_mean[1][:-1], msk_mean[0], label='mask')
    plt.plot(mean[1][:-1], mean[0], label='no mask')
    plt.plot(mean[1][:-1], msk_mean[0] - mean[0], label='mask-no_mask')
    plt.title('mean')
    plt.legend()
    plt.show()

    plt.plot(msk_std[1][:-1], msk_std[0], label='mask')
    plt.plot(std[1][:-1], std[0], label='no mask')
    plt.plot(std[1][:-1], msk_std[0] - std[0], label='mask-no_mask')
    plt.title('std')
    plt.legend()
    plt.show()
    # '''
    '''
    msk_save_name = os.path.join(folder, 'xpd_pdf_mask_' + fname)
    np.save(msk_save_name, tmsk)
    fit2d_save(tmsk, msk_save_name + '.msk')

    msk_median_out = np.nan_to_num(msk_median[0])
    save_output(msk_median[1][:-1] / 10., msk_median_out,
                os.path.join(folder, fname + '_median'), 'Q')
    msk_mean_out = np.nan_to_num(msk_mean[0])
    save_output(msk_mean[1][:-1] / 10., msk_mean_out,
                os.path.join(folder, fname + '_mean'), 'Q')
    '''