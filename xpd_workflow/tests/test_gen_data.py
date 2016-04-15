__author__ = 'christopher'
import fabio
import pyFAI
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pims.tiff_stack import TiffStack_tifffile as TiffStack
from skxray.io.save_powder_output import save_output
from xpd_workflow.mask_tools import *

geo = pyFAI.load(
    '/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/Ni_STD/Ni_PDF_60s-00000.poni')
dq = geo.deltaQ((2048, 2048))
q = geo.qArray((2048, 2048))
bins = 8000
# plt.imshow(dq)
# plt.show()
# AAA
# dq_mean = sts.binned_statistic(q.ravel(), dq.ravel(), bins=bins,
#                                           range=[0, q.max()], statistic='mean')
# dq_median = sts.binned_statistic(q.ravel(), dq.ravel(), bins=bi
#                                           range=[0, q.max()], statistic='median')
# plt.plot(dq_mean[1][:-1], dq_mean[0])
# plt.plot(dq_median[1][:-1], dq_median[0])
# plt.show()
r = geo.qArray((2048, 2048))
nr = r / np.max(r)

img = np.sin(nr * np.pi * 3) * np.exp(-10 * nr)
ideal_img = dc(img)
smax = np.max(img)
smin = np.min(img)
bad_pixels = []
'''
for i in xrange(np.random.randint(1000, 2000)):
    x, y = np.random.randint(0, 2048), np.random.randint(0, 2048)
    if np.random.random() >= .5:
        img[x, y] = smax * 3
    else:
        img[x, y] = smin * 3
    bad_pixels.append([x, y])
'''
plt.imshow(img, vmin=smin, vmax=smax)
plt.show()

# plt.imshow(idsr - dsr)
# plt.show()
# ideal_median = sts.binned_statistic(q.ravel(), ideal_img.ravel(), bins=bins,
#                                     range=[0, q.max()], statistic='median')
#
# ideal_mean = sts.binned_statistic(q.ravel(), ideal_img.ravel(), bins=bins,
#                                   range=[0, q.max()], statistic='mean')
# ideal_std = sts.binned_statistic(q.ravel(), ideal_img.ravel(), bins=bins,
#                                  range=[0, q.max()], statistic=np.std)


# median = sts.binned_statistic(q.ravel(), img.ravel(), bins=bins,
#                                     range=[0, q.max()], statistic='median')
#
# mean = sts.binned_statistic(q.ravel(), img.ravel(), bins=bins,
#                                   range=[0, q.max()], statistic='mean')
# std = sts.binned_statistic(q.ravel(), img.ravel(), bins=bins,
#                                  range=[0, q.max()], statistic=np.std)

# plt.plot(ideal_mean[1][:-1], ideal_mean[0], label='ideal mean')
# plt.plot(ideal_median[1][:-1], ideal_median[0], label='ideal median')
# plt.plot(ideal_std[1][:-1], ideal_std[0], label='ideal std')
# plt.legend()
# plt.show()

# plt.plot(mean[1][:-1], mean[0], label='mean')
# plt.plot(median[1][:-1], median[0], label='median')
# # plt.plot(std[1][:-1], std[0], label='ideal std')
# plt.legend()
# plt.show()

perfect_mask = (img - ideal_img) != 0
for i in [10,
          # 9, 8, 7, 6, 5, 4.5, 4
          ]:
    rbmsk = ring_blur_mask(img, geo, i)
    print i
    print 'good mask', np.sum(perfect_mask == rbmsk)
    print 'under masked', np.sum(perfect_mask > rbmsk)
    print 'over masked', np.sum(perfect_mask < rbmsk)
    print
# '''
plt.imshow(img, interpolation='none', origin='lower', aspect='auto')
for y, x in bad_pixels:
    plt.plot(x, y, 'ro', mfc='r', mec='r', ms=10)
for y, x in zip(
    np.where(rbmsk != 0)[0],
    np.where(rbmsk != 0)[1]
):
    plt.plot(x, y, 'go', mfc='g', mec='g', ms=5)
plt.show()
# '''
print q[1907, 173], q[173, 1907]

_, hist_bins, _ = plt.hist(img[np.where((q > 313.) & (q < 314.))], bins=50)
plt.axvline(np.mean(img[np.where((q > 313.) & (q < 314.))]), color='r')
plt.axvline(np.mean(img[np.where((q > 313.) & (q < 314.))]) + np.std(img[np.where((q > 313.) & (q < 314.))]))
plt.axvline(np.mean(img[np.where((q > 313.) & (q < 314.))]) -  np.std(img[np.where((q > 313.) & (q < 314.))]))
# plt.hist(img[np.where((q > 287.) & (q < 288.) & (rbmsk != 1))],
         # bins=50
         # bins=hist_bins
         # )
plt.show()
'''
mr = dc(q)
mr[rbmsk.astype(bool)] = -1

msk_median = sts.binned_statistic(mr.ravel(), img.ravel(), bins=bins,
                                  range=[0, mr.max()], statistic='median')
msk_mean = sts.binned_statistic(mr.ravel(), img.ravel(), bins=bins,
                                range=[0, mr.max()], statistic='mean')
msk_std = sts.binned_statistic(mr.ravel(), img.ravel(), bins=bins,
                               range=[0, mr.max()], statistic=np.std)

plt.plot(msk_mean[1][:-1], msk_mean[0], label='mean')
plt.plot(msk_median[1][:-1], msk_median[0], label='median')
# plt.plot(std[1][:-1], std[0], label='ideal std')
plt.legend()
plt.show()
# '''
