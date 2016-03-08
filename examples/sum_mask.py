__author__ = 'christopher'
import fabio
import pyFAI
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pims.tiff_stack import TiffStack_tifffile as TiffStack
from skxray.io.save_powder_output import save_output
from xpd_workflow.mask_tools import *

# load geometry
geo = pyFAI.load(
    '/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/Ni_STD/Ni_PDF_60s-00000.poni'
)
# load starting mask, mostly of the beam pole
start_mask = fabio.open(
    '/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/xpd_pdf_beamstop.msk').data

# load images
folder, name = ('/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/Au2nm_RT', '2nm_Au_RT_median_sum')
files = [os.path.join(folder, f) for f in os.listdir(folder) if
         f.endswith('.tif') and 'raw' not in f and 'dark' not in f]
img_stack = [TiffStack(f) for f in files]
imgs = [np.rot90(s[0], 3) for s in img_stack]
img = np.sum(imgs, 0)

# correct image
img /= geo.polarization(shape=img.shape, factor=.95)

# produce masks
start_mask = np.flipud(start_mask)
msk0 = mask_beamstop(img, geo, .005)
msk1 = mask_edge(img, 10)
msk2 = mask_radial_edge(img, geo, 310)

initial_mask = msk0 | msk1 | msk2 | start_mask
tmsk = msk0 | msk1 | msk2 | start_mask

for i in [10,
          9, 8, 7, 6,
          5, 4.5, 4
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

# msk_save_name = os.path.join(folder, 'xpd_pdf_mask')
# np.save(msk_save_name, tmsk)
# fit2d_save(tmsk, msk_save_name + '.msk')

# msk_median_out = np.nan_to_num(msk_median[0])
# save_output(msk_median[1][:-1] / 10., msk_median_out,
#             os.path.join(folder, name), 'Q')
