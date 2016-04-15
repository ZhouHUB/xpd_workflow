'''
This code will take a bunch of folders containing a configuration file and
data tiffs and perform:
1) read in all the not raw and not dark images in the folder

'''
import os
import ConfigParser
from xpd_workflow.mask_tools import *
from pims.tiff_stack import TiffStack_tifffile as TiffStack
from skbeam.io.save_powder_output import save_output
import fabio
import pyFAI

plot = True
save = False
target_folders = ['/path/to/folder']

for target_folder in target_folders:
    # move to the target folder
    os.chdir(target_folder)
    total_files = os.listdir('.')
    config_files = [f for f in total_files if f.startswith('config') and f.endswith('.txt')]
    for config_file in config_files:
        config = ConfigParser.ConfigParser()
        config.read(config_file)
        file_stem = config.get('sample_info', 'file_stem')


        # get all the files in the folder which are appropriate
        files = [f for f in os.listdir('.') if
                 f.endswith('.tif') and 'raw' not in f and 'dark' not in f]
        # Load the images
        img_stack = [TiffStack(f) for f in files]
        imgs = [np.rot90(s[0], 3) for s in img_stack]

        # sum the images
        img = np.sum(imgs, 0)

        # get the pyFAI geometry


        poni_file = config.get('beamline_info', 'configuration_file')
        geo = pyFAI.load(poni_file)

        start_mask_file = config.get('beamline_info', 'mask_file')
        if os.path.exists(start_mask_file):
            start_mask = fabio.open(start_mask_file).data
            start_mask = np.flipud(start_mask)
        else:
            start_mask = np.zeros(img.shape)

        img /= geo.polarization(shape=img.shape, factor=.95)

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

        msk_median_out = np.nan_to_num(msk_median[0])
        save_output(msk_median[1][:-1] / 10., msk_median_out,
                    os.path.join(folder, 'median'), 'Q')
