import os

from databroker import db, get_events
from filestore.api import retrieve
from xpd_workflow.mask_tools import *
from sidewinder_spec.utils.handlers import *
import matplotlib.pyplot as plt
from filestore.api import db_connect as fs_db_connect
from metadatastore.api import db_connect as mds_db_connect
fs_db_connect(**{'database': 'data-processing-dev', 'host': 'localhost', 'port': 27017})
mds_db_connect(**{'database': 'data-processing-dev', 'host': 'localhost', 'port': 27017})

# Get headers of interest
hdrs = [db[-1]]
for hdr in hdrs:
    time_dept_bg = True

    # Get calibrations
    geos = [retrieve(p) for p in hdr['start']['poni']]
    cal_dists = np.asarray(
        [g.dist for g in geos]) * 10  # pyFAI reports in meters
    # Get starting masks
    # start_masks = [retrieve(p) for p in hdr['start']['mask']]
    for event in get_events(hdr):
        # Pull relevent data into local vars
        data = event['data']
        img = data['img']
        detz = data['detz']

        # Find the correct calibration file, it's the one with the dist close
        # to the recorded detector dist
        cal_idx = np.argmin((detz - cal_dists) ** 2)
        geo = geos[cal_idx]
        # start_mask = start_masks[cal_idx]
        start_mask = np.zeros(img.shape, dtype=int)
        r = geo.rArray(img.shape)
        q = geo.qArray(img.shape)

        plt.imshow(img[np.invert])

        fr = r.ravel()
        fq = q.ravel()
        fimg = img.ravel()

        bins = 2000
        # Pre masking data
        bs_args = (fq, fimg)
        bs_kwargs = {'bins': bins,
                     'range': [0, fq.max()]}

        # median = sts.binned_statistic(*bs_args, statistic='median', **bs_kwargs)
        # mean = sts.binned_statistic(*bs_args, statistic='mean', **bs_kwargs)
        # std = sts.binned_statistic(*bs_args, statistic=np.std, **bs_kwargs)
        #
        # plt.plot(median[1][:-1], median[0], label='median no mask')
        # plt.plot(mean[1][:-1], mean[0], label='mean no mask')
        # plt.plot(std[1][:-1], std[0], label='std no mask')
        # plt.legend(loc='best')
        # plt.show()

        # Needs to work on a real 2D image
        msk1 = mask_edge(img, 10)
        msk0 = mask_beamstop(q, 10)
        msk2 = mask_radial_edge(img, q, 340)
        initial_mask = msk0 | msk1 | msk2 | start_mask
        tmsk = msk0 | msk1 | msk2 | start_mask
        plt.imshow(tmsk)
        plt.show()

        for i in [
            # 10,
            # 9, 8, 7, 6,
            # 5, 4.5,
            4
        ]:
            print(i)
            rbmsk = ring_blur_mask(img, geo, i, mask=tmsk)
            print('total masked pixels', tmsk.sum())
            print('new masked pixels', rbmsk.sum() - tmsk.sum())
            print('new masked pixels',
                  (rbmsk.sum() - tmsk.sum()) / tmsk.sum() * 100, '%')
            print('pixels masked',
                  (rbmsk.sum() - tmsk.sum()) / img.size * 100, '%')
            tmsk = tmsk | rbmsk

        tmsk = tmsk.astype(np.bool)
        plt.imshow(tmsk)
        plt.show()
        ftmsk = tmsk.ravel()
        fmsk_img = fimg[np.invert(ftmsk)]
        fmsk_q = fq[np.invert(ftmsk)]

        # Post masking data
        mr = dc(r)
        mr[tmsk.astype(np.bool)] = -1

        bs_args = (fmsk_q, fmsk_img)
        median = sts.binned_statistic(*bs_args, statistic='median',
                                      **bs_kwargs)
        mean = sts.binned_statistic(*bs_args, statistic='mean', **bs_kwargs)
        std = sts.binned_statistic(*bs_args, statistic=np.std, **bs_kwargs)
        plt.plot(median[1][:-1], median[0], label='median')
        plt.plot(mean[1][:-1], mean[0], label='mean ')
        plt.plot(std[1][:-1], std[0], label='std')
        plt.legend(loc='best')
        plt.show()
