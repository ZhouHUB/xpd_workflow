import os

from databroker import db, get_events
from filestore.api import retrieve
from xpd_workflow.mask_tools import *

# Get headers of interest
hdrs = []
for hdr in hdrs:
    time_dept_bg = True
    geo1 = retrieve(hdr['calibration_file1'])
    geo2 = retrieve(hdr['calibration_file2'])
    start_mask1 = retrieve(hdr['start_mask1'])
    start_mask2 = retrieve(hdr['start_mask2'])
    for event in get_events(hdr):
        img = event['data']['img']
        r = geo.rArray(img.shape)
        q = geo.qArray(img.shape)
        fimg = img.flatten()

        msk0 = mask_beamstop(geo, 1.)
        msk1 = mask_edge(img, 10)
        msk2 = mask_radial_edge(img, geo, 31)
        initial_mask = msk0 | msk1 | msk2 | start_mask
        tmsk = msk0 | msk1 | msk2 | start_mask

        for i in [
            # 10,
            9, 8, 7, 6,
            # 5, 4.5, 4
        ]:
            print i
            rbmsk = ring_blur_mask(img, geo, i, mask=tmsk)
            print 'total masked pixels', tmsk.sum()
            print 'new masked pixels', rbmsk.sum() - tmsk.sum()
            print 'new masked pixels', (
                                           rbmsk.sum() - tmsk.sum()) / tmsk.sum() * 100, '%'
            print 'pixels masked', (
                                   rbmsk.sum() - tmsk.sum()) / img.size * 100, '%'
            tmsk = tmsk | rbmsk

        tmsk = tmsk.astype(np.bool)
        bins = 2000

        # Pre masking data
        median = sts.binned_statistic(q.ravel(), img.ravel(), bins=bins,
                                  range=[0, q.max()], statistic='median')
        mean = sts.binned_statistic(q.ravel(), img.ravel(), bins=bins,
                                    range=[0, q.max()], statistic='mean')
        std = sts.binned_statistic(q.ravel(), img.ravel(), bins=bins,
                                   range=[0, q.max()], statistic=np.std)

        # Post masking data
