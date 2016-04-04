from __future__ import division, print_function

import matplotlib.pyplot as plt
from databroker import db, get_events
from metadatastore.api import db_connect as mds_db_connect

from filestore.api import db_connect as fs_db_connect
from filestore.api import retrieve
from sidewinder_spec.utils.handlers import *
from xpd_workflow.mask_tools import *
from diffpy.pdfgetx import PDFGetter

# from matplotlib.colors import LogNorm

fs_db_connect(
    **{'database': 'data-processing-dev', 'host': 'localhost', 'port': 27017})
mds_db_connect(
    **{'database': 'data-processing-dev', 'host': 'localhost', 'port': 27017})

plot = False
# Get headers of interest
hdrs = db()
# hdrs = db()
for hdr in hdrs:
    print(hdr['start']['run_folder'])
    # Get calibrations
    if not hdr['start']['is_calibration']:
        cals = [db[u]['start']['poni'][0] for u in hdr['start']['calibration']]
        print(cals)
    else:
        cals = [p for p in hdr['start']['poni']]
    geos = [retrieve(p) for p in cals]
    cal_dists = np.asarray(
        [g.dist for g in geos]) * 100  # pyFAI reports in meters
    # Get starting masks
    # start_masks = [retrieve(p) for p in hdr['start']['mask']]
    events = [e for e in get_events(hdr) if 'detz' in e['data'].keys() and e['data']['detz'] < 30.]
    # for event in get_events(hdr):
    grs = []
    rs = []
    for event in events:
        # Pull relevant data into local vars
        data = event['data']
        img = data['img']
        detz = data['detz']

        # Find the correct calibration file, it's the one with the dist close
        # to the recorded detector dist
        cal_idx = np.argmin((detz - cal_dists) ** 2)
        geo = geos[cal_idx]
        # img /= geo.polarization(img.shape, .71)
        img /= geo.polarization(img.shape, .95)

        # start_mask = start_masks[cal_idx]
        start_mask = np.zeros(img.shape, dtype=int).ravel()
        r = geo.rArray(img.shape)
        q = geo.qArray(img.shape)

        fr = r.ravel()
        fq = q.ravel()
        fimg = img.ravel()

        # a, b, c = geo.integrate2d(img, 2000, 360)
        # a, b, c = geo.integrate2d(geo.polarization(img.shape, 1.5), 2000)
        # plt.imshow(a[150:, 1000:], aspect='auto')
        # plt.show()
        # aa = np.median(a, axis=0)

        # plt.imshow((a[150:, 1000:] - aa[1000:])**2, aspect='auto', norm=LogNorm(), cmap='viridis')
        # plt.show()
        # key_list = range(1000, 1100)
        # print(key_list)
        # data_list = [(np.arange(0, 360), a[:, i]) for i in key_list]

        # app = QtGui.QApplication(sys.argv)
        # tt = LineStackExplorer(key_list, data_list)
        # tt.show()
        # sys.exit(app.exec_())

        # x = np.arange(0, 360)
        # y = a[:, 1000]
        # plt.plot(x, y)
        # plt.show()

        # spl = UnivariateSpline(x, y)
        # plt.plot(x, y)
        # xs = np.linspace(0, 360, 36*2)
        # plt.plot(xs, spl(xs))
        # plt.plot()
        # plt.show()

        bins = 2000
        # Pre masking data
        bs_args = (fq, fimg)
        bs_kwargs = {'bins': bins,
                     'range': [0, fq.max()]}

        if plot:
            f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

            median = sts.binned_statistic(*bs_args, statistic='median',
                                          **bs_kwargs)
            mean = sts.binned_statistic(*bs_args, statistic='mean',
                                        **bs_kwargs)
            std = sts.binned_statistic(*bs_args, statistic=np.std, **bs_kwargs)

            x = median[1][:-1]
            ax1.plot(x, median[0], label='median no mask')
            ax2.plot(x, mean[0], label='mean no mask')
            ax3.plot(x, std[0], label='std no mask')
        # plt.legend(loc='best')
        # plt.show()

        # Needs to work on a real 2D image
        msk0 = mask_edge(img, 30)
        # msk0 = mask_beamstop(r, .01)
        # msk2 = mask_radial_edge(img, r, .2)
        initial_mask = msk0 | start_mask
        tmsk = msk0 | start_mask
        '''
        for i in [
            # 10,
            # 9, 8, 7, 6,
            # 5, 4.5,
            # 6
            2
        ]:
            print(i)
            rbmsk = ring_blur_mask(img, geo, i, mask=tmsk)
            print('total masked pixels', tmsk.sum())
            print('new masked pixels', rbmsk.sum() - tmsk.sum())
            print('new masked pixels',
                  (rbmsk.sum() - tmsk.sum()) / tmsk.sum() * 100., '%')
            print('pixels masked',
                  (rbmsk.sum() - tmsk.sum()) / img.size * 100., '%')
            tmsk = tmsk | rbmsk
        # '''
        tmsk = tmsk.astype(np.bool)
        # plt.imshow(tmsk)
        # plt.show()
        ftmsk = tmsk.ravel()
        fmsk_img = fimg[np.invert(ftmsk)]
        fmsk_q = fq[np.invert(ftmsk)]

        # Post masking data
        bs_args = (fmsk_q, fmsk_img)
        median = sts.binned_statistic(*bs_args, statistic='median',
                                      **bs_kwargs)
        if plot:
            mean = sts.binned_statistic(*bs_args, statistic='mean',
                                        **bs_kwargs)
            std = sts.binned_statistic(*bs_args, statistic=np.std, **bs_kwargs)
            ax1.plot(x, median[0], label='median')
            ax2.plot(x, mean[0], label='mean')
            ax3.plot(x, std[0], label='std')
            ax1.set_title('median')
            ax2.set_title('mean')
            ax3.set_title('std')
            # plt.legend(loc='best')
            plt.show()
        print('start pdf')
        z = PDFGetter()
        d1 = {'qmin': 1.5, 'qmax': 34, 'qmaxinst': 34, 'rpoly': .9, 'rmax':40.,
              'composition': 'Pr2NiO4', 'dataformat': 'Qnm',
              }

        d3 = {'qmin': 1.5, 'qmax': 25, 'qmaxinst': 25, 'rpoly': .9, 'rmax':40.,
              'composition': 'Pr2NiO4', 'dataformat': 'Qnm',
              }
        d2 = {'qmin': 1.5, 'qmax': 28, 'qmaxinst': 28, 'rpoly': .9, 'rmax': 40.,
              'composition': 'Pr2NiO4', 'dataformat': 'Qnm',
              }

        if not plot:
            x = median[1][:-1]
        r, gr = z(x, median[0], **d1)
        rs.append(r)
        grs.append(gr)
        plt.plot(r, gr)
        plt.show()
