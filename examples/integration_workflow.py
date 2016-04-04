from __future__ import division, print_function

import sys
from xray_vision import QtGui, QtCore
from vis1d import LineStackExplorer
import matplotlib.pyplot as plt
from databroker import db, get_events
from metadatastore.api import db_connect as mds_db_connect

from filestore.api import db_connect as fs_db_connect
from filestore.api import retrieve
from sidewinder_spec.utils.handlers import *
from xpd_workflow.mask_tools import *
from diffpy.pdfgetx import PDFGetter
from datamuxer import DataMuxer

fs_db_connect(
    **{'database': 'data-processing-dev', 'host': 'localhost', 'port': 27017})
mds_db_connect(
    **{'database': 'data-processing-dev', 'host': 'localhost', 'port': 27017})

plot = False

# Get headers of interest
hdrs = db(
    # is_calibration=False
    is_calibration=True
)
hdrs = [db['b4875060-5a23-45da-bb8c-ef2cded62aa3']]
for hdr in hdrs:
    dm = DataMuxer()
    print(hdr['start']['run_folder'])

    # Get calibrations
    if not hdr['start']['is_calibration']:
        cals = [db[u]['start']['poni'][0] for u in hdr['start']['calibration']]
        print(cals)
    else:
        cals = [p for p in hdr['start']['poni']]

    geos = [retrieve(p) for p in cals]
    cal_dists = np.asarray([g.dist for g in geos]) * 100  # convert to meters

    # Get starting masks
    # start_masks = [retrieve(p) for p in hdr['start']['mask']]

    # Give the datamuxer our data
    dm.append_events(get_events(hdr))
    dm.to_sparse_dataframe()
    binned = dm.bin_on('img', interpolation={'T': 'linear'})
    b = binned[['I0',
                # 'T',
                'detz', 'img']]
    Ts = []
    rs = []
    grs = []

    for i, (i0,
            # T,
            detz, img) in b.iterrows():
        if detz < 30.:

            # Find the correct calibration file
            cal_idx = np.argmin((detz - cal_dists) ** 2)
            geo = geos[cal_idx]
            img /= geo.polarization(img.shape, .95)

            # start_mask = start_masks[cal_idx]
            start_mask = np.ones(img.shape, dtype=int).ravel().astype(bool)
            r = geo.rArray(img.shape)
            q = geo.qArray(img.shape)

            fr = r.ravel()
            fq = q.ravel()
            fimg = img.ravel()

            bins = 2000

            # Pre masking data
            bs_kwargs = {'bins': bins,
                         'range': [0, fq.max()]}

            # Needs to work on a real 2D image
            msk0 = mask_edge(img.shape, 30)

            tmsk = msk0 * start_mask

            '''
            for i in [2]:
                print(i)
                rbmsk = ring_blur_mask(img, r, geo.pixel1, i, mask=tmsk)
                print('total masked pixels', np.product(img.shape) - tmsk.sum())
                print('new masked pixels', rbmsk.sum() - tmsk.sum())
                print('new masked pixels',
                      (rbmsk.sum() - tmsk.sum()) / tmsk.sum() * 100., '%')
                print('pixels masked',
                      (rbmsk.sum() - tmsk.sum()) / img.size * 100., '%')
                tmsk = tmsk * rbmsk
            # '''
            fmsk_img = fimg[tmsk]
            fmsk_q = fq[tmsk]

            # Post masking data
            bs_args = (fmsk_q, fmsk_img)
            median = sts.binned_statistic(*bs_args, statistic='median',
                                          **bs_kwargs)
            plt.plot(median[0])
            plt.show()
            AAA
            print('start pdf')
            z = PDFGetter()
            d1 = {'qmin': 1.5, 'qmax': 34, 'qmaxinst': 34, 'rpoly': .9,
                  'rmax': 40.,
                  'composition': 'Pr2NiO4', 'dataformat': 'Qnm',
                  }

            d3 = {'qmin': 1.5, 'qmax': 25, 'qmaxinst': 25, 'rpoly': .9,
                  'rmax': 40.,
                  'composition': 'Pr2NiO4', 'dataformat': 'Qnm',
                  }
            d2 = {'qmin': 1.5, 'qmax': 28, 'qmaxinst': 28, 'rpoly': .9,
                  'rmax': 40.,
                  'composition': 'Pr2NiO4', 'dataformat': 'Qnm',
                  }

            if not plot:
                x = median[1][:-1]
            r, gr = z(x, median[0], **d1)
            rs.append(r)
            grs.append(gr)
            Ts.append(float(T))
    # key_list = [float(k) for k in b.I0.values]
    data_list = zip(rs, grs)
    print(data_list)
    app = QtGui.QApplication(sys.argv)
    tt = LineStackExplorer(Ts, data_list)
    tt.show()
    sys.exit(app.exec_())
