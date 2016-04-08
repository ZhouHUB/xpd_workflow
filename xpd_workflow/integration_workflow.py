from __future__ import division, print_function

from databroker import db, get_events
from datamuxer import DataMuxer
from diffpy.pdfgetx import PDFGetter
from metadatastore.api import db_connect as mds_db_connect

from filestore.api import db_connect as fs_db_connect
from filestore.api import retrieve
from sidewinder_spec.utils.handlers import *
from xpd_workflow.mask_tools import *

fs_db_connect(
    **{'database': 'data-processing-dev', 'host': 'localhost', 'port': 27017})
mds_db_connect(
    **{'database': 'data-processing-dev', 'host': 'localhost', 'port': 27017})


def main():
    plot = False
    out = {}
    # Get headers of interest
    hdrs = db(
        is_calibration=False
        # is_calibration=True
    )
    bg_hdr = db(
        run_folder='/mnt/bulk-data/research_data/USC_beamtime/APS_March_2016/'
                   'Quartz_Background/temp_exp')
    bg_dm = DataMuxer()
    bg_dm.append_events(get_events(bg_hdr))
    bg_binned = bg_dm.bin_on('img', interpolation={'T': 'linear'})
    for hdr in hdrs:
        eout = {}
        dm = DataMuxer()
        print(hdr['start']['run_folder'], hdr['start']['uid'])

        # Get calibrations
        if not hdr['start']['is_calibration']:
            cals = [db[u]['start']['poni'][0] for u in
                    hdr['start']['calibration']]
            print(cals)
        else:
            cals = [p for p in hdr['start']['poni']]

        geos = [retrieve(p) for p in cals]
        cal_dists = np.asarray(
            [g.dist for g in geos]) * 100  # convert to meters

        # Get starting masks
        # start_masks = [retrieve(p) for p in hdr['start']['mask']]

        # Give the datamuxer our data
        dm.append_events(get_events(hdr))
        df = dm.to_sparse_dataframe()
        print(df.keys())
        binned = dm.bin_on('img', interpolation={'T': 'linear'})
        b = binned[['I0',
                    'T',
                    'detz', 'img']]
        Ts = []
        rs = []
        grs = []
        bg_idx = None
        bg_chi = None
        tmsk = None
        j = 0
        for i, (i0,
                T,
                detz, img) in b.iterrows():
            if detz < 30.:

                # Find the correct calibration file
                cal_idx = np.argmin((detz - cal_dists) ** 2)
                geo = geos[cal_idx]
                img /= geo.polarization(img.shape, .95)

                r = geo.rArray(img.shape)
                q = geo.qArray(img.shape)

                bins = 2000

                # Pre masking data
                bs_kwargs = {'bins': bins,
                             'range': [0, q.max()]}

                if tmsk is None:
                    # Needs to work on a real 2D image
                    # tmsk = start_masks[cal_idx]
                    tmsk = np.ones(img.shape, dtype=int).ravel().astype(bool)
                    tmsk *= mask_edge(img.shape, 30)
                    tmsk *= ring_blur_mask(img, r, geo.pixel1, (4.5, 2.),
                                           mask=tmsk)

                fmsk_img = img[tmsk]
                fmsk_q = q[tmsk]

                # Post masking data
                bs_args = (fmsk_q, fmsk_img)
                median = sts.binned_statistic(*bs_args, statistic='median',
                                              **bs_kwargs)
                x = median[1][:-1]

                # Get background singal at correct temperature (or closest to)
                temp_metric = np.abs(T - bg_binned['T'].values)

                # Throw out the ones not at the PDF distance
                temp_metric[bg_binned['detz'].values > 30.] = 1e6
                bg_idx2 = np.argmin(temp_metric)

                # Only reload a file if needed
                if bg_idx is None or bg_idx != bg_idx2 or bg_chi is None:
                    bg_idx = bg_idx2
                    bg_img = bg_binned.img.values[bg_idx][0]

                    bg_fmsk_q = q[tmsk]
                    bg_fmsk_img = bg_img.ravel()[tmsk]
                    bs_args = (bg_fmsk_q, bg_fmsk_img)
                    print(i0 / bg_binned.I0.values[bg_idx])
                    print('optimal temp', bg_binned['T'].values[bg_idx],
                          bg_binned['detz'].values[bg_idx])
                    bg_chi = sts.binned_statistic(*bs_args, statistic='median',
                                                  **bs_kwargs)[0]
                y = median[0]
                bg_y = i0 * bg_chi / bg_binned.I0.values[bg_idx] * 1.5

                print('start pdf', T)
                z = PDFGetter()
                d1 = {'qmin': 1.5, 'qmax': 29.5, 'qmaxinst': 29.5, 'rpoly': .9,
                      'rmax': 40.,
                      'composition': 'Pr2NiO4', 'dataformat': 'Qnm',
                      }

                r, gr = z(x, y - bg_y, **d1)
                if i < 2 and plot:
                    f, (ax1, ax2, ax3) = plt.subplots(3, 1)
                    ax1.plot(x, y, 'b-')
                    ax1.plot(x, bg_y, 'g-')
                    ax1.plot(x, y - bg_y, 'r--')
                    ax2.plot(z.fq[0], z.fq[1])
                    ax3.plot(r, gr)
                    plt.show()
                rgr = np.vstack((r, gr)).T
                qfq = np.vstack(z.fq).T
                # np.savetxt(os.path.join(hdr['start']['run_folder'],
                #                         str(i).zfill(5) + '.gr'), rgr)
                # np.savetxt(os.path.join(hdr['start']['run_folder'],
                #                         str(i).zfill(5) + '.fq'), qfq)
                j += 1
                rs.append(r)
                grs.append(gr)
                Ts.append(float(T))
                eout[i] = {'T': T, 'r': r, 'gr': gr, 'fq': z.fq[1], 'q': z.fq[0]}
        out[hdr['start']['uid']] = eout
        # data_list = zip(rs, grs)
        # print(data_list)
        # app = QtGui.QApplication(sys.argv)
        # tt = LineStackExplorer(Ts, data_list)
        # tt.show()
        # sys.exit(app.exec_())
    return out


if __name__ == '__main__':
    main()
