from __future__ import (division, print_function)

import matplotlib.cm as cmx
import matplotlib.colors as colors
from matplotlib import gridspec
from metadatastore.api import db_connect as mds_db_connect

from filestore.api import db_connect as fs_db_connect

fs_db_connect(
    **{'database': 'data-processing-dev', 'host': 'localhost', 'port': 27017})
mds_db_connect(
    **{'database': 'data-processing-dev', 'host': 'localhost', 'port': 27017})

from databroker import db, get_events
from datamuxer import DataMuxer
from sidewinder_spec.utils.handlers import *
import logging
from xpd_workflow.parsers import parse_xrd_standard

logger = logging.getLogger(__name__)

if __name__ == '__main__':
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    save = True

    lam = 1.54059
    # Standard reflections for sample components
    niox_hkl = ['111', '200', '220', '311', '222', '400', '331',
                '420', '422', '511']
    niox_tth = np.asarray(
        [37.44, 43.47, 63.20, 75.37, 79.87, 95.58, 106.72, 111.84,
         129.98, 148.68])
    pr3_hkl = ['100', '001', '110', '101', '111', '200', '002', '210', '211',
               '112', '202']
    pr3_tth = np.asarray(
        [22.96, 24.33, 32.70, 33.70, 41.18, 46.92, 49.86, 52.86, 59.00, 60.91,
         70.87]
    )
    pr4_hkl = ['111', '113', '008', '117', '200', '119', '028', '0014', '220',
               '131', '1115', '0214', '317', '31Na', '2214', '040', '400']
    pr4_tth = np.asarray(
        [23.43, 25.16, 25.86, 32.62, 33.36, 37.67, 42.19, 46.11, 47.44, 53.18,
         55.55, 57.72, 59.10, 59.27, 68.25, 68.71, 70.00]
    )
    pr2_tth, pr2int, pr2_hkl = parse_xrd_standard(
        '/mnt/bulk-data/research_data/Pr2NiO4orthorhombicPDF#97-008-1577.txt')
    pr2_tth = pr2_tth[pr2int > 5.]

    prox_hkl = ['111', '200', '220', '311', '222', '400', '331', '420', '422',
                '511', '440', '531', '600']
    prox_tth = np.asarray(
        [28.25, 32.74, 46.99, 55.71, 58.43, 68.59, 75.73, 78.08, 87.27,
         94.12, 105.63, 112.90, 115.42]
    )

    standard_names = [
        # 'NiO',
        'Pr3Ni2O7',
        'Pr2NiO4',
        # 'Pr4'
        'Pr6O11'
    ]
    master_hkl = [
        # niox_hkl,
        pr3_hkl,
        pr2_hkl,
        # pr4_hkl
        prox_hkl
    ]
    master_tth = [
        # niox_tth,
        pr3_tth,
        pr2_tth,
        # pr4_tth
        prox_tth
    ]
    color_map = [
        # 'red',
        'blue',
        'black',
        'red'
    ]
    line_style = ['--', '-.', ':', ]
    ns = [1, 2, 3, 4, 5,
          # 18, 20, 22, 16, 28, 29, 27, 26
          ]
    # ns = [26]
    ns.sort()
    #
    for i in ns:
        legended_hkl = []
        print(i)
        folder = '/mnt/bulk-data/research_data/USC_beamtime/APS_March_2016/S' + str(
            i) + '/temp_exp'
        hdr = db(run_folder=folder)[0]
        dm = DataMuxer()
        dm.append_events(get_events(hdr))
        df = dm.to_sparse_dataframe()
        print(df.keys())
        binned = dm.bin_on('img', interpolation={'T': 'linear'})

        # key_list = [f for f in os.listdir(folder) if
        #             f.endswith('.gr') and not f.startswith('d')]
        key_list = [f for f in os.listdir(folder) if
                    f.endswith('.chi') and not f.startswith('d') and f.strip(
                        '0.chi') != '' and int(
                        f.lstrip('0').strip('.chi')) % 2 == 1]
        key_list.sort()
        key_list = key_list[:-1]
        # key_list2.sort()
        idxs = [int(os.path.splitext(f)[0]) for f in key_list]
        Ts = binned['T'].values[idxs]
        output = os.path.splitext(key_list[0])[-1][1:]
        if key_list[0].endswith('.gr'):
            offset = .1
            skr = 0
        else:
            skr = 8
        offset = .001
        data_list = [(np.loadtxt(os.path.join(folder, f),
                                 skiprows=skr
                                 )[:, 0],
                      np.loadtxt(os.path.join(folder, f),
                                 skiprows=skr
                                 )[:, 1])
                     for f
                     in key_list]
        ylim_min = None
        for xmax, length in zip(
                [len(data_list[0][0]) - 1, len(data_list[0][0]) - 1],
                ['short', 'full']):
            fig = plt.figure(figsize=(26, 12))
            gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1])
            ax1 = plt.subplot(gs[0])
            if length == 'short':
                ax1.set_xlim(1.5, 4.5)
            ax2 = plt.subplot(gs[1], sharey=ax1)
            plt.setp(ax2.get_yticklabels(), visible=False)
            cm = plt.get_cmap('viridis')
            cNorm = colors.Normalize(vmin=0, vmax=len(key_list))
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
            for idx in range(len(key_list)):
                xnm, y = data_list[idx]
                colorVal = scalarMap.to_rgba(idx)
                if output == 'chi':
                    x = xnm / 10.
                ax1.plot(x[:xmax], y[:xmax] + idx * offset,
                         color=colorVal)
                ax2.plot(Ts[idx], y[-1] + idx * offset, marker='o',
                         color=colorVal)
                if ylim_min is None or ylim_min > np.min(
                        y[:xmax + idx * offset]):
                    ylim_min = np.min(y[:xmax + idx * offset])
            ax2.set_xticklabels([str(f) for f in ax2.get_xticks()],
                                rotation=90)
            if output == 'gr':

                bnds = ['O-Pr', 'O-Ni', 'Ni-Ni', 'Pr-Pr', 'Ni-Pr', 'O-Pr',
                        'O-Ni',
                        'Ni-Ni-Ni', 'Pr-Ni', 'Pr-Pr', 'Pr-Ni-O', 'Ni-Pr-Ni',
                        'Pr-Pr', 'Rs:Pr-Pr', 'Rs:Pr_Pr']
                bnd_lens = [2.320, 1.955, 3.883, 3.765, 3.186, 2.771, 2.231,
                            7.767, 4.426, 6.649, 4.989, 5.404, 3.374, 3.910,
                            8.801]
                # ax1.grid(True)
                # ax2.grid(True)
                for bnd, bnd_len in zip(bnds, bnd_lens):
                    ax1.axvline(bnd_len, color='grey', linestyle='--')
                ax3 = ax1.twiny()
                ax3.set_xticks(np.asarray(bnd_lens) / x[xmax])
                ax3.set_xticklabels(bnds, rotation=90)
            else:
                std_axis = []
                for n, hkls, tths, color, ls in zip(standard_names, master_hkl,
                                                master_tth,
                                                color_map, line_style):
                    std_axis.append(ax1.twiny())
                    ax3 = std_axis[-1]
                    hkl_q = np.pi * 4 * np.sin(np.deg2rad(tths / 2)) / lam
                    for k, (hkl, q) in enumerate(zip(hkls, hkl_q)):
                        if n not in legended_hkl:
                            ax1.axvline(q, color=color, linestyle=ls,
                                        lw=2,
                                        label=n
                                        )
                            legended_hkl.append(n)
                        else:
                            ax1.axvline(q, color=color, linestyle=ls,
                                        lw=2,
                                        )
                    a = hkl_q > ax1.get_xlim()[0]
                    b = hkl_q < ax1.get_xlim()[1]
                    c = a & b
                    ax3.set_xticks(list((hkl_q[c] - ax1.get_xlim()[0]) / (
                        ax1.get_xlim()[1] - ax1.get_xlim()[0])

                                        ))
                    ax3.set_xticklabels(hkls, rotation=90, color=color)
            ax2.set_xlabel('Temperature C')
            if output == 'gr':
                fig.suptitle('S{} PDF'.format(i))
                ax1.set_xlabel(r"$r (\AA)$")
                ax1.set_ylabel(r"$G (\AA^{-2})$")
            elif output == 'chi':
                fig.suptitle('S{} I(Q)'.format(i))
                ax1.set_xlabel(r"$Q (\AA^{-1})$")
                ax1.set_ylabel(r"$I (Q) $")
            ax1.set_ylim(ylim_min)
            ax1.legend()
            gs.tight_layout(fig, rect=[0, 0, 1, .98], w_pad=1e-6)

            if save:

                fig.savefig(os.path.join('/mnt/bulk-data/Dropbox/',
                                         'S{}_{}_output_{}.png'.format(
                                             i, length, output)))
                fig.savefig(os.path.join('/mnt/bulk-data/Dropbox/',
                                         'S{}_{}_output_{}.eps'.format(
                                             i, length, output)))
            else:
                plt.show()
