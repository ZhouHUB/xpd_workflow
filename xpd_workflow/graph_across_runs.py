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

logger = logging.getLogger(__name__)

if __name__ == '__main__':
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from xpd_workflow.parsers import parse_xrd_standard

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
    line_style = ['--', '-.', ':',]
    save = True
    print('save', save)
    offset = .002
    output = 'chi'
    for ns in [
        [1, 2, 3, 4, 5],
        # [1, 16, 18, 20, 22]
    ]:
        ns.sort()
        for event_idx, event_name in zip([1, 54 // 2, -3],
                                         ['initial', 'high-temp', 'final']):
            legended_hkl = []
            for xmin, xmax, length in zip([
                # 0, 1000, 3000,
                0],
                    [
                        # 1000, 3000, 4000,
                        4000],
                    [
                        # 'short', 'medium', 'long',
                        'full']):
                fig = plt.figure(figsize=(26, 12))
                ax1 = fig.add_subplot(111)
                for i, j in enumerate(ns):
                    print(j)
                    folder = '/mnt/bulk-data/research_data/USC_beamtime/APS_March_2016/S' + str(
                        j) + '/temp_exp'
                    # key_list1 = [f for f in os.listdir(folder) if
                    #              f.endswith('.gr')]
                    key_list1 = [f for f in os.listdir(folder) if
                                 f.endswith(output) and f.startswith('0')]
                    key_list1.sort()

                    if event_idx == 54 / 2 and j == 1:
                        key = key_list1[83 // 2]
                    else:
                        key = key_list1[event_idx]
                    data = (
                        np.loadtxt(os.path.join(folder, key), skiprows=8)[:,
                        0],
                        np.loadtxt(os.path.join(folder, key), skiprows=8)[:,
                        1])

                    cm = plt.get_cmap('viridis')
                    cNorm = colors.Normalize(vmin=0, vmax=len(ns))
                    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
                    x, y = data
                    colorVal = scalarMap.to_rgba(i)
                    x /= 10.
                    ax1.plot(x[xmin:xmax], y[xmin:xmax] + i * offset,
                             # color=colorVal,
                             label='S{}'.format(j), lw=2)
                    ax1.set_xlim(1.5, 4.5)
                    if output == 'gr':
                        bnds = ['O-Pr', 'O-Ni', 'Ni-Ni', 'Pr-Pr', 'Ni-Pr',
                                'O-Pr',
                                'O-Ni',
                                'Ni-Ni-Ni', 'Pr-Ni', 'Pr-Pr', 'Pr-Ni-O',
                                'Ni-Pr-Ni',
                                'Pr-Pr']
                        bnd_lens = [2.320, 1.955, 3.883, 3.765, 3.186, 2.771,
                                    2.231,
                                    7.767, 4.426, 6.649, 4.989, 5.404, 3.374]
                        if xmin == 0:
                            for bnd, bnd_len in zip(bnds, bnd_lens):
                                ax1.axvline(bnd_len, color='grey',
                                            linestyle='--')
                            ax3 = ax1.twiny()
                            ax3.set_xticks(np.asarray(bnd_lens) / x[xmax])
                            ax3.set_xticklabels(bnds, rotation=90)
                    elif output == 'chi':
                        std_axis = []
                        for n, hkls, tths, color, ls in zip(standard_names,
                                                        master_hkl,
                                                        master_tth,
                                                        color_map, line_style):
                            std_axis.append(ax1.twiny())
                            ax3 = std_axis[-1]
                            hkl_q = np.pi * 4 * np.sin(
                                np.deg2rad(tths / 2)) / lam
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
                            ax3.set_xticks(
                                list((hkl_q[c] - ax1.get_xlim()[0]) / (
                                    ax1.get_xlim()[1] - ax1.get_xlim()[0])))
                            ax3.set_xticklabels(hkls, rotation=90, color=color)
                if output == 'chi':
                    fig_name = 'I(Q)'
                else:
                    fig_name = 'PDF'
                fig.suptitle(
                    'S{}-{} {} {}, {}'.format(min(ns), max(ns), fig_name, length,
                                               event_name))
                if output == 'gr':
                    ax1.set_xlabel(r"$r (\AA)$")
                    ax1.set_ylabel(r"$G (\AA^{-2})$")
                else:
                    ax1.set_xlabel(r"$Q (\AA^{-1})$")
                    ax1.set_ylabel(r"$I (Q)$")
                ax1.legend(loc='best', fontsize=20, framealpha=.1, fancybox=True)
                fig.tight_layout(rect=[0, 0, 1, .95], w_pad=1e-6)
                if save:
                    fig.savefig(os.path.join(
                        '/mnt/bulk-data/Dropbox/',
                        'S{}-{}_{}_{}_{}.png'.format(min(ns), max(ns), length,
                                                     event_name, output)))
                    fig.savefig(os.path.join(
                        '/mnt/bulk-data/Dropbox/',
                        'S{}-{}_{}_{}_{}.eps'.format(min(ns), max(ns), length,
                                                     event_name, output)))
                else:
                    plt.show()
