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

    # folder = '/media/sf_Data/5/S5/temp_exp'
    # folder = '/media/sf_Data/18/S18/temp_exp'
    # folder = '/media/sf_Data/S6/'
    # ns = [1, 2, 3, 4, 5, 18, 20, 22, 16, 28, 29, 27, 26]
    ns = [6]
    ns.sort()
    for i in ns:
        print(i)
        folder = '/mnt/bulk-data/research_data/USC_beamtime/APS_March_2016/S' + str(
            i) + '/temp_exp'
        hdr = db(run_folder=folder)[0]
        dm = DataMuxer()
        dm.append_events(get_events(hdr))
        df = dm.to_sparse_dataframe()
        print(df.keys())
        binned = dm.bin_on('img', interpolation={'T': 'linear'})
        key_list1 = [f for f in os.listdir(folder) if f.endswith('.gr')]
        key_list2 = [f for f in os.listdir(folder) if f.endswith('.fq')]
        # key_list = [f for f in os.listdir(folder) if f.endswith('.chi')]
        key_list1.sort()
        key_list2.sort()
        idxs = [int(os.path.splitext(f)[0]) for f in key_list1]
        Ts = binned['T'].values[idxs]
        if key_list1[0].endswith('.gr') or key_list2[0].endswith('.fq'):
            skr = 27
        else:
            skr = 4
        data_list1 = [(np.loadtxt(os.path.join(folder, f),
                                  # skiprows=skr
                                  )[:, 0],
                       np.loadtxt(os.path.join(folder, f),
                                  # skiprows=skr
                                  )[:, 1])
                      for f
                      in key_list1]
        data_list2 = [(np.loadtxt(os.path.join(folder, f),
                                  # skiprows=skr
                                  )[:, 0],
                       np.loadtxt(os.path.join(folder, f),
                                  # skiprows=skr
                                  )[:, 1])
                      for f
                      in key_list2]

        for xmax, n in zip([1000, 4000], ['short', 'full']):
            fig = plt.figure(figsize=(26, 12))
            gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1])
            ax1 = plt.subplot(gs[0])
            ax2 = plt.subplot(gs[1], sharey=ax1)
            plt.setp(ax2.get_yticklabels(), visible=False)
            cm = plt.get_cmap('viridis')
            cNorm = colors.Normalize(vmin=0, vmax=len(key_list1))
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
            for idx in range(len(key_list1)):
                x, y = data_list1[idx]
                colorVal = scalarMap.to_rgba(idx)
                ax1.plot(x[:xmax], y[:xmax] + idx * .3, color=colorVal)
                ax2.plot(Ts[idx], y[-1] + idx * .3, marker='o', color=colorVal)
            ax2.set_xticklabels([str(f) for f in ax2.get_xticks()],
                                rotation=90)
            bnds = ['O-Pr', 'O-Ni', 'Ni-Ni', 'Pr-Pr', 'Ni-Pr', 'O-Pr', 'O-Ni',
                    'Ni-Ni-Ni', 'Pr-Ni', 'Pr-Pr', 'Pr-Ni-O', 'Ni-Pr-Ni',
                    'Pr-Pr']
            bnd_lens = [2.320, 1.955, 3.883, 3.765, 3.186, 2.771, 2.231,
                        7.767, 4.426, 6.649, 4.989, 5.404, 3.374]
            # ax1.grid(True)
            # ax2.grid(True)
            for bnd, bnd_len in zip(bnds, bnd_lens):
                ax1.axvline(bnd_len, color='grey', linestyle='--')
            ax3 = ax1.twiny()
            ax3.set_xticks(np.asarray(bnd_lens) / x[xmax])
            ax3.set_xticklabels(bnds, rotation=90)

            fig.suptitle('S{} PDF'.format(i))
            ax2.set_xlabel('Temperature C')
            ax1.set_xlabel(r"$r (\AA)$")
            ax1.set_ylabel(r"$G (\AA^{-2})$")
            gs.tight_layout(fig, rect=[0, 0, 1, 1], w_pad=1e-6)
            # fig.savefig(os.path.join(folder, '{}_output.png'.format(n)))
            # fig.savefig(os.path.join(folder, '{}_output.eps'.format(n)))
            plt.show()
