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

    save = False
    print('save', save)
    offset = .002
    output = 'chi'
    for ns in [
        [1, 2, 3, 4, 5],
        # [1, 16, 18, 20, 22]
    ]:
        ns.sort()
        for event_idx, event_name in zip([1, 54 // 2, -3],
                                         ['initial', 'high temp', 'final']):
            for xmin, xmax, n in zip([
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
                    np.loadtxt(os.path.join(folder, key), skiprows=8)[:, 0],
                    np.loadtxt(os.path.join(folder, key), skiprows=8)[:, 1])

                    cm = plt.get_cmap('viridis')
                    cNorm = colors.Normalize(vmin=0, vmax=len(ns))
                    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
                    x, y = data
                    colorVal = scalarMap.to_rgba(i)
                    ax1.plot(x[xmin:xmax], y[xmin:xmax] + i * offset,
                             # color=colorVal,
                             label='S{}'.format(j))
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

                fig.suptitle(
                    'S{}-{} PDF {}, {}'.format(min(ns), max(ns), n,
                                               event_name))
                if output == 'gr':
                    ax1.set_xlabel(r"$r (\AA)$")
                    ax1.set_ylabel(r"$G (\AA^{-2})$")
                else:
                    ax1.set_xlabel(r"$Q (\AA^{-1})$")
                    ax1.set_ylabel(r"$I (Q)$")
                ax1.legend(loc='best')
                fig.tight_layout(rect=[0, 0, 1, .95], w_pad=1e-6)
                if save:
                    fig.savefig(os.path.join(
                        '/mnt/bulk-data/Dropbox/',
                        'S{}-{}_{}_{}_{}.png'.format(min(ns), max(ns), n,
                                                     event_name, output)))
                    fig.savefig(os.path.join(
                        '/mnt/bulk-data/Dropbox/',
                        'S{}-{}_{}_{}_{}.eps'.format(min(ns), max(ns), n,
                                                     event_name, output)))
                else:
                    plt.show()
