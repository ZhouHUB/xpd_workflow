from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys

from xray_vision import QtGui, QtCore
import numpy as np

from xray_vision.messenger.mpl.cross_section_2d import CrossSection2DMessenger

from xray_vision.qt_widgets import CrossSectionMainWindow

import logging

logger = logging.getLogger(__name__)


class StackExplorer(QtGui.QMainWindow):
    def __init__(self, key_list, data_list, parent=None):
        QtGui.QMainWindow.__init__(self, parent=parent)
        self.setWindowTitle('StackExplorer')
        self._main_window = CrossSectionMainWindow(data_list=data_list,
                                                   key_list=key_list,
                                                   cmap='cubehelix')

        self._main_window.setFocus()
        self.setCentralWidget(self._main_window)


if __name__ == '__main__':
    import os
    from pims.tiff_stack import TiffStack_tifffile as TiffStack
    import fabio

    # folder, name = ('/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/Au2nm_RT', '2nm_Au_RT_median_sum')
    # folder = '/mnt/bulk-data/research_data/Low T-MST-summed'
    # folder ='/mnt/bulk-data/Dropbox/BNL_Project/misc/CGO_summed/CGO_summed/Sample1_350um'
    # folder = '/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/Au2nm_Cold/'
    folder = '/mnt/bulk-data/research_data/USC_beamtime/08-05-2015/2015-08-05/SiGe/'
    file_names = [f for f in os.listdir(folder) if
                  f.endswith('.tif') and 'raw' not in f and 'dark' not in f
                  and 'Sample-1' not in f]
    print(file_names)
    files = [os.path.join(folder, f) for f in file_names]
    # mask_files = [os.path.join(folder, 'xpd_pdf_mask_'+f.split('.')[0] + '.npy') for f in
    #               file_names]
    # msks = [fabio.open(f).data for f in mask_files]
    # msks = [np.load(f) for f in mask_files]
    img_stack = [TiffStack(f) for f in files]

    imgs = [np.rot90(s[0], 3) for s in img_stack]
    a = imgs[0] - imgs[1]
    imgs = [s + np.abs(s.min()) + .1 for s in imgs]
    # imgs = [img * np.invert(msk) for (img, msk) in zip(imgs, msks)]
    # imgs = [s + .1 for s in imgs]
    key_list = files
    data_list = imgs
    key_list.append('difference')
    b = a - np.min(a) + .1
    data_list.append(b)
    app = QtGui.QApplication(sys.argv)
    tt = StackExplorer(key_list, data_list)
    tt.show()
    sys.exit(app.exec_())
