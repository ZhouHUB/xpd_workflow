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
                                                   cmap='viridis',
                                                   # norm='log',
                                                   intensity_scaling='percentile',
                                                   img_min=5.,
                                                   img_max=95.)

        self._main_window.setFocus()
        self.setCentralWidget(self._main_window)


if __name__ == '__main__':
    import os
    from pims.tiff_stack import TiffStack_tifffile as TiffStack
    import fabio

    folder = '/mnt/bulk-data/research_data/USC_beamtime/APS_March_2016/more_cells'
    # file_names = [f for f in os.listdir(folder) if
    #               f.endswith('.tif') and 'raw' not in f and 'dark' not in f
    #               and 'd95' not in f]
    file_names = ['S30_d95-00013.tif']
    print(file_names)
    files = [os.path.join(folder, f) for f in file_names]
    files.sort()
    files = files[:1]
    img_stack = [TiffStack(f) for f in files]

    imgs = [np.rot90(s[0], 3) for s in img_stack]
    # a = imgs[0] - imgs[1]
    # imgs = [s + np.abs(s.min()) + .1 for s in imgs]
    key_list = files
    data_list = imgs
    app = QtGui.QApplication(sys.argv)
    tt = StackExplorer(key_list, data_list)
    tt.show()
    sys.exit(app.exec_())
