from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys

from xray_vision import QtGui, QtCore
import numpy as np

from xray_vision.messenger.mpl.cross_section_2d import CrossSection2DMessenger

from xray_vision.qt_widgets import Stack1DMainWindow

import logging

logger = logging.getLogger(__name__)


class LineStackExplorer(QtGui.QMainWindow):
    def __init__(self, key_list, data_list, parent=None):
        QtGui.QMainWindow.__init__(self, parent=parent)
        self.setWindowTitle('StackExplorer')
        self._main_window = Stack1DMainWindow(data_list=data_list,
                                              key_list=key_list,
                                              cmap='viridis'
                                              )

        self._main_window.setFocus()
        self.setCentralWidget(self._main_window)


if __name__ == '__main__':
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    # folder = '/media/sf_Data/5/S5/temp_exp'
    # folder = '/media/sf_Data/18/S18/temp_exp'
    # folder = '/media/sf_Data/S6/'
    folder = '/mnt/bulk-data/research_data/USC_beamtime/APS_March_2016/multi_sample/test'
    # key_list = [f for f in os.listdir(folder) if f.endswith('.gr')]
    # key_list = [f for f in os.listdir(folder) if f.endswith('.fq')]
    key_list = [f for f in os.listdir(folder) if f.endswith('.chi')]
    key_list.sort()
    # key_list = key_list[::20]
    key_list = key_list
    print(len(key_list))
    if key_list[0].endswith('.gr') or key_list[0].endswith('.fq'):
        skr = 27
    else:
        skr=4
    data_list = [(np.loadtxt(os.path.join(folder, f), skiprows=skr)[:, 0],
                  np.loadtxt(os.path.join(folder, f), skiprows=skr)[:, 1]) for f
                 in key_list]

    app = QtGui.QApplication(sys.argv)
    tt = LineStackExplorer(key_list, data_list)
    tt.show()
    sys.exit(app.exec_())
