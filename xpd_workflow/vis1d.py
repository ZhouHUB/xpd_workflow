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
    def __init__(self, key_list, data_list, name='StackExplorer', parent=None):
        QtGui.QMainWindow.__init__(self, parent=parent)
        self.setWindowTitle(name)
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

    folder = '/mnt/bulk-data/research_data/USC_beamtime/APS_March_2016/S1/temp_exp'
    # key_list1 = [f for f in os.listdir(folder) if f.endswith('.gr')]
    # key_list2 = [f for f in os.listdir(folder) if f.endswith('.fq')]
    key_list1 = [f for f in os.listdir(folder) if
                 f.endswith('.chi') and not f.startswith('d') and f.strip('0.chi') != '' and int(
                     f.lstrip('0').strip('.chi')) % 2 == 1]
    key_list1.sort()
    print(key_list1)
    # if key_list1[0].endswith('.gr') or key_list2[0].endswith('.fq'):
    #     skr = 27
    # elif key_list1[0].endswith('.chi'):
    #     skr = 4
    data_list1 = [(np.loadtxt(os.path.join(folder, f),
                              skiprows=8
                              )[:, 0],
                   np.loadtxt(os.path.join(folder, f),
                              skiprows=8
                              )[:, 1])
                  for f in key_list1]

    app = QtGui.QApplication(sys.argv)
    tt = LineStackExplorer(key_list1, data_list1)
    tt.show()
    sys.exit(app.exec_())
