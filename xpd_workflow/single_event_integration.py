from __future__ import division, print_function

import sys
from xray_vision import QtGui, QtCore
from vis1d import LineStackExplorer
from databroker import db, get_events
from databroker.databroker import fill_event
from metadatastore.api import db_connect as mds_db_connect

from filestore.api import db_connect as fs_db_connect
from filestore.api import retrieve
from sidewinder_spec.utils.handlers import *
from xpd_workflow.mask_tools import *
from diffpy.pdfgetx import PDFGetter
from datamuxer import DataMuxer


def single_mask_integration(img, geo, tmsk=None):
    img /= geo.polarization(img.shape, .95)

    r = geo.rArray(img.shape)
    fq = geo.qArray(img.shape).ravel()
    fimg = img.ravel()
    bins = 2000

    # Pre masking data
    bs_kwargs = {'bins': bins,
                 'range': [0, fq.max()]}

    if not tmsk:
        tmsk = np.ones(img.shape, dtype=int).ravel().astype(bool)
    tmsk *= mask_edge(img.shape, 30)
    tmsk *= ring_blur_mask(img, r, geo.pixel1, (4.5, 2.),
                           mask=tmsk)

    fmsk_img = fimg[tmsk]
    fmsk_q = fq[tmsk]

    # Post masking data
    bs_args = (fmsk_q, fmsk_img)
    median = sts.binned_statistic(*bs_args, statistic='median',
                                  **bs_kwargs)
    x = median[1][:-1]
    y = median[0]
    return x, y


def single_event_workflow(foreground_event, background_event, geos):
    fg_x, fg_y = single_mask_integration()
    # Mask/Integrate the foreground
    # Mask/Integrate the background
    # Subtract the background
    # PDFgetx3
    pass