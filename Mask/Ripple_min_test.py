'''
Copyright (c) 2014 Brookhaven National Laboratory All rights reserved. 
Use is subject to license terms and conditions.

@author: Christopher J. Wright'''
__author__ = 'Christopher J. Wright'

from IO import load_gr_file
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
from Push_Button_Get_Figure import write_pdf
import os
import scipy.signal as signal


def Main():
    os.chdir('/home/xpdlabuser/Push_Button_Tests/ripple2')
    ripple_list = []
    qmax_min = 15
    qmax_max = 27
    qmax_step = .01
    for qmax in np.arange(qmax_min, qmax_max, qmax_step):
        gr_file, bg, q=write_pdf(
            chi_file='0_0_0_KJCarbon_300K-00010_Summedian.chi',
                  background_file='',
                  background_level=2.29,
                  pdf_format='QA', output='@r_ripple_test_'+str(qmax)+'.@o',
                  qmax=qmax,
                  composition='Ni.3Pd.7', qmaxinst=27.1,
                  qmin=0, relative_max_uncertainty=.90,
                  plot_verbose=False,)
        gr_file = os.path.splitext(gr_file)[0]+'_'+str(qmax)+'.gr'
        x, y = load_gr_file(gr_file)
        w = y - np.convolve(y, np.ones(3)/3, 'same')
        ripple_sum = np.sum(abs(w))
        ripple_list.append(ripple_sum)
    t = np.arange(qmax_min, qmax_max, qmax_step)
    ripple_array=np.array(ripple_list)
    minima_index = signal.argrelextrema(ripple_array, np.less_equal,
                                        order=15)[0]
    minima = ripple_array[minima_index]
    minima_q = t[minima_index]

    maxima_index = signal.argrelextrema(ripple_array, np.less_equal,
                                        order=15)[0]
    maxima = ripple_array[maxima_index]
    maxima_q = t[maxima_index]

    plt.plot(t, ripple_array, 'g', minima_q, minima, 'o',
             markerfacecolor='none', mec='r', maxima_q, maxima, 'o',
             markerfacecolor='none', mec='b')
    plt.title('Ripples in PDF')
    plt.xlabel('Qmax (1/A)')
    plt.ylabel('Ripple Cumulative Sum')
    plt.show()
    return np.array(ripple_list)

if __name__ == "__main__":
    Main()