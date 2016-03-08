"""
Copyright (c) 2014 Brookhaven National Laboratory All rights reserved.
Use is subject to license terms and conditions.

@author: Christopher J. Wright"""

__author__ = 'Christopher J. Wright'
import IO
import Calculate
import numpy as np
import scipy.signal as signal
import os
import matplotlib.pyplot as plt

plt.ioff()
#TODO: There is a problem here with managing how IO responds and the
# difference between not knowing the file position and not wanting to load a
#  background file, and the appropriate returns


def background_subtract(sample, background_file=None, background_level=None,
                        bg_refine_qmin_index=None, bg_refine_qmax_index=None):
    """
    Wrap the background subtraction optimization with file IO

    Parameters
    ----------
    :param bg_refine_qmin_index:
    sample: ndarray
        Sample I(Q) array
    background_file: str, optional
        The file name of the background file, if not specified give an IO
        gui to select the file, f '' no file is loaded
    background_level: float, optional
        The background level if none given, calculate via minimization

    Returns
    -------
    float:
        The background level
    array:
        The background array
    array:
        The background q array
    str:
        The name of the background file

    """
    #File IO
    background_data, background_file = IO.load_chi_file(
        background_file)
    #Unpack file IO
    if background_data is not None:
        background = background_data[:, 1]
        background_q = background_data[:, 0]
        percentile_peak_pick = np.where(
            background >= np.percentile(background,  98))[0]
        print len(percentile_peak_pick)
        print percentile_peak_pick[0]
        print percentile_peak_pick[-1]
        if len(percentile_peak_pick) >= 2:
            bg_refine_qmin_index = percentile_peak_pick[0]
            bg_refine_qmax_index = percentile_peak_pick[-1]
        else:

            if bg_refine_qmax_index is None:
                plt.plot(background)
                plt.plot(sample)
                plt.show()
                bg_refine_qmax_index = int(raw_input('High end of the '
                                                     'background to fit: '))
            if bg_refine_qmin_index is None:
                plt.plot(background)
                plt.plot(sample)
                plt.show()
                bg_refine_qmin_index = int(raw_input('Low end of the '
                                                     'background to fit: '))
        # Fit the background
        if background_level is None:
            background_level = Calculate.optimize_background_level(background,
                                                                   sample,
                                                                   bg_refine_qmin_index,
                                                                   bg_refine_qmax_index)
            print 'background level= ', background_level
    else:
        background = None
        background_q = None
    return background_level, background_q, background, background_file, \
           bg_refine_qmin_index, bg_refine_qmax_index


def qmax_statistics(sample_q, sample, uncertainty, qmax, qmin, qmaxinst,
              relative_max_uncertainty,):
    """

    """
    #TODO: Include background statistics in this?  
        
    if qmax is 'statistics':
        old_settings = np.seterr(divide='ignore')
        uncertainty_percent_array = uncertainty / sample * 100
        np.seterr(**old_settings)
    
        uncertainty_percent_array[np.isnan(uncertainty_percent_array)] = 0
        uncertainty_percent_array[np.isinf(uncertainty_percent_array)] = 0
        qmin_index = np.where(sample_q == qmin)[0]
        uncertainty_first_peak = signal.argrelmin(uncertainty)[0][0]
        if qmin_index > uncertainty_first_peak:
            starting_q = qmin_index
        else:
            starting_q = uncertainty_first_peak

        if qmaxinst is 0:
            qmaxinst = float(np.amax(sample_q))
        qmaxinst_index = np.where(sample_q == qmaxinst)[0]

        if uncertainty_percent_array[qmaxinst_index] > 100:
            qmaxinst_index = np.where(uncertainty_percent_array[starting_q:]
                                      == 100)

        print qmaxinst_index


        # find the first position where the uncertainty is above the threshold
        try:
            uncertainty_max = np.max(uncertainty_percent_array[
                                     starting_q:qmaxinst_index])
            for arr_index, value in enumerate(uncertainty_percent_array[
                                              starting_q:qmaxinst_index]):
                relative_uncertainty = value / uncertainty_max
                if relative_uncertainty >= relative_max_uncertainty:
                    qmax_index = arr_index + starting_q
                    break
        except StopIteration:
            print 'qmax_error'
            print relative_max_uncertainty
        print sample_q
        qmax = float(sample_q[qmax_index])

    print 'q=', qmax
    return qmax
    

def generate_composition():
    print 'Add composition'
    element_list = []
    amount_list = []
    while True:
        element = raw_input('Element: ')
        if element == '':
            break
        amount = raw_input('Amount: ')
        element_list.append(element)
        amount_list.append(amount)
    for element, amount in zip(element_list, amount_list):
        single_composition = str(element) + str(amount)
        if composition is None:
            composition = single_composition
        else:
            composition += single_composition
    return composition


def ripple_minima(chi_file, background_file, background_level, composition,
                  qmin, qmaxinst):
    pwd = os.getcwd()
    if os.path.exists(os.path.join(pwd, 'ripple_minima_test')) is False:
        os.mkdir(os.path.join(pwd, 'ripple_minima_test'))
    ripple_list = []
    gr_list = []
    qmax_min = 18.5
    qmax_max = qmaxinst
    qmax_step = .01
    print 'Start qmax refinement'
    for qmax in np.arange(qmax_min, qmax_max, qmax_step):
        print qmax
        gr_file, bg, q=Calculate.pdfgetx3(
            chi_file,
            background_file,
            background_level,
            composition,
            qmax,
            qmin,
            qmax,
            'QA',
            os.path.join(pwd, 'ripple_minima_test', '@b_ripple_test_'+str(
                qmax).ljust(5,'0')+'.@o'))
        gr_file = os.path.join(pwd, 'ripple_minima_test', os.path.split(
            os.path.splitext(chi_file)[0])[1]+'_ripple_test_'+str(qmax).ljust(5,'0')+'.gr')
        x, y = IO.load_gr_file(gr_file)
        w = y - np.convolve(y, np.ones(3)/3, 'same')
        ripple_sum = np.sum(abs(w))
        ripple_list.append(ripple_sum)
        gr_list.append(gr_file)
    t = np.arange(qmax_min, qmax_max, qmax_step)
    ripple_array=np.array(ripple_list)
    minima_index = signal.argrelextrema(ripple_array, np.greater_equal,
                                        order=2*len(ripple_array)/100)[0]
    minima = ripple_array[minima_index]
    minima_q = t[minima_index]

    maxima_index = signal.argrelextrema(ripple_array, np.less_equal,
                                        order=2*len(ripple_array)/100)[0]
    maxima = ripple_array[maxima_index]
    maxima_q = t[maxima_index]

    plt.plot(t, ripple_array, 'g')
    plt.plot(minima_q, minima, 'o', markerfacecolor='none', mec='r')
    plt.plot(maxima_q, maxima, 'o', markerfacecolor='none', mec='b')
    plt.title('Ripples in PDF')
    plt.xlabel('Qmax (1/A)')
    plt.ylabel('Ripple Cumulative Sum')
    plt.show()
    plt.cla()

    # data_list = []
    # key_list = []
    # for minima_q_point in minima_q:
    #     data = IO.load_gr_file([x for x in gr_list if str(qmax) in x] [0])
    #     data_list_element = data
    #     key_list_element = minima_q_point
    #     data_list.append(data_list_element)
    #     key_list.append(key_list_element)
    # plot_stack(data_list, key_list)

    ziped_minima = zip(minima_q, minima)
    while True:
        for q, value in ziped_minima:
            print 'q= ', q,'rippleness= ', value
        qmax= float(raw_input('Please pick a qmax.  qmax= '))
        gr_file = [x for x in gr_list if str(qmax) in x] [0]
        print gr_file
        x, y = IO.load_gr_file(gr_file)
        plt.plot(x, y)
        plt.show()
        if raw_input('Is this the file you wanted? y/[n] ') == 'y':
            break
    return qmax


# def plot_stack(data_list, key_list):
#     from vistools.qt_widgets import Stack1DMainWindow
#     print 'hi'
#     mw = Stack1DMainWindow(data_list=data_list, key_list=key_list)
#     print 'showing'
#     mw.show()


