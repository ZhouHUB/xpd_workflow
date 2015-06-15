"""
Copyright (c) 2014 Brookhaven National Laboratory All rights reserved.
Use is subject to license terms and conditions.

@author: Christopher J. Wright"""
__author__ = 'Christopher J. Wright'

from Push_Button_Get_Figure import *
from tkFileDialog import askdirectory
import os
import matplotlib.pyplot as plt
import time
start_time = time.time()

experiment = []
# poni_file = '/media/christopher/My Passport/July2014/121_0_0_Ni_STD/121_RT_Sample2_FinalSum.poni'
poni_file = '/media/christopher/My Passport/Mar2014/Ni/17_21_24_NiSTD_300K-00002.poni'
# mask_name = '/media/christopher/My Passport/July2014/Ni_STD-00001.msk'
mask_name = '/media/christopher/My Passport/Mar2014/Mastermask4.msk'
base_file_name = '7_7_8_NiPd_AirAnneal_438K-'
# base_file_name = 'KJ_Carbon-'
# base_file_name = '121_RT_Sample2-'
# base_file_name = 'Ni_STD-'
# tif_directory = askdirectory()
# tif_directory='/media/christopher/My Passport/July2014/121_0_0_Ni_STD'
tif_directory='/media/christopher/My Passport/July2014/121_2_2_NiPD_S_RT'
tif_directory='/media/christopher/My Passport/Mar2014/'
os.chdir(tif_directory)
i = 0
start = 0
stop = 80
# start = 23
# stop = 344
step = None
#Just sum the files
if step is None:
    for i in range(start, stop+1, 1):
        file_chunk = (base_file_name + str(i).zfill(5) + '.tif')
        if os.path.exists(file_chunk):
            print file_chunk
            experiment.append(file_chunk)
        sum_time_start = time.time()
    print '\nSumming %s through %s' % (experiment[0], experiment[-1])
    chi_name = tif_to_iq(
    poni_file=poni_file,
    image_path=experiment,
    mask=mask_name,
    out_ending=''
    )

    # plt.ioff()
    # plt.clf()
    # print chi_name
    print 'Start PDF'
    gr_file, background_level, qmax, bg_refine_qmin_index, bg_refine_qmax_index = write_pdf(
        chi_file=chi_name,
        background_file='/media/christopher/My '
                            'Passport/July2014/121_1_1_KJ_Carbon/KJ_Carbon-00344_Sum.chi',
        composition='Ni.3Pd.7', qmax='ripple', plot_verbose=True,
        qmaxinst=18.9)
    # radius, gr = IO.load_gr_file(gr_file)
    # print radius
    # print gr
    # plt.clf()
    # plt.plot(radius, gr)
    # plt.show()
    print '\n\n\n'
    print((time.time() - start_time)/60), 'minutes'

else:
    for i in range(start, stop+1, step):
        # print i
        file_chunk = []
        j = i
        for j in range(j, j + step, 1):
            file_chunk.append(base_file_name + str(j).zfill(5) + '.tif')
        experiment.append(file_chunk)
    bg_list=[]
    qmax_list=[]
    time_list=[]
    for sum_list in experiment:
        sum_time_start = time.time()
        print '\nSumming %s through %s' % (sum_list[0], sum_list[-1])
        chi_name = tif_to_iq(
        poni_file=poni_file,
        image_path=sum_list,
        mask=mask_name,
        out_ending='_q_25'
        )

        # plt.ioff()
        # plt.clf()
        # print chi_name
        print 'Start PDF'
        gr_file, bg, calculated_qmax = write_pdf(
            chi_file=chi_name,
            background_file='',
            composition='Ni.3Pd.7', qmax=25)
        bg_list.append(bg)
        qmax_list.append(calculated_qmax)

        # radius, gr = IO.load_gr_file(gr_file)
        # print radius
        # print gr
        # plt.clf()
        # plt.plot(radius, gr)
        # plt.show()
        time_list.append((time.time()-sum_time_start)/60)
    print '\n\n\n'
    print((time.time() - start_time)/60), 'minutes'

    bg_array=np.array(bg_list)
    qmax_array=np.array(qmax_list)
    time_array=np.array(time_list)

    plt.plot(time_array)
    plt.show()


    plt.plot(bg_array)
    plt.show()


    plt.plot(qmax_array)
    plt.show()