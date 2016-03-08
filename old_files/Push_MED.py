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
import fnmatch

plt.ioff()
poni_file = '/media/christopher/Seagate_Backup_Plus_Drive/Beamtime/July2014' \
            '/Static/121_0_0_Ni_STD/121_RT_Sample2_FinalSum.poni'
# integrate standard, make mask

# tif_directory = '/media/Seagate_Backup_Plus_Drive/Ni_Standard'
# os.chdir(tif_directory)

# std_chi_name, mask_name = tif_to_iq(
#     poni_file='/media/Seagate_Backup_Plus_Drive/Ni_Standard/Ni_STD-00001.poni',
#     image_path='Ni_STD-00001.tif',
#     generate_mask=True, initial_m=10,
#     plot_verbose=True, resolution=0.004
# )

# print 'Start PDF'
# gr_file, bg, calculated_qmax = write_pdf(chi_file='Ni_STD-00001.chi',
#                                          composition='Ni',
#                                          plot_verbose=True,
#                                          qmin=1, qmax='ripple')
#
# radius, gr = IO.load_gr_file(gr_file)
# plt.plot(radius, gr)
# plt.show()

mask_name = '/media/christopher/Seagate_Backup_Plus_Drive/Beamtime/July2014/121_RT_Sample2-00000biggerthresh.msk'
# integrate background, use STD mask
print 'start background'
# tif_directory = '/media/christopher/Seagate_Backup_Plus_Drive/Beamtime' \
#                 '/July2014/Static/Background_July20/'
tif_directory = askdirectory()
os.chdir(tif_directory)

includepattern = ['*.tif']
excludepattern = ['*.dark.tif', '*.raw.tif', '*.tif.metadata']

file_list = os.listdir(tif_directory)
file_set = set()

for includep in includepattern:
    file_set |= set(fnmatch.filter(file_list, includep))
for excludep in excludepattern:
    file_set -= set(fnmatch.filter(file_list, excludep))

tif_file_list = list(file_set)
print tif_file_list

bg_chi_name = tif_to_iq(
    poni_file=poni_file,
    image_path=tif_file_list,
    mask=mask_name,
    plot_verbose=True,
    resolution=0.02,
    # solid=True
)
#
print 'Start PDF'
gr_file, bg, calculated_qmax, bg_refine_qmin_index, \
        bg_refine_qmax_index = write_pdf(chi_file=bg_chi_name,
                                         composition='C SiO2',
                                         plot_verbose=True,
                                         qmin=1, qmax='ripple',
                                         qmaxinst=23.5,
                                         background_file='')

radius, gr = IO.load_gr_file(gr_file)
plt.plot(radius, gr)
plt.show()
bg_chi_name = os.path.join(tif_directory, bg_chi_name)


# #integrate all the chi files, make gr from chi files, no summing, make certain
# #  to make flat file and prep for flat file graphing
# print 'start no sum'
# tif_directory = '/media/Seagate_Backup_Plus_Drive/NiPd_MED/MED_7_7_RT'
# os.chdir(tif_directory)
# base_file_name = '7_7_NiPd_MED_RT-'
#
# includepattern = ['*.tif']
# excludepattern = ['*.dark.tif', '*.raw.tif', '*.tif.metadata']
#
# file_list = os.listdir(tif_directory)
# file_set = set()
#
# for includep in includepattern:
#     file_set |= set(fnmatch.filter(file_list, includep))
# for excludep in excludepattern:
#     file_set -= set(fnmatch.filter(file_list, excludep))
#
# tif_file_list = list(file_set)
# print tif_file_list
#
# gr_array = None
#
# for file_name in tif_file_list:
#     print file_name
#     chi_name = tif_to_iq(
#         poni_file='/media/Seagate_Backup_Plus_Drive/Ni_Standard/Ni_STD-00001.poni',
#         image_path=file_name,
#         mask=mask_name, resolution=0.004
#     )
#     #
#     print 'Start PDF'
#     gr_file, bg, calculated_qmax = write_pdf(chi_file=chi_name,
#                                              composition='Ni.3Pd.7',
#                                              qmin=1, background_file=bg_chi_name)
#     radius, gr = IO.load_gr_file(gr_file)
#     if gr_array is None:
#         gr_array = radius
#         gr_array = np.c_[gr_array, gr]
#     else:
#         gr_array = np.c_[gr_array, gr]
#
# IO.write_flat_file(gr_array, os.path.join(tif_directory,
#                                           base_file_name + 'no_sum'))


#sum and integrate chi files, make gr, get flat file

# bg_chi_name='/media/Seagate_Backup_Plus_Drive/June2014/background200mm' \
#             '/KJC_200mm_RT_1' \
# '-00020_Sum.chi'
# experiment = []
# base_file_name = '121_3_2_NiPd_MED_RT-'
#
# # tif_directory = askdirectory()
# tif_directory = '/media/christopher/Seagate_Backup_Plus_Drive/Beamtime/July2014/MED/121_3_2_NiPd_MED_RT_CORRECT'
# os.chdir(tif_directory)
# i = 0
# start = 1
# stop = 1475
# images_per_cycle = 122
#
# gr_array2 = None
# # bg_refine_qmin_index=200
# # bg_refine_qmax_index=300
# qmax = 19.8 #determined via ripple method
# for i in range(start, images_per_cycle + 1, 1):
#     # print i
#     cycle = []
#     for j in range(i, stop+1, images_per_cycle):
#         cycle.append(base_file_name + str(j).zfill(5) + '.tif')
#     experiment.append(cycle)
#
# for cycle in experiment:
#     print cycle
#     #
#     if gr_array2 is None:
#         chi_name = tif_to_iq(
#             poni_file=poni_file,
#             image_path=cycle,
#             mask=mask_name, resolution=0.004)
#         print 'Start PDF'
#         gr_file, bg, calculated_qmax, bg_refine_qmin_index, \
#         bg_refine_qmax_index = write_pdf(chi_file='7_7_NiPd_MED_RT-00365_Sum.chi',
#                                                  composition='Ni.3Pd.7',
#                                                  qmax=qmax, qmin=1,
#                                                  background_file=bg_chi_name,
#                                                  plot_verbose=True,
#                                                  bg_refine_qmin_index=200,
#                                                  bg_refine_qmax_index=300)
#         radius, gr = IO.load_gr_file(gr_file)
#         gr_array2 = radius
#         gr_array2 = np.c_[gr_array2, gr]
#     else:
#         chi_name = tif_to_iq(
#             poni_file='/media/Seagate_Backup_Plus_Drive/Ni_Standard/Ni_STD-00001.poni',
#             image_path=cycle,
#             mask=mask_name, resolution=0.004)
#         print 'Start PDF'
#         gr_file, bg, calculated_qmax, bg_refine_qmin_index, \
#         bg_refine_qmax_index = write_pdf(
#             chi_file=chi_name,composition='Ni.3Pd.7',
#             qmax=calculated_qmax, qmin=1, background_file=bg_chi_name,
#             bg_refine_qmax_index=bg_refine_qmax_index,
#             bg_refine_qmin_index=bg_refine_qmin_index)
#         print 'background scale= ', bg
#         radius, gr = IO.load_gr_file(gr_file)
#         gr_array2 = np.c_[gr_array2, gr]
#
# IO.write_flat_file(gr_array2, os.path.join(tif_directory,
#                                           base_file_name + 'no_sum'))
# # plt.imshow(gr_array2[1:])
# # plt.title('No Sum')
# # plt.show()
#
# plt.imshow(gr_array2[1:], aspect='auto')
# plt.title('Sum')
# plt.show()