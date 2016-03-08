__author__ = 'christopher'


import os

search_dir = '/media/christopher/Seagate_Backup_Plus_Drive/Beamtime/July2014/MED/121_3_2_NiPd_MED_RT'
new_path = '/media/christopher/Seagate_Backup_Plus_Drive/Beamtime/July2014/MED/121_3_2_NiPd_MED_RT_CORRECT'
base_file_name = '121_3_2_NiPd_MED_RT-'

os.chdir(search_dir)
files = filter(os.path.isfile, os.listdir(search_dir))
files = [os.path.join(search_dir, f) for f in files] # add path to each file
files.sort(key=lambda x: os.path.getmtime(x))


exclude_list = []
for j in range(0, 135, 1):
    exclude_list.append(str(j).zfill(5))

i = 0
for totalfile in files:
    afile_path, afile = os.path.split(totalfile)
    afile_base_and_number, afile_ext = afile.split('.', 1)
    afile_base, afile_number = afile_base_and_number.split('-', 1)
    if afile_number not in exclude_list and afile_ext != 'tif.metadata':
        print totalfile, '->', os.path.join(new_path, afile_base + '-' + str(
            135+i).zfill(5) + '.' + afile_ext)
        print totalfile + '.metadata', '->', os.path.join(new_path, afile_base +
                                                          '-' + str( 135+i).zfill(5) + '.' + afile_ext + '.metadata')

        os.rename(totalfile, os.path.join(new_path, afile_base + '-' + str(
            135+i).zfill(5) + '.' + afile_ext))
        os.rename(totalfile + '.metadata', os.path.join(new_path, afile_base + '-' +
                                                        str( 135+i).zfill(5) + '.' + afile_ext + '.metadata'))

        i += 1
#
# for totalfile in files:
#     path, afile = os.path.split(totalfile)
#     if os.path.splitext(afile)[1] == '':
#
#         afile_number_and_tif = afile.split('121_3_2_NiPd_MED_RT')[1]
#         afile_number = afile_number_and_tif.split('tif')[0]
#
#         new_afile = base_file_name + afile_number +'.tif'
#         print new_afile
#         os.rename(totalfile, os.path.join(path, new_afile))
#
#     if os.path.splitext(afile)[1] == '.metadata':
#         afile_no_ext, ext = os.path.splitext(afile)
#         if len(afile_no_ext.split('.')) == 1:
#             print afile_no_ext
            # afile_number_and_tif = afile_no_ext.split('121_3_2_NiPd_MED_RT')[1]
            # print afile_number_and_tif
            # afile_number = afile_number_and_tif.split('tif')[0]
            # print afile_number
            #
            # new_afile = base_file_name + afile_number +'.tif' + ext
            # print new_afile
            # os.rename(totalfile, os.path.join(path, new_afile))