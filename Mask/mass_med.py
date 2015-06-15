"""
This program is designed to take a MED data configuration file, which holds
the background, mask, base file name, directory, and other information about
the system and produce flat files.  The flat files are further used in
modulation enhanced diffraction analysis of pair distribution functions.
"""
__author__ = 'christopher'


from Push_Button_Get_Figure import *
from tkFileDialog import askdirectory
import os
import matplotlib.pyplot as plt
import time
import fnmatch
import ConfigParser
import pprint
import json


#perform the push button routiene on the first phase of each of the MED
# dictionaries, populate the qmins, qmax, in all the others, start all the
# others

#Get information on the experiments from a configuration file
cf = ConfigParser.SafeConfigParser()
cf.optionxform = str
configfile_str = '/media/christopher/Seagate_Backup_Plus_Drive/Beamtime' \
                 '/July2014/MED/MED.conf'
configfile = open(configfile_str)
cf.readfp(configfile)
print cf.sections()

#build list of MED dictionaries
multi_med_list = []

# directory, background file, mask, composition, list of phases
for section in cf.sections():
    med_dict = {}
    med_dict['directory'] = os.path.join(os.path.split(configfile_str)[0], section)
    med_dict['background_file'] = cf.get(section, 'background_file')
    med_dict['composition'] = cf.get(section, 'composition')
    med_dict['poni_file'] = cf.get(section, 'poni_file')
    med_dict['mask_file'] = cf.get(section, 'mask_file')
    med_dict['resolution'] = float(cf.get(section, 'resolution'))
    med_dict['start'] = int(cf.get(section, 'start'))
    med_dict['stop'] = int(cf.get(section, 'stop'))
    med_dict['images per cycle'] = int(cf.get(section, 'images per cycle'))
    med_dict['base_file_name'] = cf.get(section, 'base_file_name')
    med_dict['phase_dictionary_list'] = []

    # list of images in phase, qmin, qmax,qmaxinst, bgscale,
    for i in range(med_dict['start'], med_dict['start'] + med_dict['images '
                                                                   'per cycle'] + 1,
                   1):
        phase_dictionary = {}
        phase_file_list = []
        phase_dictionary['qmin'] = 2 if cf.has_option(section, 'qmin') is \
                                        False else float(cf.get(section,
                                                                'qmin'))
        phase_dictionary['qmax'] = None if cf.has_option(section, 'qmax') is \
                                        False else float(cf.get(section,
                                                                'qmax'))
        phase_dictionary['qmaxinst'] = 19.5 if cf.has_option(section,
                                                             'qmaxinst') is \
                                        False else float(cf.get(section,
                                                             'qmaxinst'))
        phase_dictionary['bg_scale'] = None if cf.has_option(section,
                                                             'bg_scale') is \
                                        False else float(cf.get(section,
                                                             'bg_scale'))
        phase_dictionary['bg_refine_qmax_index'] = None if cf.has_option(
            section, 'bg_refine_qmax_index') is False else int(
            cf.get(section, 'bg_refine_qmax_index'))
        phase_dictionary['bg_refine_qmin_index'] = None if cf.has_option(
            section, 'bg_refine_qmin_index') is False else int(
            cf.get(section, 'bg_refine_qmin_index'))
        for j in range(i, med_dict['stop'] + 1, med_dict['images per cycle']):
            file_name = med_dict['base_file_name'] + str(j).zfill(5) + '.tif'
            full_file_name = os.path.join(med_dict['directory'], file_name)
            if os.path.exists(full_file_name):
                phase_file_list.append(full_file_name)
        phase_dictionary['file_list'] = phase_file_list
        med_dict['phase_dictionary_list'].append(phase_dictionary)
    multi_med_list.append(med_dict)

#now we initialize all the experiments by examining their first phase
# pprint.pprint(multi_med_list)
#
initial_list = []
for med_dictionary in multi_med_list:
    if med_dictionary['phase_dictionary_list'][0]['qmax'] is None:
        # pprint.pprint(med_dictionary)
        chi_name = tif_to_iq(
            poni_file=med_dictionary['poni_file'],
            image_path=med_dictionary['phase_dictionary_list'][0]['file_list'],
            mask_file=med_dictionary['mask_file'],
            # plot_verbose=True,
            resolution=med_dictionary['resolution'],
            # solid=True
            )

        print 'Start PDF'

        gr_file, bg, calculated_qmax, bg_refine_qmin_index, \
        bg_refine_qmax_index = write_pdf(
            chi_file=chi_name,
            composition=med_dictionary['composition'],
            qmax='ripple',
            qmin=med_dictionary['phase_dictionary_list'][0]['qmin'],
            qmaxinst=med_dictionary['phase_dictionary_list'][0]['qmaxinst'],
            background_file=med_dictionary['background_file'],
            bg_refine_qmax_index=med_dictionary['phase_dictionary_list'][0]['bg_refine_qmax_index'],
            bg_refine_qmin_index=med_dictionary['phase_dictionary_list'][0]['bg_refine_qmin_index'],
            plot_verbose=True
        )

        out_tuple = calculated_qmax, bg_refine_qmax_index, bg_refine_qmin_index
        initial_list.append(out_tuple)

pprint.pprint(initial_list)

# deciminate the information from the first phase to the other phases in each
# experiment
for med_dictionary, phase_tuple in zip(multi_med_list, initial_list):
    for number, phase in enumerate(med_dictionary['phase_dictionary_list']):
        med_dictionary['phase_dictionary_list'][number]['qmax'] = phase_tuple[0]
        med_dictionary['phase_dictionary_list'][number][
            'bg_refine_qmax_index'] = phase_tuple[1]
        med_dictionary['phase_dictionary_list'][number][
            'bg_refine_qmin_index'] = phase_tuple[2]

#giant pair of loops that run the med for each phase
flat_file_list = []
for med_dictionary in multi_med_list:
    gr_array = None
    for phase_dictionary in med_dictionary['phase_dictionary_list']:
        chi_name = tif_to_iq(
            poni_file=med_dictionary['poni_file'],
            image_path=phase_dictionary['file_list'],
            mask_file=med_dictionary['mask_file'],
            resolution=med_dictionary['resolution'],
            # solid=True
            )

        print 'Start PDF'

        gr_file, bg, calculated_qmax, bg_refine_qmin_index, \
        bg_refine_qmax_index = write_pdf(
            chi_file=chi_name,
            composition=med_dictionary['composition'],
            qmax=phase_dictionary['qmax'],
            qmin=phase_dictionary['qmin'],
            qmaxinst=phase_dictionary['qmaxinst'],
            background_file=med_dictionary['background_file'],
            bg_refine_qmax_index=phase_dictionary['bg_refine_qmax_index'],
            bg_refine_qmin_index=phase_dictionary['bg_refine_qmin_index'],
            )
        radius, gr = IO.load_gr_file(gr_file)
        if gr_array is None:
            gr_array = radius
        gr_array = np.c_[gr_array, gr]
    IO.write_flat_file(gr_array, os.path.join(med_dictionary['directory'],
                                              med_dictionary['base_file_name']))
    flat_file_list.append(tuple([med_dictionary['base_file_name'], gr_array]))

with open('/media/christopher/Seagate_Backup_Plus_Drive/Beamtime/July2014'
           '/MED/data.json', 'w') as jsonfile:
    json.dump(multi_med_list, jsonfile, encoding='utf8')

for base_file_name, gr_array in flat_file_list:
    plt.imshow(gr_array[:, 1:].transpose(), aspect='auto',
               interpolation='nearest', cmap='cubehelix', extent=[0, 30, 0, 8])
    plt.colorbar()
    plt.title('Sum')
    plt.show()