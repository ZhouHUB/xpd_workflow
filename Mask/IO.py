"""
Copyright (c) 2014 Brookhaven National Laboratory All rights reserved.
Use is subject to license terms and conditions.

@author: Christopher J. Wright"""
__author__ = 'Christopher J. Wright'

import pyFAI
import fabio
import tkFileDialog
import os
import numpy as np
import subprocess
import math


def loadgeometry(ponifilepath=None):
    """Loads the detector geometry information from a poni file or from hard coded information
    Parameters
    ----------
    ponifilepath: str
        File path to the .poni file
    Returns
    -------
    Azimuthal integrator object:
        an object which contains the geometry information for the detector
    """
    if ponifilepath is None:
        filepath = tkFileDialog.askopenfilename()
        a = pyFAI.load(ponifilepath)
    else:
        a = pyFAI.load(ponifilepath)
    return a


def loadimage(filepath=None):
    """
    Loads images

    Parameters
    ----------
    filepath: str or list of str or tuple of str
        If str, load image, else if iterable of str load images and sum them
        together

    Returns
    -------
    ndarray:
        The image
    str:
        Name of the chi file
    str:
        File(s) loaded

    """
    if filepath is None:
        print 'Open Image'
        filepath = tkFileDialog.askopenfilenames()
        print filepath

    if type(filepath) is tuple or type(filepath) is list:
        xray_image = None
        for path in filepath:
            if path.endswith('.tif'):
                if xray_image is None:
                    if os.path.exists(path):
                        try:
                            xray_image = fabio.open(path).data
                        except AttributeError:
                            print 'File did not load, please double ' \
                                'check file.  Continuing without ' +path
                            pass
                    else:
                        print '%s not found' % (path)

                else:
                    try:
                        xray_image += fabio.open(path).data
                    except AttributeError:
                        print 'File did not load, please double ' \
                            'check ' \
                        'file.  Continuing ' \
                        'without ' +path
        if len(filepath) >1:
            file_path_name = os.path.splitext(filepath[-1])[0]+'_Sum'
        else:
            file_path_name = os.path.splitext(filepath[-1])[0]

    else:
        print filepath
        xray_image = fabio.open(filepath).data
        file_path_name = os.path.splitext(filepath)[0]
    return xray_image, file_path_name, filepath


def save_chi(integrated, geometry_object, filename, mask_file_name,
             input_file):
    """
    Save the I(Q) data as a chi file

    Parameters
    ----------
    integrated: ndarray
        The integrated data
    geometry_object: Azimuthal integrator object
        The object which holds the detector geometry used
    filename: str
        The name of the file to be saved
    mask_file_name: str
        The name of the mask file(s) used in the integration
    input_file: str
        Name(s) of the input images

    Returns
    -------
    str:
        The output file path
    """
    filepath = filename + '.chi'
    f = open(filepath, 'wb')
    f.write('pyFAI Parameters\n')
    for key, value in geometry_object.getPyFAI().items():
        f.write(key +' = ' + str(value)+'\n')

    if type(input_file) is tuple or type(input_file) is list:
        f.write('\nFile Name list:\n')
        for item in input_file:
            f.write(str(item)+'\n')
    else:
        f.write('\nFile Name = ' + str(input_file)+'\n')

    if type(mask_file_name) is tuple or type(mask_file_name) is list:
        f.write('\nMask File list:\n')
        for item in mask_file_name:
            f.write(str(item))
    else:
        f.write('\nMask File = ' + str(mask_file_name) +'\n \n')

    f.write('START DATA WITH HEADER\n')
    f.write('Q/TTH I[Q/TTH] sigma_I[Q/TTH] \n')
    np.savetxt(f, integrated, fmt='%g')
    f.close()
    return filepath


def load_chi_file(chi_file, skiplines=None):
    """
    Load chi file

    Parameters
    ----------
    chi_file: str
        Location of the chi file
    skiplines: int
        Number of lines needed to skip to get to the data, if none start
        reading at START DATA WITH HEADER

    Returns
    -------
    ndarray:
        The data
    str:
        The chi file name
    """
    if chi_file is None:
        print 'Open CHI'
        chi_file = tkFileDialog.askopenfilename()
        if chi_file is '' or chi_file is ():
            return None, None
    elif skiplines is None:
        with open(chi_file) as my_file:
            for num, line in enumerate(my_file,1):
                if 'START DATA WITH HEADER' in line:
                    skiplines=num+1
                    break
    data = np.loadtxt(chi_file, skiprows=skiplines)
    return data, chi_file


def loadmask(mask=None):
    """

    """
    if mask is None:
        print 'Open Mask'
        mask = tkFileDialog.askopenfilename()
    if mask == '':
        user_mask = None
        return None
    else:
        user_mask = fabio.fit2dmaskimage.fit2dmaskimage()
        user_mask = user_mask.read(mask).data
        flipped_user_mask = np.flipud(user_mask)
        invert_user_mask = np.abs(flipped_user_mask-1).astype('bool')
    return invert_user_mask

def fit2d_save(Mask, image_file_name):
    """
    Compresses and wraps the mask for Fit2D use
    """
    Mask = np.flipud(Mask)
    Mask = np.abs(Mask-1)
    print image_file_name
    path, imagenamext=os.path.split(image_file_name)
    image_name = os.path.splitext(imagenamext)[0]
    if path is not '':
        os.chdir(path)
    currentdir = os.getcwd()
    Maskname = os.path.join(currentdir, image_name + ".msk")
    b1 = Mask.shape[0]
    b2 = Mask.shape[1]
    b28 = b2 / 8
    C = np.zeros((b1, b28), dtype=np.uint8)
    print 'start compression'
    # bitmap compression for Fit2D io
    # TODO: add to pyFAI via github fork/pull request maybe?
    for i in range(0, b1):
        j = 0
        while j < b2:
            C[i, math.floor(j / 8)] = Mask[i, j]
            for k in range(1, 8):
                C[i, math.floor(j / 8)] += Mask[i, j + k] * (2 ** k)
            j += 8
    fout = open(Maskname, "wb")
    C.tofile(fout, "", "%x")
    print('mask compressed')
    fout.close()
    # Giant long mask header, must be there for fit2d to work with the mask
    header = 'M\x00\x00\x00A\x00\x00\x00S\x00\x00\x00K\x00\x00\x00\x00\x08\x00\x00\x00\x08\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'
    mask = open(Maskname, 'rb')
    bodypt1 = mask.read(262144)
    bodypt2 = mask.read()
    mask.close()
    # fin.close()
    body = [bodypt1, bodypt2]
    body2 = "".join(body)
    total = [header, body2]
    total2 = "".join(total)
    fout = open(Maskname, "wb")
    fout.write(total2)
    fout.close()
    print('mask finalized')
    return Maskname

def calibrate(poni_file=None, calibration_imagefile=None, dspace_file=None,
              wavelength=None, detector=None):
    """
    Load the calibration file/make a calibration file
    """

    pwd = os.getcwd()
    if poni_file is None:
        print 'Open PONI'
        poni_file = tkFileDialog.askopenfilename(defaultextension='.poni')
        if poni_file == '':
            # if prompt fail run calibration via pyFAI
            if calibration_imagefile is None:
                calibration_imagefile = tkFileDialog.askopenfilename(
                    defaultextension='.tif')
            if dspace_file is None:
                dspace_file = tkFileDialog.askopenfilename()
            detector = detector if detector is not None else raw_input(
                'Detector Name: ')
            wavelength = wavelength if wavelength is not None else raw_input(
                'Wavelength in Angstroms: ')
            calibration_image_dir, calibration_image_name = os.path.split(
                calibration_imagefile)
            os.chdir(calibration_image_dir)
            subprocess.call(['pyFAI-calib',
                             '-w', str(wavelength),
                             '-D', str(detector),
                             '-S', str(dspace_file),
                             str(calibration_image_name)])

            # get poni_file name make new poni_file
            poni_file = os.path.splitext(calibration_imagefile)[0] + '.poni'
    os.chdir(pwd)
    a = loadgeometry(poni_file)
    return a


def load_gr_file(gr_file=None, skiplines=None):
    #TODO: also give back the filename
    if gr_file is None:
        print 'Open Gr'
        gr_file = tkFileDialog.askopenfilename()
    if skiplines is None:
        with open(gr_file) as my_file:
            for num, line in enumerate(my_file,1):
                if '#### start data' in line:
                    skiplines=num+2
                    break
    data = np.loadtxt(gr_file, skiprows=skiplines)
    return data[:, 0], data[:, 1]


def write_flat_file(flat_array, filename):
    with open(filename + '.flt','w') as f:
        np.savetxt(f, flat_array)
        f.close()
    return filename