"""
Copyright (c) 2014 Brookhaven National Laboratory All rights reserved.
Use is subject to license terms and conditions.

@author: Christopher J. Wright"""
__author__ = 'Christopher J. Wright'

import subprocess
import os
import matplotlib.pyplot as plt
import IO
import Calculate
import pdf_parameters
from scipy import stats
from scipy.signal import argrelmin
import numpy as np

plt.ioff()


def tif_to_iq(
        # Calibration info
        poni_file=None, dspace_file=None, detector=None,
        calibration_imagefile=None, wavelength=None,
        # image info
        image_path=None, out_ending='',
        # Correction info
        solid=False, pol=True, polarization=0.95,
        # old_files info
        generate_mask=False, mask_file=None, initial_m=None,
        # integration info
        space='Q', resolution=None, average_method='median',
        plot_verbose=False):
    """
    Takes TIFF images to CHI files, including masking, integrating, and CHI
    writing

    Parameters
    ----------
    poni_file: str
        Location of poni file, a detector geometry file produced by pyFAI.
        If none is given a GUI prompt will ask for the file.
    dspace_file : str
        Location of the calibrant d-spacing file.  This is used to generate
        the poni file if no poni file is loaded.
    detector: str
        Name of the detector. This is used to generate
        the poni file if no poni file is loaded.
    calibration_imagefile : str
        Location of the calibration image.  This is used to generate
        the poni file if no poni file is loaded
    wavelength: float
        The wavelength of light in angstrom. This is used to generate
        the poni file if no poni file is loaded.
    image_path: str or list of str or tuple of str
        The image(s) to be masked and integrated.  If more than one image is
        passed via a list or tuple, sum the images together before proceeding.
    out_ending: str
        The ending to append to the output file name before the extension
    solid: bool
        If true, apply solid angle correction as calculated by pyFAI, else do not
    pol: bool
        If true, apply polarization angle correction
    polarization: float
        The polarization of the beam used in the polarization correction
    generate_mask: bool
        If true, generate mask via mask writing algorithum
    mask_file: str
        Location of mask file
    initial_m: float
        Initial guess at ring threshold mask scale factor
    space: {'Q'}
        Eventually will support TTH as well
    resolution:
        The resolution in Q or TTH of the detector
    average_method: str or function
        The function to use in binned statistics to calculate the average
    plot_verbose: bool
        If true, generate a bunch of plots related to statistics of the rings

    Returns
    -------
    str:
        The name of the chi file
    str:
        The name of the mask, if generated
    """

    # load the calibration file/calibrate
    resolution = resolution if resolution is not None else .02
    a = IO.calibrate(poni_file, calibration_imagefile, dspace_file, wavelength,
                     detector)

    # load image
    xray_image, bulk_image_path_no_ext, image_path = IO.loadimage(image_path)

    # mask image, expect out an array of 1s and 0s which is the mask, if the
    # mask is 0 then mask that pixel
    mask_array = None
    if generate_mask is True:

        """This should handle all mask writing functions including binning,
        applying criteria, plotting, and writing the mask out as a numpy
        array where masked pixels are zero.  This should also give the
        option to combine the calculated mask with a previous user
        generated mask."""

        mask_array = Calculate.write_mask(xray_image, a, mask_file, initial_m,
                                          resolution)
        # write mask to FIT2D file for later use
        mask_file_name = IO.fit2d_save(mask_array, bulk_image_path_no_ext)
    elif mask_file is not None:
        # LOAD THE MASK
        mask_array = IO.loadmask(mask_file)
        # mask_array = np.abs(mask_array-1)
    if mask_array is None:
        mask_array = np.ones(np.shape(xray_image))

    if plot_verbose:
        plt.imshow(mask_array, cmap='cubehelix'), plt.colorbar()
        plt.savefig('/home/christopher/mask.png', bbox_inches='tight',
                    transparent=True)
        plt.show()

        plt.imshow(mask_array*xray_image, cmap='cubehelix'), plt.colorbar()
        plt.show()

    # correct for corrections
    corrected_image = Calculate.corrections(xray_image, geometry_object=a,
                                            solid=solid, pol=pol,
                                            polarization_factor=polarization)

    position_matrix = None
    # TODO: Add options for tth
    if space == 'Q':
        position_matrix = a.qArray(np.shape(xray_image)) / 10

    max_position = np.max(position_matrix)
    masked_position = position_matrix * mask_array

    # Integrate using the global average
    independant_out_array = np.arange(0, max_position + resolution, resolution)

    integrated = Calculate.statistics1d(corrected_image,
                                        position_matrix=masked_position,
                                        max=max_position, step=resolution,
                                        statistic=average_method)
    integrated[np.isnan(integrated)] = 0

    uncertainty_array = np.sqrt(
        Calculate.variance1d(corrected_image, position_matrix=masked_position,
                             max=max_position, step=resolution))

    uncertainty_array[np.isnan(uncertainty_array)] = 0

    if plot_verbose is True:
        no_mask_integrated = Calculate.statistics1d(corrected_image,
                                                    position_matrix=position_matrix,
                                                    max=max_position,
                                                    step=resolution,
                                                    statistic=average_method)
        integrated_mean = Calculate.statistics1d(corrected_image,
                                                 position_matrix=masked_position,
                                                 max=max_position,
                                                 step=resolution,
                                                 statistic='mean')
        fig, ax = plt.subplots()
        ax.plot(independant_out_array, integrated, 'g', label='Median I[Q]')
        ax.plot(independant_out_array, integrated_mean, 'b', label='Mean I[Q]')
        ax.plot(independant_out_array, ((integrated - integrated_mean)/((
                                                                            integrated+integrated_mean)/2))*10e6,
                                        'k', label='((Median-Mean)/('
                                                   'Median+Mean)/2)x10^6')
        ax.legend(loc='upper right')
        plt.xlabel('Q (1/A)')
        plt.ylabel('Raw Counts')
        plt.savefig('/home/christopher/integrator.png', bbox_inches='tight',
                    transparent=True)
        plt.show()

        no_mask_uncertainty_percent_array = np.sqrt(
            Calculate.variance1d(corrected_image,
                                 position_matrix=position_matrix,
                                 max=max_position, step=resolution)) / no_mask_integrated
        fig, ax = plt.subplots()
        ax.plot(independant_out_array, integrated, 'g', label='I[Q]')
        ax.plot(independant_out_array, no_mask_integrated, 'k', label='No '
                                                                      'old_files I[Q]')
        ax.plot(independant_out_array, integrated - no_mask_integrated, 'r',
                label='Difference')
        legend = ax.legend(loc='upper right')
        plt.show()

        old_settings = np.seterr(divide='ignore')
        uncertainty_percent_array = uncertainty_array / integrated
        uncertainty_percent_array[np.isnan(uncertainty_percent_array)] = 0
        np.seterr(**old_settings)

        fig, ax = plt.subplots()
        ax.plot(independant_out_array, uncertainty_percent_array * 100, 'g',
                label='Sigma I[Q]')
        plt.xlabel('Q (1/A)')
        plt.ylabel('% uncertainty')
        legend = ax.legend(loc='upper right')
        plt.show()

        ring_std = Calculate.statistics1d(corrected_image,
                                          position_matrix=masked_position,
                                          max=max_position, step=resolution,
                                          statistic=np.std) / integrated

        ring_std_no_mask = Calculate.statistics1d(corrected_image,
                                                  position_matrix=position_matrix,
                                                  max=max_position,
                                                  step=resolution,
                                                  statistic=np.std) / no_mask_integrated
        ring_s = Calculate.statistics1d(corrected_image,
                                        position_matrix=masked_position,
                                        max=max_position, step=resolution,
                                        statistic=stats.skew)

        ring_s_no_mask = Calculate.statistics1d(corrected_image,
                                                position_matrix=position_matrix,
                                                max=max_position,
                                                step=resolution,
                                                statistic=stats.skew)
        fig, ax = plt.subplots()
        ax.plot(independant_out_array, ring_std, 'g',
                label='Ring Standard Deviation')
        ax.plot(independant_out_array, ring_std_no_mask, 'b',
                label='Ring Standard Deviation No old_files')
        legend = ax.legend(loc='upper right')
        plt.show()

        fig, ax = plt.subplots()
        ax.plot(independant_out_array, ring_s, 'g',
                label='Ring Skew')
        ax.plot(independant_out_array, ring_s_no_mask, 'b',
                label='Ring Skew No old_files')
        legend = ax.legend(loc='upper right')
        plt.show()

        # plt.plot(independant_out_array, test_uncertainty / integrated * 100)
        # plt.show()

        plt.plot(
            independant_out_array, integrated, 'g',
            independant_out_array, integrated - uncertainty_array, 'b',
            independant_out_array, integrated + uncertainty_array, 'r'
        )
        plt.show()

        # plt.plot(
        # independant_out_array, integrated, 'g',
        #     independant_out_array, integrated - test_uncertainty, 'b',
        #     independant_out_array, integrated + test_uncertainty, 'r'
        # )
        # plt.show()
    # Save file as chi file with uncertainties
    out_array = np.c_[independant_out_array, integrated, uncertainty_array]
    # print out_array.shape
    chi_filename = os.path.splitext(bulk_image_path_no_ext)[0] + out_ending
    out_file_name = IO.save_chi(out_array, a, chi_filename, mask_file, image_path)
    if generate_mask is True:
        return out_file_name, mask_file_name
    return out_file_name


def write_pdf(chi_file=None, background_file=None, background_level=None,
              pdf_format='QA', output='@r.@o', qmax='statistics',
              composition=None, qmaxinst=0.0,
              qmin=0.0, relative_max_uncertainty=.90,
              plot_verbose=False, bg_refine_qmin_index=None,
              bg_refine_qmax_index=None):
    """
    Generates the G(r) and F(q) files, potentially using optimized values
    for the background subtraction and qmax determination

    Parameters
    ----------
    chi_file: str
        Location of the chi file to use as the foreground to generate the PDF
    background_file: str
        Location of the chi file to use as the background
    background_level: float
        The background level, if not set automatically generate the background
        scale
    pdf_format: str
        The format of the chi file
    out_put: str
        The naming convention for the output PDF
    qmax: float or {'statistics', 'ripple'}
        If float, the qmax value to be used in the FFT.  If statistics,
        generate qmax by examining the statistical uncertainty of I(Q).  If
        ripple generate the PDFs at various qmax and compare the rippleness
        between the various qmax values
    composition: str
        The compostion of the material to be studied
    qmaxinst: float
        The maximum reliable qmax generated by the detector, if 0.0 the max
        of the chi file
    qmin: float
        The minimum q to be used in the FFT
    relative_max_uncertainty: float [0,1]
        If qmax is statistics the percentile of the uncertainty to except
    plot_verbose: bool
        If true, generate plots
    bg_refine_qmin_index: int
        The lower bound on the range of values to fit between the background
        and foreground for determining the background level.
    bg_refine_qmax_index: int
        The upper bound on the range of values to fit between the background
        and foreground for determining the background level.

    Returns
    -------
    str:
        The file name of the output
    float:
        The background level
    float:
        The qmax
    int:
        The background refine lower bound
    int:
        The background refine upper bound
    """
    #Get composition
    if composition is None:
        composition = pdf_parameters.generate_composition()
    # load up sample I[Q]
    chi, chi_file = IO.load_chi_file(chi_file)
    if os.path.split(chi_file)[0] != '':
        os.chdir(os.path.split(chi_file)[0])
    sample_q, sample, uncertainty = chi[:, 0], chi[:, 1], chi[:, 2]

    if qmaxinst > np.amax(sample_q):
        qmaxinst = np.amax(sample_q)
    elif qmaxinst == 0:
        qmaxinst = np.amax(sample_q)
    # perform background subtraction
    if background_file is not '' and background_level is None:
        # implies that you want to load a background file
        background_level, background_q, background, background_file,  \
        bg_refine_qmin_index, bg_refine_qmax_index = \
            pdf_parameters.background_subtract(sample, background_file,
                                               background_level,
                                               bg_refine_qmin_index=bg_refine_qmin_index,
                                               bg_refine_qmax_index=bg_refine_qmax_index)

        if plot_verbose is True:
            fig, ax = plt.subplots()
            ax.plot(sample_q, sample, 'g', label='Sample')
            ax.legend(loc='upper right')
            plt.title('Sample I(Q)')
            plt.xlabel('Q (1/A)')
            plt.ylabel('Raw Counts')
            plt.savefig('/home/christopher/sample.png',
                        bbox_inches='tight',
                        transparent=True)
            plt.show()

            fig, ax = plt.subplots()
            ax.plot(sample_q, background, 'b',
                    label='Background')
            ax.legend(loc='upper right')
            plt.title('Background (IQ)')
            plt.xlabel('Q (1/A)')
            plt.ylabel('Raw Counts')
            plt.savefig('/home/christopher/background.png',
                        bbox_inches='tight',
                    transparent=True)
            plt.show()

            fig, ax = plt.subplots()
            ax.plot(sample_q, sample, 'g', label='Sample')
            ax.plot(sample_q, background * background_level, 'b',
                    label='Scaled Background')
            ax.plot(sample_q, sample - background * background_level, 'k',
                    label='Sample-Background')
            ax.legend(loc='upper right')
            plt.title('Background Subtraction')
            plt.xlabel('Q (1/A)')
            plt.ylabel('Raw Counts')
            plt.savefig('/home/christopher/integrated_minus_background.png',
                        bbox_inches='tight',
                    transparent=True)
            plt.show()
    # Determine Qmax
    if qmax is 'statistics':
        qmax = pdf_parameters.qmax_statistics(sample_q, sample, uncertainty,
                                        relative_max_uncertainty)
    elif qmax is 'ripple':
        qmax = pdf_parameters.ripple_minima(chi_file, background_file,
                                            background_level, composition,
                                            qmin, qmaxinst)
    if plot_verbose is True:
        plt.plot(sample_q, uncertainty/sample * 100)
        plt.show()
    # PDFGETX3 call/ Generate PDF
    gr_file, background_level, qmax = Calculate.pdfgetx3(chi_file,
                                                     background_file, background_level,
                       composition, qmax, qmin, qmaxinst, pdf_format, output)
    return gr_file, background_level, qmax, bg_refine_qmin_index, bg_refine_qmax_index

    # MODELS


if __name__ == "__main__":
    plt.ioff()
    chi_name = tif_to_iq(
        poni_file=None,
        image_path=None,
        mask_file=None,
        plot_verbose=True,
    )

    print 'Start PDF'
    gr_file, bg, calculated_qmax, bg_refine_qmin_index, bg_refine_qmax_index = \
        write_pdf(
            chi_file=None,
            background_file=None,
            composition='Ni.3Pd.7',
            plot_verbose=True,
            qmax='ripple',
            qmaxinst=27.1,
            bg_refine_qmin_index=None,
            bg_refine_qmax_index=None)

    radius, gr = IO.load_gr_file(gr_file)
    plt.plot(radius, gr)
    plt.show()