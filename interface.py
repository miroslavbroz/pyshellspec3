#!/usr/bin/env python3

import copy
import os
import sys
import nlopt
import numpy as np
import shutil
import subprocess

from astropy import units
from matplotlib import pyplot as plt
from multiprocessing import cpu_count
from multiprocessing import Queue
from multiprocessing import Process
from scipy import ndimage
from scipy.interpolate import interpn
from scipy.interpolate import RectBivariateSpline
from scipy.io import FortranFile
from .auxilliary import DFT2d
from .auxilliary import get_offset
from .auxilliary import get_mulfac
from .auxilliary import interpolate1d
from .auxilliary import interpolate1d_hermite
from .auxilliary import quad_phase
from .auxilliary import write_file
from .auxilliary import append_file
from .auxilliary import rotateX
from .auxilliary import rotateZ
from .definitions import fitter_definitions
from .definitions import filters, eff_wave, eff_band, calibration_flux
from .definitions import shsp_template   as default_template
from .definitions import shsp_executable as default_executable
from .definitions import shsp_params     as default_params
from .definitions import shsp_abundance  as default_abundance
from .model import Model
from .objects import CentralObject
from .objects import Companion
from .objects import Disk
from .objects import Nebula
from .objects import Distance
from .objects import Envelope
from .objects import Orbit
from .objects import Spot
from .objects import Ufo
from .objects import Flow
from .objects import Jet
from .observations import Data
from .plotting import plot_light_curve
from .plotting import plot_squared_visibility
from .plotting import plot_triple_product

from . import limcof

class Interface(object):
    """
    This class should serve as an interface to the 
    Pyshellspec
    """
    def __init__(self, data=None, model=None, image_size=None, shellspec_abundance=None,
                 shellspec_executable=None, shellspec_params=None, 
                 shellspec_template=None, debug=False, ncpu=1,
                 if_phase_precision=3, df_phase_precision=3, lc_phase_precision=2, sed_phase_precision=2, spe_phase_precision=2,
                 if_ew_precision=9, df_ew_precision=10, lc_ew_precision=9, sed_ew_precision=10, spe_ew_precision=10,
                 exclude_cp=False, exclude_t3amp=False, exclude_vis2=False, exclude_visamp=False, exclude_visphi=False,
                 use_offset=True, use_differential=False, dry_run=False, overwrite=True):
        """
        Constructs the class.
        :param data: All data files wrapped in pyshellspec.Data type.
        :param model: All pyshellspec models wrapped in pyshellspec.Model type.
        :param image_size: Size of the zero-padded image. This parameter is important only for FFT.
        :param shellspec_abundance: Path to shellspec 'abundances' file, if None, built-in is used
        :param shellspec_executable:  Path to shellspec executable file, if None, built-in is used
        :param shellspec_params:  Path to shellspec parameter (param.inc) file, if None, built-in is used
        :param shellspec_template: Path to shellspec control (shellspec.in) file, if None, built-in is used
        :param if_phase_precision: The phase rounding order for interferometry,
                                   phase = np.around(phase, if_phase_precision)
        :param df_phase_precision: The same for differential interferometry.
        :param lc_phase_precision: The same for photometric data.
        :param sed_phase_precision: The same for SED data.
        :param spe_phase_precision: The same for spectral data.
        :param if_ew_precision: The effective wavelength rounding order for interferometry,
                                ew = np.around(ew, if_ew_precision)
        :param df_ew_precision: The same for differential interferometry.
        :param lc_ew_precision: The same for photometric data.
        :param sed_ew_precision: The same for SED data.
        :param spe_ew_precision: The same for spectral data.
        :param debug: If true, more information is printed / plotted.
        :param ncpu: number of cpus that will be employd in computation
        :param exclude_cp: excludes closure phase from chi^2 fitting
        :param exclude_t3amp: excludes triple product amplitude from chi^2 fitting
        :param exclude_vis2: excludes squared visibility from chi^2 fitting
        :param exclude_visamp: excludes visibility amplitude from chi^2 fitting
        :param exclude_visphi: excludes visibility phase from chi^2 fitting
        :param use_offset: offset synthetic lightcurves to better match the observed ones
        :param use_differential: assume differential visibilities and shift |V|
        :param dry_run: do NOT call shellspec, read previous output files; very useful for debugging
        :param overwrite: do NOT overwrite files; very useful for restarts of hi-res simulations
        :return:
        """
    
        # stores data
        self.__data = data
    
        # stores model
        self.__model = model

        # stores the size of images from which FTis computed
        self.__image_size = image_size

        # sets which interferometric observables should not be fitted
        self.__exclude_cp = exclude_cp
        self.__exclude_t3amp = exclude_t3amp
        self.__exclude_vis2 = exclude_vis2
        self.__exclude_visamp = exclude_visamp
        self.__exclude_visphi = exclude_visphi
        self.__use_offset = use_offset
        self.__use_differential = use_differential

        # sets address of all shellspec files
        # abundance
        if shellspec_abundance is None:
            self.__shellspec_abundance_file = default_abundance
        else:
            self.__shellspec_abundance_file = shellspec_abundance
        
        # executable
        if shellspec_executable is None:
            self.__shellspec_executable_file = default_executable
        else:
            self.__shellspec_executable_file = shellspec_executable
        
        # grid parameters
        if shellspec_params is None:
            self.__shellspec_params_file = default_params
        else:
            self.__shellspec_params_file = shellspec_params
            
        # template
        if shellspec_template is None:
            self.__shellspec_template_file = default_template
        else:
            self.__shellspec_template_file = shellspec_template

        # stores precision in phase and effective wavelength for
        # photometric and interferometric data in 10 ** prec
        self.__if_phase_precision = if_phase_precision
        self.__df_phase_precision = df_phase_precision
        self.__lc_phase_precision = lc_phase_precision
        self.__sed_phase_precision = sed_phase_precision
        self.__spe_phase_precision = spe_phase_precision
        self.__if_ew_precision = if_ew_precision
        self.__df_ew_precision = df_ew_precision
        self.__lc_ew_precision = lc_ew_precision
        self.__sed_ew_precision = sed_ew_precision
        self.__spe_ew_precision = spe_ew_precision

        # read the template file
        self.__read_template()

        # update the model from  the template file -- if some was passed
        # and model was passed
        if shellspec_template is not None and model is not None:
            self.__update_model_from_initial_template()
        
        # update the model from  the template file -- if some was passed
        # but no model was passed
        elif shellspec_template is not None and model is None:
            
            # creates the model from the shellspec 
            # template
            objects = []
            if int(self.__read_parameter_from_template('istar')) > 0:
                objects.append(CentralObject())
            if int(self.__read_parameter_from_template('icomp')) > 0:
                objects.append(Companion())
            if int(self.__read_parameter_from_template('idisc')) > 0:
                objects.append(Disk())
            if int(self.__read_parameter_from_template('inebl')) > 0:
                objects.append(Nebula())
            if int(self.__read_parameter_from_template('ienv')) > 0:
                objects.append(Envelope())
            if int(self.__read_parameter_from_template('ispot')) > 0:
                objects.append(Spot())
            if int(self.__read_parameter_from_template('iufo')) > 0:
                objects.append(Ufo())
            if int(self.__read_parameter_from_template('iflow')) > 0:
                objects.append(Flow())
            if int(self.__read_parameter_from_template('ijet')) > 0:
                objects.append(Jet())
                
            # always append orbit -- this is necessary now
            # because pyshellspec needs that always
            objects.append(Orbit())
            
            # construct the model
            model = Model(objects=objects)
            
            # assign the model to the object
            self.__model = model
            
            # update the model from template
            self.__update_model_from_initial_template()
        
        # extract dependent and independent variables from data
        # and create arrays for consequent comparison
        if self.__data is not None:
            self.__cp_comparison = self.__ready_comparison('cp')
            self.__lc_comparison, self.__lcsyn = self.__ready_comparison('lc')
            self.__vis2_comparison = self.__ready_comparison('vis2')
            self.__vis_comparison = self.__ready_comparison('vis')
            self.__sed_comparison, self.__sedsyn = self.__ready_comparison('sed')
            self.__spe_comparison, self.__spesyn = self.__ready_comparison('spe')

        # attributes for effective wavelength and phases 
        # at which the model will be evaluated
        self.__has_phase = False

        # number of employed cpu
        # correct for situation where we exceed the number of
        # cpus in the computer
        ncpu = min([cpu_count(), ncpu])
        self.__ncpu = ncpu
        print("ncpu = ", ncpu)

        # set debug mode
        self.debug = debug
        self.dry_run = dry_run
        self.overwrite = overwrite

        # empty array for iteration
        self.__iterations = []

        self.limcofst = None
        self.limcofcp = None

    def compute_chi2(self, pars=[], fitparams=None, verbose=False):
        """
        Comutes chi^2
        :param pars:
        :param fitparams:
        :param verbose:
        :return:
        """
        # just in case the function is called withou
        # parameters
        if fitparams is None:
            fitparams = self.get_fitted_parameters()

        # update parameters if the passed array
        # is the same length as that of fitted parameters
        # and then update parameters based on the orbital
        # model
        if len(pars) == len(fitparams):
            self.__update_fitted_parameters(pars, fitparams)
        self.__update_from_orbit()  # modified by MB, Nov 27th 2017

        # propagate them into the shellspec control file
        self.set_model_to_shellspec()

        # prepare limb-darkening coef. for both (central) star and companion
        self.limcofst= limcof.Limcof(Teff=self.__model['central_object']['tstar'].to('Kelvin').value, R=self.__model['central_object']['rstar'].to('solRad').value, M=self.__model['central_object']['emstar'].to('solMass').value)
        self.limcofcp= limcof.Limcof(Teff=self.__model['companion']['tempcp'].to('Kelvin').value, R=self.__model['companion']['rcp'].to('solRad').value, M=self.__model['companion']['qq'].value*self.__model['central_object']['emstar'].to('solMass').value)

        # compute light curves
        self.compute_lc()

        # compute interferometric observables
        self.compute_if_DFT(dtype='if')
        self.compute_if_DFT(dtype='df')

        # compute SED
        self.compute_sed()

        # compute spectra
        self.compute_spe()

        # compute partial chi-squares
        # light curves
        self.__lc_comparison['chi2'] = ((self.__lc_comparison['magnitude'] - self.__lc_comparison['magnitudesyn']) / self.__lc_comparison['error']) ** 2
        chi2_lc = np.sum(self.__lc_comparison['chi2'])
        n_lc = self.__lc_comparison['magnitude'].size

        # squared visibilities
        chi2_vis2 = 0.0
        n_vis2 = 0
        if not self.__exclude_vis2:
            for i in range(len(self.__vis2_comparison['vis2data'])):
                if self.__vis2_comparison['vis2data'][i] == None:
                    self.__vis2_comparison['chi2'].append(None)
                else:
                    weight = self.__vis2_comparison['weight'][i]
                    self.__vis2_comparison['chi2'][i] = weight * ((self.__vis2_comparison['vis2data'][i] - self.__vis2_comparison['vis2syn'][i]) / self.__vis2_comparison['vis2err'][i]) ** 2
                    chi2_vis2 += self.__vis2_comparison['chi2'][i]
                    n_vis2 += weight

        # closure phases
        chi2_cp = 0.0
        n_cp = 0
        if not self.__exclude_cp:
            for i in range(len(self.__cp_comparison['t3phi'])):
                if self.__cp_comparison['t3phi'][i] == None:
                    self.__cp_comparison['chi2phi'].append(None)
                else:
                    self.__cp_comparison['chi2phi'][i] = ((self.__cp_comparison['t3phi'][i] - self.__cp_comparison['t3phisyn'][i]) / self.__cp_comparison['t3phierr'][i]) ** 2
                    chi2_cp += self.__cp_comparison['chi2phi'][i]
                    n_cp += 1

        # triple product amplitudes
        chi2_t3amp = 0.0
        n_t3amp = 0
        if not self.__exclude_t3amp:
            for i in range(len(self.__cp_comparison['t3amp'])):
                if self.__cp_comparison['t3amp'][i] == None:
                    self.__cp_comparison['chi2amp'][i] = None
                else:
                    weight = self.__cp_comparison['weight'][i]
                    self.__cp_comparison['chi2amp'][i] = weight * ((self.__cp_comparison['t3amp'][i] - self.__cp_comparison['t3ampsyn'][i]) / self.__cp_comparison['t3amperr'][i]) ** 2
                    chi2_t3amp += self.__cp_comparison['chi2amp'][i]
                    n_t3amp += weight

        # visibility amplitudes
        chi2_visamp = 0.0
        n_visamp = 0
        if not self.__exclude_visamp:
            for i in range(len(self.__vis_comparison['visamp'])):
                if self.__vis_comparison['visamp'][i] == None:
                    self.__vis_comparison['chi2amp'].append(None)
                else:
                    weight = self.__vis_comparison['weight'][i]
                    self.__vis_comparison['chi2amp'][i] = weight * ((self.__vis_comparison['visamp'][i] - self.__vis_comparison['visampsyn'][i]) / self.__vis_comparison['visamperr'][i]) ** 2
                    chi2_visamp += self.__vis_comparison['chi2amp'][i]
                    n_visamp += weight

        # visibility phases
        chi2_visphi = 0.0
        n_visphi = 0
        if not self.__exclude_visphi:
            for i in range(len(self.__vis_comparison['visphi'])):
                if self.__vis_comparison['visphi'][i] == None:
                    self.__vis_comparison['chi2phi'].append(None)
                else:
                    weight = self.__vis_comparison['weight'][i]
                    self.__vis_comparison['chi2phi'][i] = weight * ((self.__vis_comparison['visphi'][i] - self.__vis_comparison['visphisyn'][i]) / self.__vis_comparison['visphierr'][i]) ** 2
                    chi2_visphi += self.__vis_comparison['chi2phi'][i]
                    n_visphi += weight

        # spectral-energy distribution
        self.__sed_comparison['chi2'] = ((self.__sed_comparison['flux'] - self.__sed_comparison['fluxsyn']) / self.__sed_comparison['error']) ** 2
        chi2_sed = np.sum(self.__sed_comparison['chi2'])
        n_sed = self.__sed_comparison['flux'].size

        # spectra
        self.__spe_comparison['chi2'] = ((self.__spe_comparison['flux'] - self.__spe_comparison['fluxsyn']) / self.__spe_comparison['error']) ** 2
        chi2_spe = np.sum(self.__spe_comparison['chi2'])
        n_spe = self.__spe_comparison['flux'].size

        # compute total chi-square
        n = n_lc + n_vis2 + n_cp + n_t3amp + n_visamp + n_visphi + n_sed + n_spe
        chi2 = chi2_lc + chi2_vis2 + chi2_cp + chi2_t3amp + chi2_visamp + chi2_visphi + chi2_sed + chi2_spe

        ns = [n_lc, n_vis2, n_cp, n_t3amp, n_visamp, n_visphi, n_sed, n_spe, n]
        chi2s = [chi2_lc, chi2_vis2, chi2_cp, chi2_t3amp, chi2_visamp, chi2_visphi, chi2_sed, chi2_spe, chi2]
        r_chi2s = [chi2 / (n+1.0e-8) for chi2, n in zip(chi2s, ns)]

        # append info on each iteration
        if len(pars) == 0:
            one_iter = [rec['value'] for rec in fitparams]
        else:
            one_iter = [par for par in pars]
        for rec in ns:
            one_iter.append(rec)
        for rec in chi2s:
            one_iter.append(rec)
        for rec in r_chi2s:
            one_iter.append(rec)
        self.__iterations.append(one_iter)

        self.write_iterations()
        if len(self.__iterations) % 10 == 0:
            self.write_iterations('fit.intermediate.log')

        self.write_model()  # dbg

        # return the output
        if verbose:
            return chi2s.extend(r_chi2s)
        else:
            return chi2

    def compute_if_DFT(self,dtype='if'):
        """
        Computes synthetic visibilities.
        """
        if dtype == 'df':
            phase_precision = self.__df_phase_precision
            ew_precision = self.__df_ew_precision
        else:
            phase_precision = self.__if_phase_precision
            ew_precision = self.__if_ew_precision

        # first copy the model to shellspec
        self.set_model_to_shellspec()

        # get phases and ews at which we shall compute
        ifsyn = self.get_ew_and_phase(dtype)

        # get the position of the barycentre, inclination
        # ascending node longitude
        if self.__model.has_object('orbit'):
            bar_pos = self.__model['orbit'].get_barycentre()
            omega = self.__model['orbit']['omega_an']
            incl = self.__model['orbit']['dinc']
        else:
            raise TypeError('The tool is currently working with binaries only '
                            'i.e. for models, where orbit has been defined.')

        # get the physical resolution along the
        # x-axis and size of the grid in pixels
        if self.__model.has_object('grid'):
            npx, npy, npz = self.__model['grid'].get_pixel_size()
        else:
            raise KeyError('Object grid was not found among defined objects.'
                           'It should be alway defined')

        # get the size of the final image - it is padded to get a
        # square format, which is convenient for rotation of the
        # image. Nonetheless padding of the image is not
        # necessary for DFT.
        newsize = self.__image_size
        if newsize is None:
            if npx > npy:
                newsize = npx
            else:
                newsize = npy

        # extract the resolution along x-axis
        phys_res = self.__model['grid']['stepx']

        # compute spatial frequencies
        # get parallax and stepsize
        plx = self.__model['distance']['dd'].to('m')
        step_physical = phys_res.to('m')

        # compute angular step
        step_angular = step_physical.value / plx.value

        # axis of the image
        if newsize % 2 > 1:
            xscale = np.arange(-newsize // 2 + 1, newsize // 2 + 1, 1) * step_angular
        else:
            xscale = np.arange(-newsize // 2, newsize // 2, 1) * step_angular
        yscale = xscale.copy()

        if self.__ncpu > 1:

            # get number of cpu
            ncpu = self.__ncpu

            # get the number of wavelengths
            newave = len(ifsyn['eff_wave'])
            print("newave = " + str(newave))

            # creates output files (different directory for each wavelength)
            directories = []
            for i in range(newave):

                # get directory
                directory = os.path.join(os.getcwd(), 'temp' + dtype + str(i).zfill(2))
                if not os.path.isdir(directory):
                    os.mkdir(directory)
                directories.append(directory)

            # go over each wavelength
            for i in range(0, newave, ncpu):

                # empty list for threads
                threads = []
                queues = []
                results = []

                # spawn processes
                nprocess = min([ncpu, len(ifsyn['eff_wave']) - i])
                for j in range(nprocess):

                    # get phase and wavelength
                    ew = ifsyn['eff_wave'][i + j] * units.m
                    phase = ifsyn['phase'][i + j]

                    # spawn queue
                    queues.append(Queue())

                    # define process args
                    args = (directories[i+j], ew, phase, npx, npy, phys_res, bar_pos, incl, omega, newsize, queues[j], True)

                    # spawn processes
                    threads.append(Process(target=self.__pp_get_images_one_wavelength, args=args))

                # run processes
                for j in range(len(threads)):
                    threads[j].start()

                # readout processes and join threads
                for j in range(len(threads)):
                    results.append(queues[j].get())
                    queues[j].close()
                    threads[j].join()
                    if threads[j].exitcode != 0:
                        print("compute_if_DFT: Error running shellspec! exitcode = " + str(threads[j].exitcode))
                        sys.exit(1)

                # compute DFT over the data
                for res in results:
                    ew, phases, images = res
                    for phase, img in zip(phases, images):
                        self.__get_if_observables_DFT(img, ew.value, phase, xscale, yscale, ew_precision, phase_precision)
        else:  # 1 cpu

            # get the number of wavelengths
            newave = len(ifsyn['eff_wave'])
            print("newave = " + str(newave))

            # creates output files (different directory for each wavelength)
            directories = []
            for i in range(newave):

                # get directory
                directory = os.path.join(os.getcwd(), 'temp' + dtype + str(i).zfill(2))
                if not os.path.isdir(directory):
                    os.mkdir(directory)
                directories.append(directory)

            # go over each wavelength
            for i in range(newave):
                print("directory = ", directories[i])

                # set directory and phase
                ew = ifsyn['eff_wave'][i] * units.m
                phase = ifsyn['phase'][i]

                # process one wavelength
                self.__process_wavelength_DFT(directories[i], ew, phase, npx, npy, phys_res, bar_pos,
                                              incl, omega, newsize, xscale, yscale, image_only=True)

        if dtype == 'df':
            self.zero_slips()
            self.adjust_differential()
            self.adjust_slips()
            self.adjust_differential()
            self.adjust_slips()


    def zero_slips(self):
        """
        Set mulfac to one, offset and slips to zero.
        """
        n = len(self.__vis_comparison['mulfac'])
        self.__vis_comparison['mulfac'][:] = np.ones(n)
        self.__vis_comparison['offset'][:] = np.zeros(n)
        self.__vis_comparison['slips'][:] = np.zeros(n)
        self.__vis_comparison['visamp'][:] = self.__vis_comparison['visamp_']
        self.__vis_comparison['visphi'][:] = self.__vis_comparison['visphi_']


    def adjust_differential(self):
        """
        Multiply synthetic (differential) visibility amplitude |V|,
        and shift phase arg V to match the observed values.
        """
        
        filenames = self.__vis_comparison['filename']
        for filename in np.unique(filenames):
            find = np.where(filenames == filename)[0]
            visamp = self.__vis_comparison['visamp'][find]
            visampsyn = self.__vis_comparison['visampsyn'][find]

            if self.__use_differential:
                mulfac = get_mulfac(visamp, visampsyn)
                mulfac = mulfac.x
            else:
                mulfac = 1.0
            self.__vis_comparison['visampsyn'][find] = visampsyn * mulfac
            self.__vis_comparison['mulfac'][find] *= np.ones(len(visampsyn)) * mulfac

            visphi = self.__vis_comparison['visphi'][find]
            visphisyn = self.__vis_comparison['visphisyn'][find]

            if self.__use_differential:
                offset = get_offset(visphi, visphisyn)
                offset = offset.x
            else:
                offset = 0.0
            self.__vis_comparison['visphisyn'][find] = visphisyn + offset
            self.__vis_comparison['offset'][find] += np.ones(len(visphisyn)) * offset


    def adjust_slips(self):
        """
        Adjust slips of the differential phase by +-360 deg.
        """

        visphi = self.__vis_comparison['visphi']
        visphisyn = self.__vis_comparison['visphisyn']

        find = np.where(visphisyn-visphi > 180.0)
        slips = -360.0
        self.__vis_comparison['visphisyn'][find] = visphisyn[find] + slips
        self.__vis_comparison['slips'][find] = np.ones(len(find))*slips

        find = np.where(visphisyn-visphi < -180.0)
        slips = 360.0
        self.__vis_comparison['visphisyn'][find] = visphisyn[find] + slips
        self.__vis_comparison['slips'][find] = np.ones(len(find))*slips


    def compute_if_FFT(self):
        """
        Computes synthetic visibilities.
        """
        # first copy the model to shellspec
        self.set_model_to_shellspec()

        # get phases and ews at which we shall compute
        ifsyn = self.get_ew_and_phase('if')

        # get the position of the barycentre, inclination
        # ascending node longitude
        if self.__model.has_object('orbit'):
            bar_pos = self.__model['orbit'].get_barycentre()
            omega = self.__model['orbit']['omega_an']
            incl = self.__model['orbit']['dinc']
        else:
            raise TypeError('The tool is currently working with binaries only '
                            'i.e. for models, where orbit has been defined.')

        # get the physical resolution along the
        # x-axis and size of the grid in pixels
        if self.__model.has_object('grid'):
            npx, npy, npz = self.__model['grid'].get_pixel_size()
        else:
            raise KeyError('Object grid was not found among defined objects.'
                           'It should be alway defined')

        # get the size of the padded image
        newsize = self.__image_size

        # extract the resolution along x-axis
        phys_res = self.__model['grid']['stepx']

        # compute spatial frequencies
        # get parallax and stepsize
        plx = self.__model['distance']['dd'].to('m')
        step_physical = phys_res.to('m')

        # compute angular step
        step_angular = step_physical.value / plx.value

        # compute the frequency
        fu_img = np.fft.fftshift(np.fft.fftfreq(newsize, d=step_angular))
        fv_img = fu_img.copy()

        # set directory
        directory = os.path.join(os.getcwd(), 'tempfft')
        if not os.path.isdir(directory):
            os.mkdir(directory)

        for i in range(len(ifsyn['eff_wave'])):

            # readout one effective wavelength and corresponding phases
            ew = ifsyn['eff_wave'][i] * units.m
            phase = ifsyn['phase'][i]

            # set the effective wavelength to template
            self.set_wavelength(w0=ew, wn=ew + 1.0 * units.nm, step=1.0 * units.nm)

            # write template and phases
            phase_file = os.path.join(directory, 'phases')
            control_file = os.path.join(directory, 'shellspec.in')
            self.write_phase(phase, filename=phase_file)
            self.write_template(filename=control_file)

            # run shellspec
            cwd = os.getcwd()
            os.chdir(directory)
            self.__run_shellspec()
            os.chdir(cwd)

            # read the images - one for each phase
            for j in range(0, len(phase)):

                # set the filename
                filename = os.path.join(directory, '2Dimage_%03d' % (j + 1))

                # read the image
                img = self.__read_image(filename)

                # compute its Fourier transform
                fftimg = self.__FT_one_image(img, npx, npy, phys_res, bar_pos, phase[j], incl, omega, newsize)

                # in case of debug mode, save both - the image and its Fourier transform
                if self.debug:
                    # the image
                    tmp = '%.2f' % (ew.to('Angstrom').value)
                    figname = '.'.join(['fftimg', tmp, str(phase[j]), 'png'])
                    plt.imshow(self.debug_image, cmap='plasma')
                    plt.xlabel('x [pxl]')
                    plt.ylabel('y [pxl]')
                    plt.colorbar()
                    plt.tight_layout()
                    plt.savefig(figname)
                    plt.close()

                    # and its Fourier transform
                    figname = '.'.join(['fftabs', tmp, str(phase[j]), 'png'])
                    plt.imshow(np.abs(fftimg), cmap='plasma')
                    plt.xlabel('u [pxl]')
                    plt.ylabel('v [pxl]')
                    plt.colorbar()
                    plt.tight_layout()
                    plt.savefig(figname)
                    plt.close()

                    figname = '.'.join(['fftang', tmp, str(phase[j]), 'png'])
                    plt.imshow(np.angle(fftimg, deg=True), cmap='plasma')
                    plt.xlabel('u [pxl]')
                    plt.ylabel('v [pxl]')
                    plt.colorbar()
                    plt.tight_layout()
                    plt.savefig(figname)
                    plt.close()

                # transpose AFTER debugging and BEFORE observables! otherwise, it's misleading...
                fftimg = fftimg.T

                # compute observables
                self.__get_if_observables_FFT(fftimg, ew.value, phase[j], fu_img, fv_img, self.__if_ew_precision, self.__if_phase_precision)

        return ifsyn

    def compute_lc(self):
        """
        Computes synthetic light curve.
        """
        # first copy the model to shellspec
        self.set_model_to_shellspec()

        # free synthetic light curves from previous
        # computation
        self.__lcsyn['magnitude'] = []

        # get ews and phases at which we shall compute
        lcsyn = self.get_ew_and_phase('lc')

        # get number of cpu
        ncpu = self.__ncpu

        # get the number of wavelengths
        newave = len(lcsyn['eff_wave'])
        print("newave = " + str(newave))

        # creates output files
        directories = []
        for i in range(newave):

            # get directory
            directory = os.path.join(os.getcwd(), 'templc' + str(i).zfill(2))
            if not os.path.isdir(directory):
                os.mkdir(directory)
            directories.append(directory)

        # go over each wavelength
        for i in range(0, newave, ncpu):

            # empty list for threads
            threads = []
            queues = []
            results = []

            # spawn processes
            nprocess = min([ncpu, newave - i])
            print("nprocess = " + str(nprocess))

            # effective wavelength and phase for one passband
            for j in range(nprocess):
                ew = lcsyn['eff_wave'][i + j]
                phase = lcsyn['phase'][i + j]
                # print "ew = " + str(ew)  # dbg
                # print "phase = " + str(phase)  # dbg

                # spawn queue
                queues.append(Queue())

                # define process args
                args = (directories[i+j], ew, phase, queues[j])

                # spawn processes
                threads.append(Process(target=self.__pp_get_lc_one_wavelength, args=args))

            # run processes
            for j in range(nprocess):
                threads[j].start()

            # readout processes and  join threads
            for j in range(nprocess):
                results.append(queues[j].get())
                queues[j].close()
                threads[j].join()
                if threads[j].exitcode != 0:
                    print("compute_lc: Error running shellspec! exitcode = " + str(threads[j].exitcode))
                    sys.exit(1)

            # append to the output structure
            for j in range(nprocess):
                mags = results[j]
                lcsyn['magnitude'].append(mags)

        # compute observables
        self.__get_lc_observables(lcsyn)

        # overwrite lcsyn
        self.__lcsyn = copy.deepcopy(lcsyn)

    def get_ew_and_phase(self, dtype):
        """
        Determines at which phase and effective wavelength will the model be evaluated
        :param dtype:
        :return:
        """

        # compute phases
        if not self.__has_phase:
            self.__set_phase()

        # photometric data
        if dtype == 'lc':

            # copy the object with effective wavelengths and passbands
            lc_phase_precision = self.__lc_phase_precision
            lc_ew_precision = self.__lc_ew_precision
            ew_ph = copy.deepcopy(self.__lcsyn)

            # add list for phases
            ew_ph['phase'] = []

            # determine phases for each passband
            for i, pband in enumerate(ew_ph['passband']):

                # extract all phases with the coresponding passband
                # ind = np.where(self.__lc_comparison['passband'] == pband)[0]

                # print self.__lc_comparison['phase']
                # phase = self.__lc_comparison['phase'][ind]

                # get extrema
                step = 10 ** -lc_phase_precision
                # pmin = max([np.around(phase.min() - 2 * step, lc_phase_precision), 0.0])
                # pmax = min([np.around(phase.max() + 2 * step, lc_phase_precision), 1.0 - 10 ** -lc_phase_precision])
                pmin = 0.0
                pmax = 1.0 - step

                # create array of phases
                compute_phase = np.arange(pmin, pmax + step / 2., step)

                # apend it to the structure
                ew_ph['phase'].append(compute_phase)

            return ew_ph

        # visibilities and closure phases
        elif dtype == 'if':
            # first join visibility and closure phase data
            if_phase_precision = self.__if_phase_precision
            if_ew_precision = self.__if_ew_precision
            ew = self.__vis2_comparison['eff_wave'].copy()
            ew = np.append(ew, self.__cp_comparison['eff_wave'])
#            ew = np.append(ew, self.__vis_comparison['eff_wave'])
            phase = self.__vis2_comparison['phase'].copy()
            phase = np.append(phase, self.__cp_comparison['phase'])
#            phase = np.append(phase, self.__vis_comparison['phase'])

            # around phases
            phase = np.around(phase, if_phase_precision)

            # create a dictionary of ews and phases
            ew_ph = dict(eff_wave=np.unique(np.around(ew, if_ew_precision)).tolist(), phase=[])

            # attach the corresponding phases
            for eff_wave in ew_ph['eff_wave']:
                diff = np.absolute(ew - eff_wave)
                ind = np.where(diff < 10 ** -if_ew_precision)[0]
                ew_ph['phase'].append(phase[ind].tolist())

            # make all phases unique
            for i in range(0, len(ew_ph['eff_wave'])):
                ew_ph['phase'][i] = np.unique(ew_ph['phase'][i]).tolist()

            return ew_ph

        # differential visibilities (separate)
        elif dtype == 'df':
            phase_precision = self.__df_phase_precision
            ew_precision = self.__df_ew_precision
            ew = self.__vis_comparison['eff_wave'].copy()
            phase = self.__vis_comparison['phase'].copy()

            phase = np.around(phase, phase_precision)
            ew_ph = dict(eff_wave=np.unique(np.around(ew, ew_precision)).tolist(), phase=[])

            for eff_wave in ew_ph['eff_wave']:
                diff = np.absolute(ew - eff_wave)
                ind = np.where(diff < 10 ** -ew_precision)[0]
                ew_ph['phase'].append(phase[ind].tolist())

            for i in range(0, len(ew_ph['eff_wave'])):
                ew_ph['phase'][i] = np.unique(ew_ph['phase'][i]).tolist()

            return ew_ph

        elif dtype == 'sed' or dtype == 'spe':
            if dtype == 'sed':
                phase_precision = self.__sed_phase_precision
                ew_precision = self.__sed_ew_precision
                ew = self.__sed_comparison['eff_wave'].copy()
                ph = self.__sed_comparison['phase'].copy()
                ew_ph = copy.deepcopy(self.__sedsyn)
            elif dtype == 'spe':
                phase_precision = self.__spe_phase_precision
                ew_precision = self.__spe_ew_precision
                ew = self.__spe_comparison['eff_wave'].copy()
                ph = self.__spe_comparison['phase'].copy()
                ew_ph = copy.deepcopy(self.__spesyn)

            ph_ = np.unique(np.around(ph, phase_precision)).tolist()
            ew = np.floor(ew/10**(-ew_precision))*10**(-ew_precision)  # np.around() w. float
            ew_ph['phase'] = ph_

            # attach the corresponding wavelengths
            for phase in ph_:
                diff = np.absolute(ph - phase)
                idx = np.where(diff < 10**(-phase_precision))[0]
                ew_ph['eff_wave'].append(ew[idx].tolist())

            for i in range(0,len(ph_)):
                ew_ph['eff_wave'][i] = np.unique(ew_ph['eff_wave'][i]).tolist()

            return ew_ph

    def get_fitted_parameters(self, attr=None):
        """
        Returns a list of all fitted parameters.
        :param attr:
        :return:
        """
        # empty array to store output
        fitparams = []

        # go over each object and parameter
        for objname in list(self.__model.keys()):
            for parname in list(self.__model[objname].keys()):
                if self.__model[objname].get_parameter(parname, 'fitted'):
                    # append the whole parameter
                    if attr is None:
                        value = self.__model[objname].get_parameter(parname, 'value')
                        vmin = self.__model[objname].get_parameter(parname, 'vmin')
                        vmax = self.__model[objname].get_parameter(parname, 'vmax')
                        fitparams.append(dict(object=objname,
                                              parameter=parname,
                                              value=value.value,
                                              vmin=vmin.value,
                                              vmax=vmax.value))
                    # append attribute of a parameter
                    else:
                        fitparams.append(self.__model[objname].get_parameter(parname, attr))
        return fitparams

    def get_image(self, phase, eff_wave, figname=None, img_only=False, **img_kwargs):
        """
        Creates an image of the system.
        :param phase:
        :param eff_wave: effective wavelength in m
        :param figname:
        :param img_only:
        :param img_kwargs:
        :return:
        """
        # convert to units
        ew = eff_wave * units.m

        # set the name
        if figname is None:
            figname = '.'.join(['image', str(eff_wave), str(phase), 'png'])

        # first copy the model to shellspec
        self.set_model_to_shellspec()

        # get the position of the barycentre, inclination
        # ascending node longitude
        if self.__model.has_object('orbit'):
            bar_pos = self.__model['orbit'].get_barycentre()
            omega = self.__model['orbit']['omega_an']
            incl = self.__model['orbit']['dinc']
        else:
            raise TypeError('The tool is currently working with binaries only '
                            'i.e. for models, where orbit has been defined.')

        # get the physical resolution along the
        # x-axis and size of the grid in pixels
        if self.__model.has_object('grid'):
            npx, npy, npz = self.__model['grid'].get_pixel_size()
        else:
            raise KeyError('Object grid was not found among defined objects.'
                           'It should be alway defined')

        # get the size of the final image - it is padded to get a
        # square format, which is convenient for rotation of the
        # image. Nonetheless padding of the image is not
        # necessary for DFT.
        newsize = self.__image_size
        if newsize is None:
            if npx > npy:
                newsize = npx
            else:
                newsize = npy

        # extract the resolution along x-axis
        phys_res = self.__model['grid']['stepx']

        # compute spatial frequencies
        # get parallax and stepsize
        plx = self.__model['distance']['dd'].to('m')
        step_physical = phys_res.to('m')

        # compute angular step
        step_angular = step_physical.value / plx.value

        # axis of the image
        if newsize % 2 > 1:
            xscale = np.arange(-newsize / 2 + 1, newsize / 2 + 1, 1) * step_angular
        else:
            xscale = np.arange(-newsize / 2, newsize / 2, 1) * step_angular

        # set the effective wavelength to template
        self.set_wavelength(w0=ew.to('nm'), wn=ew.to('nm') + 1.0 * units.nm, step=1.0 * units.nm)

        # runs shellspec
        directory = os.getcwd()
        phase_file = os.path.join(directory, 'phases')
        control_file = os.path.join(directory, 'shellspec.in')
        self.write_phase([phase], filename=phase_file)
        self.write_template(filename=control_file)

        # run shellspec
        self.__run_shellspec()

        # filename with the image
        filename = os.path.join(directory, '2Dimage_001')

        # read the image
        img = self.__read_image(filename)

        # here it only rotates the image and returns it
        img = self.__FT_one_image(img, npx, npy, phys_res, bar_pos, phase, incl, omega, newsize, image_only=True)

        if img_only:
            return img, xscale

        # get image kwargs
        fontsize = img_kwargs.get('fontsize', 16)
        title = img_kwargs.get('title', '')
        vmin = img_kwargs.get('vmin', None)
        vmax = img_kwargs.get('vmax', None)

        # plot the image
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, aspect='equal')
        ax.imshow(img, cmap='gray', extent=[xscale.min(), xscale.max(), xscale.min(), xscale.max()], vmin=vmin,
                  vmax=vmax)
        ax.set_xlim(xscale.min(), xscale.max())
        ax.set_ylim(xscale.min(), xscale.max())
        ax.set_title(title, fontsize=fontsize)
        ax.set_xlabel(r'$\alpha$(rad)', fontsize=fontsize)
        ax.set_ylabel(r'$\delta$(rad)', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        # plt.axes().set_aspect('equal', 'datalim')
        plt.tight_layout()
        plt.savefig(figname)
        plt.close()

    def get_synthetic_light_curve(self, phase, passband=None, eff_wave=None):
        """
        Computes a synthetic light curve.
        :return:
        """
        if eff_wave is not None:
            eff_wave = eff_wave * units.nm

        # set effective wavelength from passband
        if passband is not None:
            idx = filters.index(self.__passband)
            eff_wave = lambda_eff[idx]

        magnitude = self.__pp_get_lc_one_wavelength(directory='.', phase=phase, ew=eff_wave)

        return magnitude

    def plot_comparison(self, filename=None, observable='all', figname=None):
        """
        Plots comparison of data and model.
        :param filename:
        :param observable:
        :param figname:
        :return:
        """
        # make it case insensitive
        observable = observable.lower()

        # we can either plot file by file
        if filename is not None:
            print(filename)
            # first plot squared visibility
            if filename in self.__vis2_comparison['filename'] and observable in ['vis2', 'all']:

                # select corresponding data
                ind = np.where(self.__vis2_comparison['filename'] == filename)[0]

                # extract observables
                ucoord = self.__vis2_comparison['ucoord'][ind]
                vcoord = self.__vis2_comparison['vcoord'][ind]
                eff_wave = self.__vis2_comparison['eff_wave'][ind]
                vis2data = self.__vis2_comparison['vis2data'][ind]
                vis2err = self.__vis2_comparison['vis2err'][ind]
                vis2syn = self.__vis2_comparison['vis2syn'][ind]
                phase = self.__vis2_comparison['phase'][ind]

                # get mean phase and ew
                phmean = phase.mean()
                ewmean = eff_wave.mean()
                # print phase, phmean, ewmean
                img, scale = self.get_image(phase=phmean, eff_wave=ewmean, img_only=True)

                # setup some default figure name
                if figname is None:
                    figname = filename

                outputname = '.'.join([figname, 'vis2.png'])

                # do the plotting
                plot_squared_visibility(outputname, ucoord / eff_wave, vcoord / eff_wave, vis2data, vis2err, vis2syn,
                                        image=img, scale=scale)

            # first plot triple product
            if filename in self.__cp_comparison['filename'] and observable in ['cp', 'all']:

                # select corresponding data
                ind = np.where(self.__cp_comparison['filename'] == filename)[0]

                # extract observables
                u1coord = self.__cp_comparison['u1coord'][ind]
                v1coord = self.__cp_comparison['v1coord'][ind]
                u2coord = self.__cp_comparison['u2coord'][ind]
                v2coord = self.__cp_comparison['v2coord'][ind]
                eff_wave = self.__cp_comparison['eff_wave'][ind]
                t3amp = self.__cp_comparison['t3amp'][ind]
                t3amperr = self.__cp_comparison['t3amperr'][ind]
                t3phi = self.__cp_comparison['t3phi'][ind]
                t3phierr = self.__cp_comparison['t3phierr'][ind]
                t3ampsyn = self.__cp_comparison['t3ampsyn'][ind]
                t3phisyn = self.__cp_comparison['t3phisyn'][ind]

                # setup some default figure name
                if figname is None:
                    figname = filename

                outputname = '.'.join([figname, 'cp.png'])

                # do the plotting
                plot_triple_product(outputname, u1coord / eff_wave, v1coord / eff_wave,
                                    u2coord / eff_wave, v2coord / eff_wave,
                                    t3amp, t3phi, t3amperr, t3phierr, t3ampsyn, t3phisyn)

            # finally plot the light curve
            if filename in self.__lc_comparison['filename'] and observable in ['lc', 'all']:

                # select the corresponding data
                ind = np.where(self.__lc_comparison['filename'] == filename)[0]

                # extract the data
                obsmag = self.__lc_comparison['magnitude'][ind]
                synmag = self.__lc_comparison['magnitudesyn'][ind]
                passband = self.__lc_comparison['passband'][ind]
                obsphase = self.__lc_comparison['phase'][ind]
                error = self.__lc_comparison['error'][ind]
                offset = self.__lc_comparison['offset'][ind]

                # take just name of passband and only one value of
                # the offset
                passband = passband[0]
                offset = offset[0]

                # extract whole lc
                idx = self.__lcsyn['passband'].index(passband)
                whole_lc_phase = self.__lcsyn['phase'][idx]
                whole_lc_mag = self.__lcsyn['magnitude'][idx] + offset

                # setup some default figure name
                if figname is None:
                    figname = filename
                outputname = '.'.join([figname, 'lc.png'])

                plot_light_curve(outputname, obsphase, obsmag, error, synmag,
                                 whole_lc_phase, whole_lc_mag, passband=passband)

    def run_fit(self, fitter='nlopt_nelder_mead', **fit_kwargs):
        """
        Runs the fitting.
        :param fitter:
        :param fit_kwargs:
        :return:
        """
        # make input case-insensitive
        fitter = fitter.lower()

        # get fitted parameters
        fitparams = self.get_fitted_parameters()

        # get the number of fitted parameters
        npar = len(fitparams)
        # print npar

        # check that they do not lie outside
        # the fitted region
        for p in fitparams:
            objname = p['object']
            self.__model[objname].check_boundaries()

        # extract their value and bounds
        values = [x['value'] for x in fitparams]
        vmins = [x['vmin'] for x in fitparams]
        vmaxs = [x['vmax'] for x in fitparams]

        # get the properties of the fitter
        fitter_props = fitter_definitions[fitter]

        if fitter.find('nlopt') > -1:

            # define objective function
            def obj_func(p, grad):
                return self.compute_chi2(p, fitparams)

            # define fitter
            opt = nlopt.opt(fitter_props['environment'], npar)

            # set boundaries
            if fitter_props['uses_bounds']:
                opt.set_lower_bounds(vmins)
                opt.set_upper_bounds(vmaxs)

            # set initial step
            stepsize = (np.array(vmaxs) - np.array(vmins)) / 4.
            opt.set_initial_step(stepsize.tolist())

            # set objective function
            opt.set_min_objective(obj_func)

            # set fitting criteria
            for key in list(fit_kwargs.keys()):
                key = key.lower()
                if key == 'xtol':
                    opt.set_xtol_rel(fit_kwargs[key])
                if key == 'ftol':
                    opt.set_ftol_rel(fit_kwargs[key])
                if key == 'maxfun':
                    opt.set_maxeval(fit_kwargs[key])

            # empty iterations
            self.__iterations = []

            # run the fitting
            result = opt.optimize(values)

        # SciPy fitters
#        elif fitter.find('sp_diff_evol') > -1:
        elif fitter.find('sp_') > -1:

            # extract fitting environment
            opt = fitter_props['environment']

            # define bounds
            if fitter_props['uses_bounds']:
                bounds = [[vmin, vmax] for vmin, vmax in zip(vmins, vmaxs)]
                result = opt(self.compute_chi2, bounds, **fit_kwargs)
                result = result.x
            else:
                raise NotImplementedError('Fitters not using bounds were not implemented.')

        else:
            raise NotImplementedError('The fitter %s is not implemented' % fitter)

        return result

    def set_model_to_shellspec(self):
        """
        Propagates the model values to the template.
        :return:
        """

        # go over all objects
        for objname in list(self.__model.keys()):

            # resolved internal constraints within
            # each object
            # print objname
            self.__model[objname].resolve_constraints()

            # if it is the object - set only the inclination
            if objname in ['orbit']:
                value = self.__model[objname]['dinc']
                self.__set_parameter_to_template('dinc', value)
                continue

            # go over each parameter
            for parname in list(self.__model[objname].keys()):

                # get is value
                value = self.__model[objname][parname]
                # print parname, value, type(value)

                # write it to shellspec control file
                self.__set_parameter_to_template(parname, value)

    def set_parameter(self, name, **kwargs):
        """
        Sets parameter in model to a certain value.
        :param name:
        :param kwargs:
        :return:
        """
        # go over all objects
        for objname in list(self.__model.keys()):
            if name.lower() in list(self.__model[objname].keys()):
                self.__model[objname].set_parameter(name, **kwargs)
                return

        raise KeyError('The parameter %s does not belong to any'
                       ' defined objects. Defined objects are %s.' %
                       (name, str(list(self.__model.keys()))))
    
    def get_parameter(self, name, **kwargs):
        """
        Gets parameter from model.
        :param name:
        :param kwargs:
        :return:
        """
        # go over all objects
        for objname in list(self.__model.keys()):
            if name.lower() in list(self.__model[objname].keys()):
                return self.__model[objname].get_parameter(name, **kwargs)

        raise KeyError('The parameter %s does not belong to any'
                       ' defined objects. Defined objects are %s.' %
                       (name, str(list(self.__model.keys()))))
    
    def set_wavelength(self, w0, wn, step):
        """
        Sets wavelength at which we will compute, 
        wavelengths have to be set in angstrom
        :param w0:
        :param wn:
        :param step:
        """
        # if astropy quatities are pass they 
        # are transformed to angstrom
        if isinstance(w0, units.Quantity):
            w0 = w0.to('Angstrom').value
        if isinstance(wn, units.Quantity):
            wn = wn.to('Angstrom').value
        if isinstance(step, units.Quantity):
            step = step.to('Angstrom').value

        # set the wavelengths to the template
        self.__set_parameter_to_template('alam1', w0)
        self.__set_parameter_to_template('alamn', wn)
        self.__set_parameter_to_template('alams', step)

    def write_iterations(self, filename='fit.log'):
        """
        Writes all iterations.
        :param filename: name of the output file
        :return:
        """
        string = ''
        # write header
        fitparams = self.get_fitted_parameters(attr='name')
        string += '%6s' % 'Niter'
        for name in fitparams:
            string += '%22s' % name
        for name in ['nlc', 'nvis2', 'ncp', 'nt3amp', 'nvisamp', 'nvisphi', 'nsed', 'nspe', 'ntotal']:
            string += '%8s' % name
        for name in ['chi2lc', 'chi2vis2', 'chi2cp', 'chi2t3amp', 'chi2visamp', 'chi2visphi', 'chi2sed', 'chi2spe', 'chi2total', 'r_chi2lc', 'r_chi2vis2', 'r_chi2cp', 'r_chi2t3amp', 'r_chi2visamp', 'r_chi2visphi', 'r_chi2sed', 'r_chi2spe', 'r_chi2total']:
            string += '%16s' % name
        string += '\n'

        # write iterations
        nparams = len(fitparams)
        nterms = 9
        for i in range(len(self.__iterations)):
            string += '%6i' % i
            for j in range(nparams):
                string += '%22.16g' % self.__iterations[i][j]
            for j in range(nterms):
                string += '%8d' % self.__iterations[i][nparams+j]
            for j in range(2*nterms):
                string += '%16.8g' % self.__iterations[i][nparams+nterms+j]
            string += '\n'

        write_file(filename, string)

    def write_model(self):
        """
        Writes the model for each file.
        :return:
        """
        data_structs = [self.__lc_comparison,
                        self.__vis2_comparison,
                        self.__cp_comparison,
                        self.__vis_comparison,
                        self.__sed_comparison,
                        self.__spe_comparison]
        data_types = ['lc', 'vis2', 'cp', 'vis', 'sed', 'spe']

        for ds, dt in zip(data_structs, data_types):

            # extract filename for each file
            filenames = np.unique(ds['filename'])

            for filename in filenames:

                # select correct keys
                if dt == 'lc':
                    keys = ['hjd', 'eff_wave', 'eff_band', 'magnitude', 'error', 'magnitudesyn', 'offset', 'chi2']
                elif dt == 'vis2':
#                    keys = ['ucoord', 'vcoord', 'hjd', 'eff_wave', 'eff_band', 'vis2data', 'vis2err', 'vis2syn', 'chi2']  # we don't have it (yet)
                    keys = ['ucoord', 'vcoord', 'hjd', 'eff_wave', 'vis2data', 'vis2err', 'vis2syn', 'chi2']
                elif dt == 'cp':
#                    keys = ['u1coord', 'v1coord', 'u2coord', 'v2coord', 'hjd', 'eff_wave', 'eff_band', 't3amp', 't3amperr', 't3ampsyn', 't3phi', 't3phierr', 't3phisyn', 'chi2amp', 'chi2phi']
                    keys = ['u1coord', 'v1coord', 'u2coord', 'v2coord', 'hjd', 'eff_wave', 't3amp', 't3amperr', 't3ampsyn', 't3phi', 't3phierr', 't3phisyn', 'chi2amp', 'chi2phi']
                elif dt == 'vis':
                    keys = ['ucoord', 'vcoord', 'hjd', 'eff_wave', 'eff_band', 'visamp', 'visamperr', 'visampsyn', 'visphi', 'visphierr', 'visphisyn', 'mulfac', 'offset', 'slips', 'chi2amp', 'chi2phi']
                elif dt == 'sed' or dt == 'spe':
                    keys = ['hjd', 'eff_wave', 'eff_band', 'flux', 'error', 'fluxsyn', 'dataset', 'phase', 'chi2']

                # select correct data
                ind = np.where(ds['filename'] == filename)[0]

                # group them into block
                block = []
                for key in keys:
#                    print("filename = ", filename)  # dbg
#                    print("key = ", key)  # dbg
#                    print("ind = ", ind)  # dbg
                    block.append(ds[key][ind])

                # substitute all None's by np.nan
                for i in range(len(block)):
                    for j in range(len(block[i])):
                        if block[i][j] == None:
                            block[i][j] = np.nan

                # set filename
                filename = '.'.join([filename, dt, 'syn.dat'])

                # set header
                header = ''.join(['%20s' % key for key in keys])

                # save the file
                np.savetxt(filename, np.column_stack(block), header=header, fmt='%20.12e')

    def write_phase(self, phase, filename='phases'):
        """
        Writes a column of phases.
        :param phase:
        :param filename:
        """
        string = ''
        for ph in phase:
            string += "%s\n" % str(ph)
        write_file(filename, string)

    def write_template(self, filename='shellspec.in'):
        """
        Writes the shellspec template file.
        :param filename:
        """
        write_file(filename, self.__shellspec_template)

    def __compute_visibility(self, fftimg, fu_img, fv_img, fu, fv):
        """
        Computes synthetic visibility.
        :param fftimg: Fourier transform of an image
        :param fu_img: Spatial frequencies along the first axis of the image.
        :param fv_img: Spatial frequencies along the second axis of the image.
        :param fu: Spatial frequencies at which visibility should be computed -- first axis.
        :param fv: Spatial frequencies at which visibility should be computed -- second axis.
        :return: vis: visibility V(fu,fv)
        """

        # group baseline projections
        # obsbase = np.column_stack([fu, fv])

        # split the image into modulus and phase
        imgmod = np.abs(fftimg)
        imgphase = np.angle(fftimg)

        # do the interpolation in modulus and phase angle separately
        # mod = interpn((fu_img, fv_img), imgmod, obsbase, method='splinef2d')
        # phase = interpn((fu_img, fv_img), imgphase, obsbase, method='splinef2d')

        # alternative computation
        f_mod = RectBivariateSpline(fu_img, fv_img, imgmod)
        f_phase = RectBivariateSpline(fu_img, fv_img, imgphase)
        mod = f_mod.ev(fu, fv)
        phase = f_phase.ev(fu, fv)

        # get complex visibility
        vis = mod * np.exp(1j * phase)

        return vis

    def __FT_one_image(self, img, npx, npy, phys_res, bar_pos, phase, incl, omega, newsize, image_only=False, order=1):
        """
        Computes Fast Fourier Transform of one image.
        :param img: Original image
        :param npx:
        :param npy:
        :param phys_res:
        :param bar_pos:
        :param phase:
        :param incl:
        :param omega:
        :param newsize:
        :param image_only: 
        :param order: Interpolation order, 0 .. nearest-neighbor, 1 .. linear, 3 .. cubic spline.
        :return: img: Synthetic image.
        :return: fftimg: Fourier Transform.

        Note: Beware of artefacts if order=3!
        """
        # first reformat the image --- 1 column -> 2d array
        img = img.reshape((npx, npy))

        # get omega in degrees and inclination in radians
        incl = incl.to('rad').value
        omega = omega.to('deg').value

        # shift the centre of the image into barycentre
        # if its position is nonzero
        if bar_pos > 0.:

            # convert the shift into pixel scale
            bar_pos_pix = bar_pos / phys_res

            # convert phase to shellspec angle
            angle = ((phase * 2 - 0.5) * np.pi) % (2 * np.pi)

            # position of the barycentre in the rotated frame
            xc, yc, zc = rotateZ(bar_pos_pix, 0., 0., -angle)
            xcc, ycc, zcc = rotateX(xc, yc, zc, -incl)
            print("__FT_one_image: xcc = ", xcc, " ycc = ", ycc, " bar_pos_pix = ", bar_pos_pix)  # dbg

            # move the centre of the image into barycentre
            img = ndimage.shift(img, (-xcc, -ycc), order=order)

        # pad the image with zeros to transform it into square
        # odd number of pixels along x-axis
        if npx < newsize:
            if npx % 2 > 0:
                padx = ((newsize - npx) // 2 + 1, (newsize - npx) // 2)
            # even number
            else:
                padx = ((newsize - npx) // 2, (newsize - npx) // 2)
        else:
            padx = (0, 0)
        # the same y-axis for odd
        if npy < newsize:
            if npy % 2 > 0:
                pady = ((newsize - npy) // 2 + 1, (newsize - npy) // 2)
            # even
            else:
                pady = ((newsize - npy) // 2, (newsize - npy) // 2)
        else:
            pady = (0, 0)

        # pad the image with zeros
        img = np.lib.pad(img, (padx, pady), 'constant', constant_values=(0.0,0.0))

        # shift the image so the barycentre is in the exact centre, i.e. between zero pixels
        # use up-scaling (by 2) and nearest-neighbor (0) to prevent blurring, i.e. fake limb darkening!
        # img = ndimage.shift(img, (-0.5, -0.5), order=order)
        img = ndimage.zoom(img, 2, order=0)
        img = ndimage.shift(img, (-1, -1), order=order)

        # rotate the image according to the ascending node
        img = ndimage.rotate(img, omega, order=order, reshape=False)

        # move the barycentre back to the centre of the pixel newsize // 2 in both axes
        # use down-scaling; we prefer linear (1) here due to ~45 deg angles 
        # img = ndimage.shift(img, (0.5, 0.5), order=order)
        img = ndimage.shift(img, (1, 1), order=order)
        img = ndimage.zoom(img, 0.5, order=order)

        # if we want only prepared image --- for DFT for example
        if image_only:
            return img / img.sum()

        # do NOT transpose the image here; otherwise it's misleading...
        #img = img.T

        # store the image for consequent plotting/debugging
        if self.debug:
            self.debug_image = img

        # compute the normalized Fourier transform
        fftimg = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(img)))

        return fftimg / np.max(np.abs(fftimg))

    def __get_if_observables_DFT(self, img, ew, phase, xscale, yscale, if_ew_precision, if_phase_precision):
        """
        Computes interferometric observables.
        :param img: The object image at certain phase and wavelength.
        :param ew: The effective wavelength
        :param phase: The orbital phase
        :param xscale: The angular scale of the image along East-West axis.
        :param yscale: The angular scale of the image along North-South axis.
        :param if_ew_precision: Resolution in effective wavelength.
        :param if_phase_precision: Resolution in phase.
        :return:
        """

        # select observed deps and indeps corresponding to ew and phase
        # for squared visibility
        if len(self.__vis2_comparison['eff_wave']) == 0:
            visind = []
        else:
            diff_ph = np.absolute(self.__vis2_comparison['phase'] - phase)
            diff_ew = np.absolute(self.__vis2_comparison['eff_wave'] - ew)
            visind = np.where((diff_ph < 10 ** -if_phase_precision) & (diff_ew < 10 ** -if_ew_precision))[0]

        # compute squared visibility
        if len(visind) > 0:
            # extract independent variable baseline
            ucoord = self.__vis2_comparison['ucoord'][visind]
            vcoord = self.__vis2_comparison['vcoord'][visind]

            # effective wavelength
            eff_wave = self.__vis2_comparison['eff_wave'][visind]

            # compute visibility
            vis = np.zeros(ucoord.size, dtype='complex64')
            for i in range(0, ucoord.size):
                vis[i] = DFT2d(img, xscale, yscale, ucoord[i] / eff_wave[i], vcoord[i] / eff_wave[i])

            # compute squared visibility
            vis2syn = np.abs(vis) ** 2

            # assign it to list of comparisons
            self.__vis2_comparison['vis2syn'][visind] = vis2syn

        # select observed deps and indeps corresponding to ew and phase
        # for closure phase
        if len(self.__cp_comparison['eff_wave']) == 0:
            cpind = []
        else:
            diff_ph = np.absolute(self.__cp_comparison['phase'] - phase)
            diff_ew = np.absolute(self.__cp_comparison['eff_wave'] - ew)
            cpind = np.where((diff_ph < 10 ** -if_phase_precision) & (diff_ew < 10 ** -if_ew_precision))[0]

        # compute closure phase
        # extract independent variables
        if len(cpind) > 0:

            # first baseline
            u1coord = self.__cp_comparison['u1coord'][cpind]
            v1coord = self.__cp_comparison['v1coord'][cpind]

            # second baseline
            u2coord = self.__cp_comparison['u2coord'][cpind]
            v2coord = self.__cp_comparison['v2coord'][cpind]

            # sthird baseline
            u3coord = -(u2coord + u1coord)
            v3coord = -(v2coord + v1coord)

            # effective wavelength
            eff_wave = self.__cp_comparison['eff_wave'][cpind]

            # empty arrays for the output
            visb1 = np.zeros(u1coord.size, dtype='complex64')
            visb2 = np.zeros(u1coord.size, dtype='complex64')
            visb3 = np.zeros(u1coord.size, dtype='complex64')

            # carry out the computation for each baseline
            for i in range(0, u1coord.size):
                visb1[i] = DFT2d(img, xscale, yscale, u1coord[i] / eff_wave[i], v1coord[i] / eff_wave[i])
                visb2[i] = DFT2d(img, xscale, yscale, u2coord[i] / eff_wave[i], v2coord[i] / eff_wave[i])
                visb3[i] = DFT2d(img, xscale, yscale, u3coord[i] / eff_wave[i], v3coord[i] / eff_wave[i])

            # compute triple product -0.5286E+13  0.6955E+11  0.00000E+000
            t3syn = visb1 * visb2 * visb3

            # get its amplitude and phase
            t3ampsyn = np.abs(t3syn)
            t3phisyn = np.angle(t3syn)

            # unwrap phase
            mean_phase = np.mean(t3phisyn)
            t3phisyn = np.angle(np.exp(1j * t3phisyn) * np.exp(-1j * mean_phase)) + mean_phase
            t3phisyn = np.degrees(t3phisyn)

            # assign it into comparison array
            self.__cp_comparison['t3ampsyn'][cpind] = t3ampsyn
            self.__cp_comparison['t3phisyn'][cpind] = t3phisyn


        # visibility amplitude and phase
        if len(self.__vis_comparison['eff_wave']) == 0:
            visind = []
        else:
            diff_ph = np.absolute(self.__vis_comparison['phase'] - phase)
            diff_ew = np.absolute(self.__vis_comparison['eff_wave'] - ew)
            visind = np.where((diff_ph < 10 ** -if_phase_precision) & (diff_ew < 10 ** -if_ew_precision))[0]

        if len(visind) > 0:
            ucoord = self.__vis_comparison['ucoord'][visind]
            vcoord = self.__vis_comparison['vcoord'][visind]

            eff_wave = self.__vis_comparison['eff_wave'][visind]

            vis = np.zeros(ucoord.size, dtype='complex64')
            for i in range(0, ucoord.size):
                vis[i] = DFT2d(img, xscale, yscale, ucoord[i] / eff_wave[i], vcoord[i] / eff_wave[i])

            visampsyn = np.abs(vis)
            visphisyn = np.angle(vis)

            # 2DO: unwrap phase <- shall we use it or NOT? <- NOT (already unwrapped in the pipeline)
#            mean_phase = np.mean(visphisyn)
#            visphisyn = np.angle(np.exp(1j * visphisyn) * np.exp(-1j * mean_phase)) + mean_phase
            visphisyn = np.degrees(visphisyn)

            # 2DO: compute differential visibility?! <- if YES, we need offsets of individual files (see LC below)

            self.__vis_comparison['visampsyn'][visind] = visampsyn
            self.__vis_comparison['visphisyn'][visind] = visphisyn

    def __get_if_observables_FFT(self, fftimg, ew, phase, fu_img, fv_img, if_ew_precision, if_phase_precision):
        """
        ****************
        NOT WORKING WELL
        ****************
        Computes interferometric observables.
        :param fftimg: Fourier transform of the object image at certain phase and wavelength.
        :param ew: The effective wavelength
        :param phase: The orbital phase
        :param fu_img: The spatial frequencies corresponding to each point of the image along East-West axis.
        :param fv_img: The spatial frequencies corresponding to each point of the image along North-South axis.
        :param if_ew_precision: Resolution in effective wavelength.
        :param if_phase_precision: Resolution in phase.
        :return:
        """

        # select observed deps and indeps corresponding to ew and phase
        # for squared visibility
        if len(self.__vis2_comparison['eff_wave']) == 0:
            visind = []
        else:

            diff_ph = np.absolute(self.__vis2_comparison['phase'] - phase)
            diff_ew = np.absolute(self.__vis2_comparison['eff_wave'] - ew)
            visind = np.where((diff_ph < 10 ** -if_phase_precision) & (diff_ew < 10 ** -if_ew_precision))[0]

        # compute squared visibility
        if len(visind) > 0:
            # extract independent variables
            # baseline
            ucoord = self.__vis2_comparison['ucoord'][visind]
            vcoord = self.__vis2_comparison['vcoord'][visind]

            # effective wavelength
            eff_wave = self.__vis2_comparison['eff_wave'][visind]

            # compute visibility
            vis = self.__compute_visibility(fftimg, fu_img, fv_img, ucoord / eff_wave, vcoord / eff_wave)

            # compute squared visibility
            vis2syn = np.abs(vis) ** 2

            # assign it to list of comparisons
            self.__vis2_comparison['vis2syn'][visind] = vis2syn

            # select observed deps and indeps corresponding to ew and phase
            # for closure phase
        if len(self.__cp_comparison['eff_wave']) == 0:
            cpind = []
        else:
            diff_ph = np.absolute(self.__cp_comparison['phase'] - phase)
            diff_ew = np.absolute(self.__cp_comparison['eff_wave'] - ew)
            cpind = np.where((diff_ph < 10 ** -if_phase_precision) & (diff_ew < 10 ** -if_ew_precision))[0]

        # compute closure phase
        # extract independent variables
        if len(cpind) > 0:

            # first baseline
            u1coord = self.__cp_comparison['u1coord'][cpind]
            v1coord = self.__cp_comparison['v1coord'][cpind]

            # second baseline
            u2coord = self.__cp_comparison['u2coord'][cpind]
            v2coord = self.__cp_comparison['v2coord'][cpind]

            # sthird baseline
            u3coord = -(u2coord + u1coord)
            v3coord = -(v2coord + v1coord)

            # effective wavelength
            eff_wave = self.__cp_comparison['eff_wave'][cpind]

            # carry out the computation for each baseline
            visb1 = self.__compute_visibility(fftimg, fu_img, fv_img, u1coord / eff_wave, v1coord / eff_wave)
            visb2 = self.__compute_visibility(fftimg, fu_img, fv_img, u2coord / eff_wave, v2coord / eff_wave)
            visb3 = self.__compute_visibility(fftimg, fu_img, fv_img, u3coord / eff_wave, v3coord / eff_wave)

            # compute triple product
            t3syn = visb1 * visb2 * visb3

            # get its amplitude and phase
            t3ampsyn = np.abs(t3syn)
            t3phisyn = np.angle(t3syn)

            # unwrap phase
            mean_phase = np.mean(t3phisyn)
            t3phisyn = np.angle(np.exp(1j * t3phisyn) * np.exp(-1j * mean_phase)) + mean_phase
            t3phisyn = np.degrees(t3phisyn)

            # assign it into comparison array
            self.__cp_comparison['t3ampsyn'][cpind] = t3ampsyn
            self.__cp_comparison['t3phisyn'][cpind] = t3phisyn

            # 2DO: we do NOT have visibility amplitude and phase here (see __get_if_observables_DFT)...

    def __get_lc_observables(self, lcsyn):
        """
        Evaluates synthetic light curve at desired  phase points.
        :return:
        """

        # go over each passband
        for i, band in enumerate(lcsyn['passband']):

            # synthetic phase and magnitude
            synmag = lcsyn['magnitude'][i]
            synphase = lcsyn['phase'][i]

            # mirror the flux for phase = 0.0
            synmag = np.append(synmag, synmag[0])
            synphase = np.append(synphase, 1.0)

            # extract all data for a given passband
            ind = np.where(self.__lc_comparison['passband'] == band)[0]

            # get filenames and observed phases
            filenames = self.__lc_comparison['filename'][ind]
            obsphase = self.__lc_comparison['phase'][ind]
            obsmag = self.__lc_comparison['magnitude'][ind]

            # interpolate to observed phases
            synmag_itp = interpolate1d(obsphase, synphase, synmag)

            # now go over each file and determine the offset
            for filename in lcsyn['filename'][i]:

                # select data from one file
                find = np.where(filenames == filename)[0]

                # get synthetic and observed magnitudes
                obsmag_one_file = obsmag[find]
                synmag_one_file = synmag_itp[find]

                # convert magnitude monochromatic to passband (m_lambda to m_band)
                # in lightcurve, unit of F_lambda is erg s^-1 cm^-2 cm^-1
                # in definitions.py, unit of F_calib is J s^-1 m^-2 m^-1
                m_lambda = synmag_one_file  
                idx = filters.index(band)

                erg = 1.e-7
                cm = 1.e-2
                F_lambda = 10.**(-0.4*m_lambda) * erg*cm**(-2)*cm**(-1)
                F_band = F_lambda * eff_band[idx]
                m_band = -2.5 * np.log10(F_band/(calibration_flux[idx] * eff_band[idx]))

#                if band == 'johnson_v':
#                    print("band = " + str(band))
#                    print("m_lambda = " + str(m_lambda[0]))
#                    print("F_lambda = " + str(F_lambda[0]))
#                    print("idx = " + str(idx))
#                    print("eff_band = " + str(eff_band[idx]))
#                    print("calibration_flux = " + str(calibration_flux[idx]))
#                    print("F_band = " + str(F_band[0]))
#                    print("m_band = " + str(m_band[0]))
#                    sys.exit(1)

                # compute additonal offset
                if self.__use_offset:
                    offset = get_offset(obsmag_one_file, m_band)
                    offset = offset.x
                else:
                    offset = 0.0

                # assign synthetic data to comparison list
                self.__lc_comparison['magnitudesyn'][ind[find]] = m_band + offset
                self.__lc_comparison['offset'][ind[find]] = np.ones(len(m_band)) * offset

    def __pp_get_images_one_wavelength(self, directory, ew, phase, npx, npy, phys_res, bar_pos,
                                       incl, omega, newsize, queue, image_only=True):
        """
        Runs shellspec for one wavelength
        :param ew:
        :param phase:
        :param npx:
        :param npy:
        :param phys_res:
        :param bar_pos:
        :param incl:
        :param omega:
        :param newsize:
        :param queue:
        :param image_only:
        :return:
        """

        # set the effective wavelength to template
        self.set_wavelength(w0=ew, wn=ew + 1.0 * units.nm, step=1.0 * units.nm)

        # set limb-darkening coefficients
        dlst = self.limcofst.interp(ew.to('m').value)
        dlcp = self.limcofcp.interp(ew.to('m').value)
        self.__set_parameter_to_template('dlst', dlst)
        self.__set_parameter_to_template('dlcp', dlcp)

        # runs shellspec
        phase_file = os.path.join(directory, 'phases')
        control_file = os.path.join(directory, 'shellspec.in')
        self.write_phase(phase, filename=phase_file)
        self.write_template(filename=control_file)

        # run shellspec
        cwd = os.getcwd()
        os.chdir(directory)
        exitcode = self.__run_shellspec()
        os.chdir(cwd)

        if exitcode != 0:
            if queue is not None:
                queue.put(None)
            sys.exit(exitcode)

        # image length
        # imglen = npx * npy

        images = []
        # read the images - one for each phase
        for j in range(0, len(phase)):
            # set the filename
            filename = os.path.join(directory, '2Dimage_%03d' % (j + 1))
            # read the image
            img = self.__read_image(filename)
            # img = self.__read_binary_image(filename, imglen)

            # here it only rotates the image and returns it
            img = self.__FT_one_image(img, npx, npy, phys_res, bar_pos,
                                      phase[j], incl, omega, newsize, image_only=image_only)
            # in case of debug mode, save both - the image and its Fourier transform
            if self.debug:
                # the image
                figname = "img_%.2f_%.4f.png" % (ew.to('Angstrom').value, phase[j])
                plt.imshow(img, cmap='gray')
                plt.xlabel('x [pxl]')
                plt.ylabel('y [pxl]')
                plt.colorbar()
                plt.tight_layout()
                plt.savefig(figname)
                plt.close()
                figname = "img_%.2f_%.4f.dat" % (ew.to('Angstrom').value, phase[j])
                np.savetxt(figname, img)

            images.append(img)

        # save the result into queue
        queue.put([ew, phase, images])

    def __pp_get_lc_one_wavelength(self, directory, ew, phase, queue=None):
        """
        Computes one light curve.
        :param directory:
        :param ew:
        :param phase:
        :param queue:
        :return:
        """

        # set the effective wavelength to template
        self.set_wavelength(w0=ew, wn=ew + 1.0 * units.nm, step=1.0 * units.nm)

        # set limb-darkening coefficients
        dlst = self.limcofst.interp(ew.to('m').value)
        dlcp = self.limcofcp.interp(ew.to('m').value)
        self.__set_parameter_to_template('dlst', dlst)
        self.__set_parameter_to_template('dlcp', dlcp)

        # write template and phases
        phase_file = os.path.join(directory, 'phases')
        control_file = os.path.join(directory, 'shellspec.in')
        self.write_phase(phase, filename=phase_file)
        self.write_template(filename=control_file)

        # run shellspec
        cwd = os.getcwd()
        os.chdir(directory)
        exitcode = self.__run_shellspec()
        os.chdir(cwd)

        if exitcode != 0:
            if queue is not None:
                queue.put(None)
            sys.exit(exitcode)

        # read the light curve
        mags = self.__read_light_curve(filename=os.path.join(directory, 'lightcurve'))

        # store the result of parallel computing
        if queue is not None:
            queue.put(mags)
        else:
            return mags

    def __process_wavelength_DFT(self, directory, ew, phase, npx, npy, phys_res, bar_pos,
                                 incl, omega, newsize, xscale, yscale, image_only=True):
        """
        Runs shellspec for one wavelength
        :param ew:
        :param phase:
        :param npx:
        :param npy:
        :param phys_res:
        :param bar_pos:
        :param incl:
        :param omega:
        :param newsize:
        :param xscale:
        :param yscale:
        :param image_only:
        :return:
        """

        # set the effective wavelength to template
        self.set_wavelength(w0=ew, wn=ew + 1.0 * units.nm, step=1.0 * units.nm)

        # set limb-darkening coefficients
        dlst = self.limcofst.interp(ew.to('m').value)
        dlcp = self.limcofcp.interp(ew.to('m').value)
        self.__set_parameter_to_template('dlst', dlst)
        self.__set_parameter_to_template('dlcp', dlcp)

        # runs shellspec
        phase_file = os.path.join(directory, 'phases')
        control_file = os.path.join(directory, 'shellspec.in')
        self.write_phase(phase, filename=phase_file)
        self.write_template(filename=control_file)

        # run shellspec
        cwd = os.getcwd()
        os.chdir(directory)
        exitcode = self.__run_shellspec()
        os.chdir(cwd)

        # image size
        # imglen = npx*npy

        # read the images - one for each phase
        for j in range(0, len(phase)):
            # set the filename
            filename = os.path.join(directory, '2Dimage_%03d' % (j + 1))
            # read the image
            img = self.__read_image(filename)
            # img = self.__read_binary_image(filename, imglen)
            # here it only rotates the image and returns it
            img = self.__FT_one_image(img, npx, npy, phys_res, bar_pos,
                                      phase[j], incl, omega, newsize, image_only=image_only)
            # in case of debug mode, save both - the image and its Fourier transform
            if self.debug:
                # the image
                figname = '.'.join(['img', '%.2f' % (ew.to('Angstrom').value), str(phase[j]), 'png'])
                plt.imshow(img, cmap='gray')
                plt.savefig(figname)
                plt.close()

            # compute observables
            self.__get_if_observables_DFT(img, ew.value, phase[j], xscale, yscale, self.__if_ew_precision,
                                          self.__if_phase_precision)

    def __read_binary_image(self, f, imglen):
        """
        Reads shellspec image in binary format
        :param f:
        :param imglen:
        :return:
        """
        # open the image file
        ifile = FortranFile(f, 'r')

        # empty image list
        img = []

        # get the image column
        for i in range(imglen):
            rec = ifile.read_record('f8')
            img.append(rec[-1])

        return np.array(img)

    def __read_image(self, f):
        """
        Reads an image produced by shellspec.
        :param f: file produced by shellspec, where the image is stored
        :return: img -- the image
        """
        # read the file
        ifile = open(f, 'r')
        lines = ifile.readlines()
        ifile.close()

        # translate into image
        img = []
        for i in range(len(lines)):
            data = lines[i].split()

            # throw away empty lines
            if len(data) == 3:
                img.append(float(data[-1]))

        return np.array(img)

    def __read_light_curve(self, filename='lightcurve'):
        """
        Reads light curve produced by shellspec, but 
        only records for the first wavelength.

        lightcurve format:
        0 ... phase
        1 ... radial velocity [km s^-1]
        2 ... magnitude [mag]; m = -2.5 log_10 F_lambda, F_lambda = 10^{-0.4 m}
        3 ... lambda [Ang]
        4 ... monochromatic flux F_nu [erg s^-1 cm^-2 Hz^-1]; F_lambda = c/lambda^2 F_nu

        """

        # open the light curve
        ifile = open(filename, 'r')
        lines = ifile.readlines()
        ifile.close()

        lc = []
        for l in lines:
            d = l.split()
            # discard empty lines
            if len(d) < 5:
                continue
            # read the full ones
            else:
                lc.append(float(d[2]))

        lc = np.array(lc)
        # print lc
        
        return lc[0:lc.size:2]

    def __read_parameter_from_template(self, name):
        """
        Updates the model based on the value in
        the template.
        :return:
        """
        for i, l in enumerate(self.__shellspec_template):
            # first find row with the record
            if l.lower().find(name) > - 1:

                # second find the record within the row
                l = l.replace('#', '').lower().split()
                for j in range(len(l)):

                    # go over each record
                    if l[j].find(name) > -1:

                        # adapt for cases where unit in parentheses
                        # is given
                        if l[j].find('[') > -1:
                            strmax = l[j].find('[')
                        else:
                            strmax = len(l[j])

                        # try to adapt for the cases where one
                        # parameter name is a subset of another one
                        if name != l[j][:strmax].lower():
                            continue

                        # split the line below
                        line = self.__shellspec_template[i+1].split()

                        return line[j]

    def __read_template(self):
        """
        Reads the template and saves two copies.
        One will be kept in its original state.
        """

        # read the file
        ifile = open(self.__shellspec_template_file, 'r')
        self.__shellspec_template = ifile.readlines()
        ifile.close()

        # make a backup copy
        self.__shellspec_template_default = copy.deepcopy(self.__shellspec_template)

    def __ready_comparison(self, observable):
        """
        Prepares a dictionary where the observations and the model
        will be stored.
        :return:
        """
        # closure phase
        if observable == 'cp':
            # create empty dictionary
            compdict = dict(filename=[],
                            u1coord=np.array([]),
                            v1coord=np.array([]),
                            u2coord=np.array([]),
                            v2coord=np.array([]),
                            hjd=np.array([]),
                            phase=np.array([]),
                            eff_wave=np.array([]),
                            eff_band=np.array([]),
                            t3amp=np.array([]),
                            t3amperr=np.array([]),
                            t3phi=np.array([]),
                            t3phierr=np.array([]),
                            t3ampsyn=np.array([]),
                            t3phisyn=np.array([]),
                            chi2amp=np.array([]),
                            chi2phi=np.array([]),
                            weight=np.array([])
                            )

            # readout all data files
            for i in range(self.__data.get_observation_number('if')):
                # stop if closure phase is unwanted
                if self.__exclude_cp and self.__exclude_t3amp:
                    break

                # extract data from one file
                one_file = self.__data.get_observation('if', i).get_data(dtype='cp')
                filename = self.__data.get_observation('if', i).get_filename().split('/')[-1]

                # just in case there are no data
                if one_file is None:
                    continue

                # append all data
                for key in list(one_file.keys()):
                    compdict[key] = np.append(compdict[key], one_file[key])

                # data length
                one_file_lenght = len(one_file['hjd'])

                # append filenames
                compdict['filename'].extend([filename] * one_file_lenght)

                # append empty array for synthetic triple product
                compdict['t3ampsyn'] = np.append(compdict['t3ampsyn'], np.zeros(one_file_lenght))
                compdict['t3phisyn'] = np.append(compdict['t3phisyn'], np.zeros(one_file_lenght))
                compdict['phase'] = np.append(compdict['phase'], np.zeros(one_file_lenght))
                compdict['chi2amp'] = np.append(compdict['chi2amp'], np.zeros(one_file_lenght))
                compdict['chi2phi'] = np.append(compdict['chi2phi'], np.zeros(one_file_lenght))
                compdict['weight'] = np.append(compdict['weight'], np.ones(one_file_lenght) * self.__data.get_observation('if', i).weight_t3amp)

            # make filenames array
            compdict['filename'] = np.array(compdict['filename'])
            return compdict

        # lightcurves
        elif observable == 'lc':

            # create empty comparison dictionary
            compdict = dict(filename=[],
                            passband=[],
                            hjd=np.array([]),
                            magnitude=np.array([]),
                            error=np.array([]),
                            eff_wave=np.array([]),
                            eff_band=np.array([]),
                            offset=np.array([]),
                            magnitudesyn=np.array([]),
                            phase=np.array([]),
                            chi2=np.array([])
                            )

            # empty dictionary for synthetic light curves
            lcsyn = dict(eff_wave=[],
                         magnitude=[],
                         passband=[],
                         filename=[]
                         )

            # readout all data files
            for i in range(self.__data.get_observation_number('lc')):

                # extract observational data from one file
                # (including eff_wave, eff_band, error, hjd, magnitude - see observations.py)
                one_file = self.__data.get_observation('lc', i).get_data()
                filename = self.__data.get_observation('lc', i).get_filename().split('/')[-1]
                passband = self.__data.get_observation('lc', i).get_passband()

                # just in case there are no data
                if one_file is None:
                    continue

                # append the passband
                if passband not in lcsyn['passband']:
                    lcsyn['passband'].append(passband.lower())
                    lcsyn['eff_wave'].append(one_file['eff_wave'][0])
                    lcsyn['filename'].append([filename])

                # if it exists just append filenames - it is necessary
                # for computation of offsets
                else:
                    idx = lcsyn['passband'].index(passband)
                    lcsyn['filename'][idx].append(filename)

                # append all observational data
                for key in list(one_file.keys()):
                    value = one_file[key][0]
                    if isinstance(value, units.Quantity):
                        values = list(map(lambda x: x.value, one_file[key]))
                    else:
                        values = one_file[key]
                    compdict[key] = np.append(compdict[key], values)
#                    print("key = ", key)  # dbg
#                    print("one_file = ", one_file[key])  # dbg
#                    print("compdict = ", compdict[key])  # dbg

                # data length
                one_file_lenght = len(one_file['hjd'])

                # append filenames and passbands
                compdict['filename'].extend([filename] * one_file_lenght)
                compdict['passband'].extend([passband.lower()] * one_file_lenght)

                # append empty array for synthetic magnitudes
                compdict['magnitudesyn'] = np.append(compdict['magnitudesyn'], np.zeros(one_file_lenght))

                # append empty array for offsets between synthetic and observed curve
                compdict['offset'] = np.append(compdict['offset'], np.zeros(one_file_lenght))

                # append empty array for phase
                compdict['phase'] = np.append(compdict['phase'], np.zeros(one_file_lenght))

                compdict['chi2'] = np.append(compdict['chi2'], np.zeros(one_file_lenght))

            # make passbands array
            compdict['passband'] = np.array(compdict['passband'])

            # make filenames array
            compdict['filename'] = np.array(compdict['filename'])
            return compdict, lcsyn

        # squared visibility
        elif observable == 'vis2':

            # create empty dictionary
            compdict = dict(filename=[],
                            ucoord=np.array([]),
                            vcoord=np.array([]),
                            hjd=np.array([]),
                            phase=np.array([]),
                            eff_wave=np.array([]),
                            eff_band=np.array([]),
                            vis2data=np.array([]),
                            vis2err=np.array([]),
                            vis2syn=np.array([]),
                            chi2=np.array([]),
                            weight=np.array([])
                            )

            # readout all data files
            for i in range(self.__data.get_observation_number('if')):

                # break if square visibility is unwanted
                if self.__exclude_vis2:
                    break

                # extract data from one file
                one_file = self.__data.get_observation('if', i).get_data(dtype='vis2')
                filename = self.__data.get_observation('if', i).get_filename().split('/')[-1]

                # just in case there are no data
                if one_file is None:
                    continue

                # append all data
                for key in list(one_file.keys()):
                    compdict[key] = np.append(compdict[key], one_file[key])

                # data length
                one_file_lenght = len(one_file['hjd'])

                # append filenames
                compdict['filename'].extend([filename] * one_file_lenght)

                # append empty array for squared visibility
                compdict['vis2syn'] = np.append(compdict['vis2syn'], np.zeros(one_file_lenght))
                compdict['phase'] = np.append(compdict['phase'], np.zeros(one_file_lenght))
                compdict['chi2'] = np.append(compdict['chi2'], np.zeros(one_file_lenght))
                compdict['weight'] = np.append(compdict['weight'], np.ones(one_file_lenght) * self.__data.get_observation('if', i).weight_vis2)

            # make filenames array
            compdict['filename'] = np.array(compdict['filename'])
            return compdict

        # visibility amplitude and phase
        elif observable == 'vis':

            # create empty dictionary
            compdict = dict(filename=[],
                            ucoord=np.array([]),
                            vcoord=np.array([]),
                            hjd=np.array([]),
                            phase=np.array([]),
                            eff_wave=np.array([]),
                            eff_band=np.array([]),
                            visamp=np.array([]),
                            visamperr=np.array([]),
                            visphi=np.array([]),
                            visphierr=np.array([]),
                            visampsyn=np.array([]),
                            visphisyn=np.array([]),
                            mulfac=np.array([]),
                            offset=np.array([]),
                            slips=np.array([]),
                            chi2amp=np.array([]),
                            chi2phi=np.array([]),
                            weight=np.array([]),
                            visamp_=np.array([]),
                            visphi_=np.array([])
                            )

            # readout all data files
            for i in range(self.__data.get_observation_number('if')):

                # break if square visibility is unwanted
                if self.__exclude_visamp and self.__exclude_visphi:
                    break

                # extract data from one file
                one_file = self.__data.get_observation('if', i).get_data(dtype='vis')
                filename = self.__data.get_observation('if', i).get_filename().split('/')[-1]

                # just in case there are no data
                if one_file is None:
                    continue

                # append all data
                for key in list(one_file.keys()):
                    compdict[key] = np.append(compdict[key], one_file[key])

                # data length
                one_file_lenght = len(one_file['hjd'])

                # append filenames
                compdict['filename'].extend([filename] * one_file_lenght)

                # append empty array for visibility amplitude and phase
                compdict['visampsyn'] = np.append(compdict['visampsyn'], np.zeros(one_file_lenght))
                compdict['visphisyn'] = np.append(compdict['visphisyn'], np.zeros(one_file_lenght))
                compdict['phase'] = np.append(compdict['phase'], np.zeros(one_file_lenght))
                compdict['mulfac'] = np.append(compdict['mulfac'], np.ones(one_file_lenght))
                compdict['offset'] = np.append(compdict['offset'], np.zeros(one_file_lenght))
                compdict['slips'] = np.append(compdict['slips'], np.zeros(one_file_lenght))
                compdict['chi2amp'] = np.append(compdict['chi2amp'], np.zeros(one_file_lenght))
                compdict['chi2phi'] = np.append(compdict['chi2phi'], np.zeros(one_file_lenght))
                compdict['weight'] = np.append(compdict['weight'], np.ones(one_file_lenght) * self.__data.get_observation('if', i).weight_vis2)

            # make filenames array
            compdict['filename'] = np.array(compdict['filename'])

            # backup of visamp and visphi (due to mulfac, offset and slips)
            compdict['visamp_'] = np.append(compdict['visamp_'], compdict['visamp'])
            compdict['visphi_'] = np.append(compdict['visphi_'], compdict['visphi'])
            return compdict

        # SED and spectra
        elif observable == 'sed' or observable == 'spe':

            compdict = dict(filename=[],
                            hjd=np.array([]),
                            eff_wave=np.array([]),
                            eff_band=np.array([]),
                            flux=np.array([]),
                            error=np.array([]),
                            dataset=np.array([]),
                            fluxsyn=np.array([]),
                            phase=np.array([]),
                            chi2=np.array([])
                            )

            dtypesyn = dict(eff_wave=[],
                         eff_wavesyn=[],
                         fluxsyn=[],
                         filename=[]
                         )

            dtype = observable
            for i in range(self.__data.get_observation_number(dtype)):

                one_file = self.__data.get_observation(dtype, i).get_data()
                filename = self.__data.get_observation(dtype, i).get_filename().split('/')[-1]

                if one_file is None:
                    continue
                for key in list(one_file.keys()):
                    compdict[key] = np.append(compdict[key], one_file[key])

                one_file_lenght = len(one_file['hjd'])
                compdict['filename'].extend([filename] * one_file_lenght)
                compdict['fluxsyn'] = np.append(compdict['fluxsyn'], np.zeros(one_file_lenght))
                compdict['phase'] = np.append(compdict['phase'], np.zeros(one_file_lenght))
                compdict['chi2'] = np.append(compdict['chi2'], np.zeros(one_file_lenght))

            compdict['filename'] = np.array(compdict['filename'])
            return compdict, dtypesyn
        

    def __symlink(self, src, dst):
        if not os.path.exists(dst):
            os.symlink(src, dst)

    def __remove(self, dst):
        if os.path.exists(dst):
            os.remove(dst)

    def __run_shellspec(self):
        """
        Runs the shellspec.
        """

        if self.dry_run == True:
            print("Warning: This is a dry-run, no shellspec computation is performed!")
            return 0

        if self.overwrite == False and os.path.exists('shellspec.out'):
            print("Warning: shellspec.out already exists and overwrite is set to False!")
            return 0

        # get the current working directory
        cwd = os.getcwd()
       
        # remove old files (some of)
        self.__remove('lightcurve')
        self.__remove('shellspec.out')
        self.__remove('shellspectrum')

        # copy executable and abundances to the current directory
        shutil.copy2(self.__shellspec_executable_file, cwd)
        shutil.copy2(self.__shellspec_abundance_file, os.path.join(cwd, 'abundances'))

        # symlink files for synthetic spectra and line transfer
        self.__symlink('../pyterpoldata', 'pyterpoldata')
        self.__symlink('../starspec1', 'starspec1')
        self.__symlink('../starspec2', 'starspec2')
        self.__symlink('../line.dat', 'line.dat')
        
        # run shellspec
        cmd = [os.path.join(cwd, 'shellspec')]
        exitcode = subprocess.call(cmd)

        return exitcode

    def __set_parameter_to_template(self, name, value):
        """
        Adjusts a parameter within the template.
        """
        for i, l in enumerate(self.__shellspec_template):
            # first find row with the record
            if l.lower().find(name) > - 1:

                # second find the record within the row
                l = l.replace('#', '').lower().split()
                for j in range(len(l)):

                    # go over each record
                    if l[j].find(name) > -1:

                        # adapt for cases where unit in parentheses
                        # is given
                        if l[j].find('[') > -1:
                            strmax = l[j].find('[')
                        else:
                            strmax = len(l[j])

                        # try to adapt for the cases where one
                        # parameter name is a subset of another one
                        if name != l[j][:strmax].lower():
                            continue
                            
                        # split the line below    
                        line = self.__shellspec_template[i+1].split()
                        
                        # insert the value
                        if isinstance(value, units.Quantity):
                            line[j] = str(value.value)
                        else:
                            line[j] = str(value)
                        
                        # put it back into the template 
                        self.__shellspec_template[i+1] = '   '.join(line) + '\n'
                        
                        return
                    
    def __set_phase(self):
        """
        Sets orbital phase to all data.
        """
        
        # if there is orbit it will 
        # set phase throughout the data
        if self.__model.has_object('orbit'):
            
            # get the ephemeris from the data
            t0, period, qeph = self.__model['orbit'].get_ephemeris()
            
            # compute phases for all data
            self.__vis2_comparison['phase'] = quad_phase(self.__vis2_comparison['hjd'],
                                                         t0.to('d').value,
                                                         period.to('d').value,
                                                         qeph.value) % 1.0
            self.__cp_comparison['phase'] = quad_phase(self.__cp_comparison['hjd'],
                                                       t0.to('d').value,
                                                       period.to('d').value,
                                                       qeph.value) % 1.0
            self.__lc_comparison['phase'] = quad_phase(self.__lc_comparison['hjd'],
                                                       t0.to('d').value,
                                                       period.to('d').value,
                                                       qeph.value) % 1.0
            self.__vis_comparison['phase'] = quad_phase(self.__vis_comparison['hjd'],
                                                       t0.to('d').value,
                                                       period.to('d').value,
                                                       qeph.value) % 1.0
            self.__sed_comparison['phase'] = quad_phase(self.__sed_comparison['hjd'],
                                                       t0.to('d').value,
                                                       period.to('d').value,
                                                       qeph.value) % 1.0
            self.__spe_comparison['phase'] = quad_phase(self.__spe_comparison['hjd'],
                                                       t0.to('d').value,
                                                       period.to('d').value,
                                                       qeph.value) % 1.0
            self.__has_phase = True

        else:
            raise NotImplementedError('No orbit has been attached, phases have to be '
                                      'computed manually.')

    def __update_fitted_parameters(self, pars, fitparams):
        """
        Assigns a list of fitted parameters to
        list of dictionaries describing the parameters.
        :param pars: a list of values
        :param fitparams: a list of dictionaries generated by self.get_fitted_parameters()
        :return:
        """
        assert len(fitparams) == len(pars)

        for i in range(len(fitparams)):
            objname = fitparams[i]['object']
            parname = fitparams[i]['parameter']
            self.__model[objname][parname] = pars[i]

    def __update_from_orbit(self):
        """
        If the object has an orbit attached -- some
        parameters are updated based on the orbital
        elements.
        :return:
        """

        # if there is no orbit to update, just end
        if not self.__model.has_object('orbit'):
            return
        else:
            self.__model['orbit'].resolve_constraints()  # added by MB, Nov 27th 2018

            K1 = self.__model['orbit'].get_semiamplitude('primary'); print("K1 = ", K1, " km/s")  # dbg
            K2 = self.__model['orbit'].get_semiamplitude('secondary'); print("K2 = ", K2, " km/s")  # dbg

            if self.__model.has_object('central_object'):
                # from orbit get mass
                m = self.__model['orbit'].get_mass('primary').to('solMass')

                # assign it to central object
                self.__model['central_object']['emstar'] = m.value

                K1 = self.__model['orbit'].get_semiamplitude('primary')
                self.__model['central_object']['vyst'] = -K1

                # if istar = 2, update radius from filling factor
                if self.__model['central_object']['istar'] == 2:
                    fill_fact = self.__model['central_object']['ffst']
                    radius = self.__model['orbit'].filling_factor_to_rpole(fill_fact, 'primary')
                    self.__model['central_object']['rstar'] = radius.to('solRad').value

            if self.__model.has_object('companion'):
                # from orbit get semimajor axis and mass ratio
                q = self.__model['orbit']['q']
                sma = self.__model['orbit']['sma'].to('solRad').value

                # assign it to companion
                self.__model['companion']['qq'] = q
                self.__model['companion']['xcp'] = sma

                K2 = self.__model['orbit'].get_semiamplitude('secondary')
                self.__model['companion']['vycp'] = -K2

                # if icomp = 2, update radius from filling factor
                if self.__model['companion']['icomp'] == 2:
                    fill_fact = self.__model['companion']['ffcp']
                    radius = self.__model['orbit'].filling_factor_to_rpole(fill_fact, 'secondary')
                    self.__model['companion']['rcp'] = radius.to('solRad').value

            if self.__model.has_object('envelope'):
                # get mass of the central object mass ratio and sma
                m = self.__model['orbit'].get_mass('primary').to('solMass').value
                q = self.__model['orbit']['q']
                sma = self.__model['orbit']['sma'].to('solRad').value

                # set it to envelope
                self.__model['envelope']['emen'] = m
                self.__model['envelope']['qqen'] = q
                self.__model['envelope']['aen'] = sma

            if self.__model.has_object('disk'):
                # get mass and radius of the central object
                m = self.__model['orbit'].get_mass('primary').to('solMass').value
                r = self.__model['central_object']['rstar'].to('solRad').value
                K1 = self.__model['orbit'].get_semiamplitude('primary')

                # assign it to disk
                self.__model['disk']['emdc'] = m
                self.__model['disk']['rdc'] = r
                self.__model['disk']['vydc'] = -K1

            if self.__model.has_object('nebula'):
                # get mass, radius and RV of the central object
                m = self.__model['orbit'].get_mass('primary').to('solMass').value
                r = self.__model['central_object']['rstar'].to('solRad').value
                K1 = self.__model['orbit'].get_semiamplitude('primary')

                # assign it to nebula
                self.__model['nebula']['emnb'] = m
                self.__model['nebula']['rnb'] = r
                self.__model['nebula']['vynb'] = -K1

            if self.__model.has_object('flow'):
                sma = self.__model['orbit']['sma']
                qq = self.__model['companion']['qq']
                P = self.__model['flow']['pfw']
                omega = 2.*np.pi/P
                r0 = qq/(1.+qq)*sma
                vyfw = (-r0*omega).to('km / s')

                print("sma = ", sma)  # dbg
                print("qq = ", qq)  # dbg
                print("P = ", P)  # dbg
                print("r0 = ", r0)  # dbg
                print("vyfw = ", vyfw)  # dbg

                self.__model['flow']['vxfw'] = 0.0
                self.__model['flow']['vyfw'] = vyfw
                #sys.exit(1)

    def __update_model_from_initial_template(self):
        """
        Updates model parameters based on the template that was
        passed to Interface class.
        :return:
        """

        # go over all model parameters
        for objname in list(self.__model.keys()):
            # orbit is a "speciual case"
            if objname == 'orbit':
                for parname in ['dinc']:
                    # read them from the template and assign
                    value = self.__read_parameter_from_template(parname)
                    dtype = self.__model[objname].get_parameter(parname, 'dtype')
                    self.__model[objname][parname] = dtype(value)
            else:
                for parname in list(self.__model[objname].keys()):

                    # exception for two parameter of the spot, that are not in
                    # the shellspec.in
                    if objname == 'spot' and parname in ['rpolsp', 'vpolsp', 'pangsp']:
                        continue
                    if objname == 'jet' and parname in ['rpoljt', 'vpoljt', 'pangjt']:
                        continue
                    if objname == 'flow' and parname in ['rpolfw', 'pangfw', 'r12fw', 'z12fw', 'v12fw']:
                        continue
                    # read them from template and assign
                    value = self.__read_parameter_from_template(parname)
                    dtype = self.__model[objname].get_parameter(parname, 'dtype')
                    try:
                        self.__model[objname][parname] = dtype(value)
                    except TypeError as ex:
                        print("objname = ", objname)
                        print("parname = ", parname)
                        print("value = ", value)
                        raise ex

    def __str__(self):
        """
        String representation of the class.
        :return:
        """
        string = ""
        string += '==================================================================================================\n'
        string += 'THE MODEL\n'
        string += '==================================================================================================\n'
        # append string of individual objects
        for objname in list(self.__model.keys()):
            string += str(self.__model[objname])

        return string


    def compute_sed(self, dtype='sed'):
        """
        Computes spectral-energy distribution (SED).
        """
        self.set_model_to_shellspec()

        if dtype == 'sed':
            self.__sedsyn['flux'] = []
        else:
            self.__spesyn['flux'] = []

        sedsyn = self.get_ew_and_phase(dtype)

        ncpu = self.__ncpu
        newphase = len(sedsyn['phase'])
        print("newphase = " + str(newphase))

        directories = []
        for i in range(newphase):
            directory = os.path.join(os.getcwd(), 'temp' + dtype + str(i).zfill(2))
            if not os.path.isdir(directory):
                os.mkdir(directory)
            directories.append(directory)

        sedsyn['eff_wavesyn'] = []
        sedsyn['fluxsyn'] = []
        for i in range(0, newphase, ncpu):

            threads = []
            queues = []
            results = []
            nprocess = min([ncpu, newphase - i])
            print("nprocess = " + str(nprocess))

            for j in range(nprocess):
                phase = sedsyn['phase'][i + j]
                ew = sedsyn['eff_wave'][i + j]
                # print "phase = " + str(phase)  # dbg
                # print "ew = " + str(ew)  # dbg

                queues.append(Queue())
                args = (directories[i+j], ew, phase, queues[j], dtype)
                threads.append(Process(target=self.__pp_get_sed_one_spectrum, args=args))

            for j in range(nprocess):
                threads[j].start()

            for j in range(nprocess):
                results.append(queues[j].get())
                queues[j].close()
                threads[j].join()
                if threads[j].exitcode != 0:
                    print("compute_sed: Error running shellspec! exitcode = " + str(threads[j].exitcode))
                    sys.exit(1)

            for j in range(nprocess):
                wave, flux = results[j]
                sedsyn['eff_wavesyn'].append(wave)
                sedsyn['fluxsyn'].append(flux)

        self.__get_sed_observables(sedsyn, dtype=dtype)
        if dtype == 'sed':
            self.__sedsyn = copy.deepcopy(sedsyn)
        else:
            self.__spesyn = copy.deepcopy(sedsyn)

    def __pp_get_sed_one_spectrum(self, directory, ew, phase, queue=None, dtype='sed'):
        """
        Computes one spectrum.
        :param directory:
        :param ew:
        :param phase:
        :param queue:
        :return:
        """

        # set wavelength range
        ang = 1.e-10
        w0 = min(ew)/ang
        wn = max(ew)/ang
        step = (wn-w0)/len(ew)
        self.set_wavelength(w0=w0, wn=wn, step=step)
        self.__set_parameter_to_template('loglam', 2)  # we compute only required lambdas
        lambda_file = os.path.join(directory, 'lambda')
        self.write_phase(np.array(ew)/ang, filename=lambda_file)

        # set limb-darkening coefficients
        w1 = 0.5*(min(ew)+max(ew))
        dlst = self.limcofst.interp(w1)
        dlcp = self.limcofcp.interp(w1)
        self.__set_parameter_to_template('dlst', dlst)
        self.__set_parameter_to_template('dlcp', dlcp)

        phase_file = os.path.join(directory, 'phases')
        control_file = os.path.join(directory, 'shellspec.in')
        self.write_phase([phase], filename=phase_file)
        self.write_template(filename=control_file)

        cwd = os.getcwd()
        os.chdir(directory)
        exitcode = self.__run_shellspec()
        os.chdir(cwd)

        if exitcode != 0:
            if queue is not None:
                queue.put(None)
            sys.exit(exitcode)

        wave, flux = self.__read_spectrum(filename=os.path.join(directory, 'shellspectrum'), dtype=dtype)

        if queue is not None:
            queue.put([wave, flux])
        else:
            return wave, flux


    def __read_spectrum(self, filename='shellspectrum', dtype='sed'):
        """
        Reads spectrum by shellspec:

        0 ... lambda [Ang]
        1 ... radial velocity [km s^-1]
        2 ... monochromatic flux F_nu [erg s^-1 cm^-2 Hz^-1]
        3 ... F_lambda [erg s^-1 cm^-2 cm^-1]
        4 ... normalized F_nu []
        5 ... dtto F_nu [] shifted for subsequent phases by 'offset'

        """

        ifile = open(filename, 'r')
        lines = ifile.readlines()
        ifile.close()

        wave = []
        flux = []
        if dtype == 'sed':
            j = 3
        else:
            j = 5
        for l in lines:
            d = l.split()
            if len(d) < 5:
                continue
            else:
                wave.append(float(d[0]))
                flux.append(float(d[j]))

        wave = np.array(wave)
        flux = np.array(flux)

        ang = 1.e-10
        wave = wave * ang
        if dtype == 'sed':
            erg = 1.e-7
            cm = 1.e-2
            flux = flux * erg*cm**(-2)*cm**(-1)

        return wave, flux

    def __get_sed_observables(self, sedsyn, dtype='sed'):
        """
        Evaluates synthetic spectrum at desired wavelength points.
        :return:
        """

        for i, phase in enumerate(sedsyn['phase']):

            synwave = sedsyn['eff_wavesyn'][i]
            synflux = sedsyn['fluxsyn'][i]

            if dtype == 'sed':
                ph = self.__sed_comparison['phase']
                phase_precision = self.__sed_phase_precision
            else:
                ph = self.__spe_comparison['phase']
                phase_precision = self.__spe_phase_precision

            # Note: this has to be slightly more precise than precision
            # (otherwise, we would obtain two phases for one observed)
            diff = np.absolute(ph - phase)
            idx = np.where((diff <= 0.5*10**(-phase_precision)))[0]

            if dtype == 'sed':
                obswave = self.__sed_comparison['eff_wave'][idx]
                obsflux = self.__sed_comparison['flux'][idx]
            else:
                obswave = self.__spe_comparison['eff_wave'][idx]
                obsflux = self.__spe_comparison['flux'][idx]

            synflux_itp = interpolate1d_hermite(obswave, synwave, synflux)

            if dtype == 'sed':
                self.__sed_comparison['fluxsyn'][idx] = synflux_itp
            else:
                self.__spe_comparison['fluxsyn'][idx] = synflux_itp


    def compute_spe(self):
        """
        Computes synthetic spectrum (normalized).
        """
        self.compute_sed(dtype='spe')


