#!/usr/bin/env python

import astropy.constants
import astropy.coordinates
import astropy.time
import astropy.units as units
import copy
import numpy as np
import re

from .auxilliary import load_ascii
from .definitions import filters, eff_wave, eff_band, calibration_flux
from .definitions import location_definitions
from .include.oifits import oifits


class Data:
    """
    This class wraps several kinds of observations.
    """

    def __init__(self, obs):
        """
        """

        # intialize structures containing the data
        self.__lc_observations = []
        self.__if_observations = []
        self.__df_observations = []
        self.__sed_observations = []
        self.__spe_observations = []

        # adds the passed observations
        self.add_observation(obs)

    def add_observation(self, obs):
        """
        Add observations to either light curves, or interferometry.
        :param obs: an observation wrapped either in LCData or IFData
        """

        # in it is passed in form of a list
        if isinstance(obs, list):
            for i in range(len(obs)):
                if isinstance(obs[i], LCData):
                    self.__lc_observations.append(obs[i])
                elif isinstance(obs[i], IFData):
                    self.__if_observations.append(obs[i])
                elif isinstance(obs[i], DFData):
                    self.__df_observations.append(obs[i])
                elif type(obs[i]) == SEDData:
                    self.__sed_observations.append(obs[i])
                elif type(obs[i]) == SPEData:
                    self.__spe_observations.append(obs[i])

        # and the remaining options
        elif isinstance(obs, LCData):
            self.__lc_observations.append(obs)
        elif isinstance(obs, IFData):
            self.__if_observations.append(obs)
        elif isinstance(obs, DFData):
            self.__df_observations.append(obs)
        elif type(obs) == SEDData:
            self.__sed_observations.append(obs)
        elif type(obs) == SPEData:
            self.__spe_observations.append(obs)

    def get_filenames(self, dtype):
        """
        Returns list of all filenames.
        :param dtype:
        :return:
        """
        filenames = []
        for i in range(0, len(self.get_observation_number(dtype))):
            filename = self.get_observation(dtype, i).get_filename()
            filenames.append(filename)

        return filenames

    def get_observation(self, dtype, idx):
        """
        Returns idx-th observation of the dtype.
        :param dtype:
        :param idx:
        :return:
        """
        if dtype == 'lc':
            return self.__lc_observations[idx]
        elif dtype == 'if':
            return self.__if_observations[idx]
        elif dtype == 'df':
            return self.__df_observations[idx]
        elif dtype == 'sed':
            return self.__sed_observations[idx]
        elif dtype == 'spe':
            return self.__spe_observations[idx]

    def get_observation_number(self, dtype):
        """
        Returns number of registered data for selected type.
        :param dtype: either 'lc' or 'if'
        :return:
        """
        if dtype == 'lc':
            return len(self.__lc_observations)
        elif dtype == 'if':
            return len(self.__if_observations)
        elif dtype == 'df':
            return len(self.__df_observations)
        elif dtype == 'sed':
            return len(self.__sed_observations)
        elif dtype == 'spe':
            return len(self.__spe_observations)


class IFData(object):
    """
    Container for the interferometric data. Relies hugely on 
    oifits library by Paul Boley --- http://astro.ins.urfu.ru/pages/~pboley/oifits/.
    """

    def __init__(self, filename, dec=None, ra=None, location=None, format='FITS', exclude_cp=False, exclude_t3amp=False, weight_vis2=1.0, weight_t3amp=1.0):
        """
        
        :param filename: 
        :param dec: declination in deg
        :param ra: right ascension in deg
        :param location: name of a location predefined in definitions.py
        :param format: input format, either FITS or ASCII
        :param exclude_cp: excludes closure phase from reading (for ASCII)
        :param exclude_t3amp: excludes triple product amplitude from reading (dtto)
        :param weight_vis2: optional weigth for |V^2| data
        :param weight_t3amp: weigth for |T_3| data
        """

        # initialization of arguments
        self.__dec = dec
        self.__filename = filename
        self.__location = location
        self.__ra = ra
        self.weight_vis2 = weight_vis2
        self.weight_t3amp = weight_t3amp

        # initialization of triggers
        self.__has_closure_phase = False
        self.__has_visibility = False
        self.__has_phase = False
        self.__has_visibility1 = False

        # initialization of data-holding structures
        self.array = None
        # self.phase = None
        self.target = None
        self.t3 = None
        self.vis2 = None
        self.wavelength = None

        # convert declination and right ascension tu units
        if dec is not None:
            self.__dec *= units.deg
        if ra is not None:
            self.__ra *= units.deg

        # if the file is in fits
        if format.lower() == 'fits':
            # loads the data
            self.load_data()

            # if location was passed --- it is either taken by name
            # or one has to pass the two locations. There is
            # third option, that it will be read from oifits.
            if self.__location is not None:
                self.__set_location()

            # transform utc --> HJD and assign effective
            # wavelength
            if self.__has_visibility:
                self.__utc_to_hjd(self.vis2)
                self.__set_effective_wavelength(self.vis2)

            if self.__has_closure_phase:
                self.__utc_to_hjd(self.t3)
                self.__set_effective_wavelength(self.t3)

            if self.__has_visibility1:
                self.__utc_to_hjd(self.vis)
                self.__set_effective_wavelength(self.vis)
                self.__barycentric_correction(self.vis)

            # data has to be read from structures
            self.__data = None

        # if the file is in ascii format
        elif format.lower() == 'ascii':
            self.__data = self.load_data_ascii(exclude_cp, exclude_t3amp)
        else:
            raise ValueError('Only eligible formats are \'ascii\' and \'fits\'.')

    def fits_to_ascii(self, filename=None):
        """
        If the file is in fits format it is written
        into ascii in a format that is readable for this class.
        :return:
        """
        if filename is None:
            def_filename = True
        else:
            def_filename = False

        # first extract visibility
        if self.__has_visibility and self.vis2 is not None:
            if def_filename:
                filename = '.'.join(self.__filename.split('.')[:-1]) + '.vis2.dat'

            # get data
            data = self.get_data(dtype='vis2')
            attributes = ['hjd', 'ucoord', 'vcoord', 'vis2data', 'vis2err', 'eff_wave', 'eff_band']

            # get header
            header = ''.join(['%20s' % (attr) for attr in attributes])

            # get data in a list
            data = [data[attr] for attr in attributes]

            # write it
            np.savetxt(filename, np.column_stack(data), fmt='%20.12e', header=header)

        # second extract closure phase
        if self.__has_closure_phase and self.t3 is not None:
            if def_filename:
                filename = '.'.join(self.__filename.split('.')[:-1]) + '.cp.dat'

            # get data
            data = self.get_data(dtype='cp')
            attributes = ['hjd', 'u1coord', 'v1coord', 'u2coord', 'v2coord', 't3amp', 't3amperr', 't3phi', 't3phierr',
                          'eff_wave', 'eff_band']

            # get header
            header = ''.join(['%20s' % (attr) for attr in attributes])

            # get data in a list
            data = [data[attr] for attr in attributes]

            # write it
            np.savetxt(filename, np.column_stack(data), fmt='%20.12e', header=header)

        # thrid extract visibility (not squared)
        if self.__has_visibility1 and self.vis is not None:
            if def_filename:
                filename = '.'.join(self.__filename.split('.')[:-1]) + '.vis.dat'

            # get data
            data = self.get_data(dtype='vis')
            attributes = ['hjd', 'ucoord', 'vcoord', 'visamp', 'visamperr', 'visphi', 'visphierr', 'eff_wave', 'eff_band']

            # get header
            header = ''.join(['%20s' % (attr) for attr in attributes])

            # get data in a list
            data = [data[attr] for attr in attributes]

            # write it
            np.savetxt(filename, np.column_stack(data), fmt='%20.12e', header=header)

    def get_data(self, dtype='vis2'):
        """
        Returns all observed data in the file.
        :param dtype: 
        """

        if dtype == 'vis2':
            if self.__has_visibility == False:
                return
            attributes = ['hjd', 'ucoord', 'vcoord', 'vis2data', 'vis2err', 'eff_wave', 'eff_band', 'flag']
        elif dtype == 'cp':
            if self.__has_closure_phase == False:
                return
            attributes = ['hjd', 'u1coord', 'v1coord', 'u2coord', 'v2coord', 't3amp', 't3amperr', 't3phi', 't3phierr',
                          'eff_wave', 'eff_band', 'flag']
        elif dtype == 'vis':
            if self.__has_visibility1 == False:
                return
            attributes = ['hjd', 'ucoord', 'vcoord', 'visamp', 'visamperr', 'visphi', 'visphierr', 'eff_wave', 'eff_band', 'flag']
        else:
            raise ValueError('dtype not understood --- allowed values are \'vis2\', \'cp\' and \'vis\'.')

        # if the wrapped filename was ascii its contents are just returned
        if self.__data is not None:
            return self.__data

        # builds a dictionary 
        data = {attr: self.get_property(attr, dtype) for attr in attributes}

        # finalizes the output data set
        data = self.__verify_data(data, dtype)

        # transform into arrays
        data = {k: np.array(data[k]) for k in list(data.keys())}

        return data


    def get_filename(self):
        """
        Returns the filenam.
        :return:
        """

        return self.__filename

    def get_property(self, attribute, dtype='vis2'):
        """
        Returns an attribute.
        :param attribute:
        :param dtype: 
        """

        if dtype == 'vis2':
            obj = self.vis2
        elif dtype == 'cp':
            obj = self.t3
        elif dtype == 'vis':
            obj = self.vis
        else:
            raise ValueError('dtype not understood --- allowed values are \'vis2\', \'cp\' and \'vis\'.')

        # just take the property from the array
        try:
            prop = [getattr(obj[i], attribute.lower()) for i in range(len(obj))]
        except AttributeError:
            prop = [getattr(obj[i], attribute.upper()) for i in range(len(obj))]
        return prop

    def load_data(self):
        """
        Loads the oifits from a fits file using the 
        oifits package.
        """

        # opens the structure
        hdu = oifits.open(self.__filename)

        # target
        self.target = hdu.target

        # array
        self.array = hdu.array

        # wavelengths
        self.wavelength = hdu.wavelength
#        print self.wavelength
#        print self.wavelength.keys()
#        print self.wavelength['VEGA']
#        print self.wavelength['VEGA'].eff_wave
#        print self.wavelength['VEGA'].eff_band

        # readout data - squared visibility
        if len(hdu.vis2) > 0:
            self.vis2 = hdu.vis2
            self.__has_visibility = True

        # closure phase
        if len(hdu.t3) > 0:
            self.t3 = hdu.t3
            self.__has_closure_phase = True

        # visibility amplitude and phase
        if len(hdu.vis) > 0:
            self.vis = hdu.vis
            self.__has_visibility1 = True

    def load_data_ascii(self, exclude_cp=False, exclude_t3amp=False):
        """
        Loads data in ascii format with each column
        preceeded by a column name
        :return:
        """
        # structure for the output
        data = {}

        # open and read the file
        ifile = open(self.__filename, 'r')
        lines = ifile.readlines()
        ifile.close()

        # read the header
        header = lines[0].split()
        if header[0].find('#') > -1:
            if len(header[0]) == 1:
                header = header[1:]
            else:
                header[0] = header[0].lstrip('#')
        header = [k.lower() for k in header]

        # determine data type
        vis2_attrs = ['hjd', 'ucoord', 'vcoord', 'vis2data', 'vis2err', 'eff_wave', 'eff_band']
        cp_attrs = ['hjd', 'u1coord', 'v1coord', 'u2coord', 'v2coord', 't3amp', 't3amperr', 't3phi', 't3phierr', 'eff_wave', 'eff_band']
        vis_attrs = ['hjd', 'ucoord', 'vcoord', 'visamp', 'visamperr', 'visphi', 'visphierr', 'eff_wave', 'eff_band']

        vis2count = 0
        cpcount = 0
        viscount = 0
        for k in header:
            if k in vis2_attrs:
                vis2count += 1
            if k in cp_attrs:
                cpcount += 1
            if k in vis_attrs:
                viscount += 1
            if k not in vis2_attrs + cp_attrs + vis_attrs:
                raise KeyError('Column %s in file %s was not understood.' % (k, self.__filename))

        # recognise the data type
        if viscount > vis2count and viscount > cpcount:
            self.__has_visibility = False
            self.__has_closure_phase = False
            self.__has_visibility1 = True
        elif vis2count > cpcount:
            self.__has_visibility = True
            self.__has_closure_phase = False
            self.__has_visibility1 = False
        else:
            self.__has_visibility = False
            self.__has_closure_phase = True
            self.__has_visibility1 = False

        # start the lists for each record
        data = {k.lower(): [] for k in header}

        # read out the data
        for line in lines[1:]:
            if not re.match("^#", line) and len(line) > 1:
                temp = list(map(float, line.split()))
                for i, k in enumerate(header):
                    if ((exclude_t3amp and (k == 't3amp' or k == 't3amperr'))
                        or (exclude_cp and (k == 't3phi' or k == 't3phierr'))):
                        data[k].append(None)
                    else:
                        data[k].append(temp[i])
        return data


    def set_filename(self, f):
        """
        Assigns filename to the class
        :param f: new filename
        :return:
        """
        self.__filename = f

    def __set_effective_wavelength(self, obj):
        """
        Assigns effective wavelength to each data record...
        it is kinda redundant though, since the data are 
        there...
        """

        for i in range(len(obj)):
            setattr(obj[i], 'eff_wave', obj[i].wavelength.eff_wave)
            setattr(obj[i], 'eff_band', obj[i].wavelength.eff_band)

    def __set_location(self):
        """
        Sets location of the observatory
        """
        if isinstance(self.__location, str):
            idx = location_definitions[0].index(self.__location)
            self.__latitude = location_definitions[1][idx] * units.deg
            self.__longitude = location_definitions[2][idx] * units.deg
        elif isinstance(self.__location, (tuple, list)):
            self.__latitude = self.__location[0]
            self.__longitude = self.__location[1]

    def __utc_to_hjd(self, obj):
        """
        For given object, it transforms the universal time 
        to heliocentric Julian date. 
        
        :param obj: oifits object - either t3, or vis2
        """
        # coordinates of the star on the sky

        # extract coordinates of the target
        if self.__ra == None or self.__dec == None:
            if len(self.target) > 1:
                raise ValueError(
                    'The fits file contains data for more targets. Currently only one target is allowed! Sry...')
            else:
                # extract the only target
                target = self.target[0]

                # extract coordinates on the sky
                ra = target.raep0
                dec = target.decep0

                # just a test that there aren't any wrong values...
                if abs(ra.angle) < 1e-6 and abs(dec.angle) < 1e-6:
                    raise ValueError(
                        'Coordinates could have not been read from the fits file %s. '
                        'They have to be supplied by the user.' % (self.__filename))
                else:
                    self.__ra = ra.angle * units.deg
                    self.__dec = dec.angle * units.deg

        # sets the coordinates for the computation of Julian date
        # on the sky and on the Earth
        skypos = astropy.coordinates.SkyCoord(self.__ra, self.__dec, unit=(units.deg, units.deg), frame='icrs')
        earloc = astropy.coordinates.EarthLocation(lat=self.__latitude, lon=self.__longitude)

        for i in range(len(obj)):
            # computes Julian date
            time = astropy.time.Time(obj[i].timeobs.isoformat(' '), format='iso', scale='utc')

            # computes the heliocentric correction
            helcor = time.light_travel_time(skypos, 'heliocentric', earloc)

            # computes the heliocentric times
            timecor = time.utc + helcor

            # converts to HJD
            hjd = timecor.jd

            # append hjd to the structure
            setattr(obj[i], 'hjd', hjd)

    def __barycentric_correction(self, obj):
        """Apply barycentric RV correction to differential IF data"""

        skypos = astropy.coordinates.SkyCoord(self.__ra, self.__dec, unit=(units.deg, units.deg), frame='icrs')
        earloc = astropy.coordinates.EarthLocation(lat=self.__latitude, lon=self.__longitude)

        for i in range(len(obj)):
            time = astropy.time.Time(obj[i].timeobs.isoformat(' '), format='iso', scale='utc')
            barcor = skypos.radial_velocity_correction(kind='barycentric', obstime=time, location=earloc)

            print("time, barcor = ", time.jd, barcor)  # dbg

            tmp = barcor / astropy.constants.c
            for j in range(len(obj[i].eff_wave)):
                obj[i].eff_wave[j] += obj[i].eff_wave[j] * tmp

    def __verify_data(self, data, dtype):
        """
        Removes flagged data and transfers them into 1-column-per-data-type format.
        
        :param data:
        """
        # records where one value applies to a set of visibilities
        common_dtypes = ['hjd', 'ucoord', 'vcoord', 'u1coord', 'v1coord', 'u2coord', 'v2coord']

        # extract keys defining the data
        keys = list(data.keys())

        # set the reference data type
        if dtype == 'vis2':
            refdt = 'vis2data'
        elif dtype == 'cp':
            refdt = 't3phi'
        elif dtype == 'vis':
            refdt = 'visamp'

        # output structure -- remove dtype flag -- all data
        # are correct in the output structure
        output_data = {k: [] for k in keys}
        output_data.pop('flag')


        # go over each row
        for i in range(len(data[refdt])):

            # remove flagged data 
            idx_not_flagged = np.where(data['flag'][i] == False)[0]

            for k in keys:
                if k in common_dtypes:
                    continue
                else:
                    # if isinstance(data[k][i], list):
                    if data != None:
                        if len(data[k][i]) > 1:
                            data[k][i] = [data[k][i][iok] for iok in idx_not_flagged]
                        else:
                            if not idx_not_flagged[0]:
                                data[k][i] = [data[k][i]]
                            else:
                                data[k][i] = []

            # extend lists that are common for a given row of data
            if isinstance(data[refdt][i], list):
                ref_data_length = len(idx_not_flagged)
                for k in common_dtypes:
                    if k in list(data.keys()):
                        data[k][i] = [data[k][i] for rdl in range(ref_data_length)]

            # append data into output structure - where there is one
            # column for each data type
            for k in list(output_data.keys()):
                output_data[k].extend(copy.deepcopy(data[k][i]))

        return output_data


class DFData(IFData):
    """Differential visibility data."""


class LCData(object):
    """
    Container for the photometric in one pass-band. 
    Each set of observations is given by a passband, 
    filename and the set of corresponding observations.
    """

    def __init__(self, filename, passband='johnson_v', columns=['hjd', 'magnitude', 'error'], global_error=None,
                 comment='#', delimiter=None):
        """
        Constructs the class.
       
        :param filename: Ascii file containing columns defined in 'columns'.
        :param passband: Passband for which teh data were acquired.
        :param columns: Name of columns.
        :param global_error: error that is the same for all records in the file.
        :param comment:
        :param delimiter:
        """

        # Save the arguments
        self.__filename = filename
        self.__passband = passband.lower()
        self.__columns = [c.lower() for c in columns]
        self.__global_error = global_error

        # checks columns names
        self.__check_columns()

        # empty arrays for the data
        self.hjd = None
        self.magnitude = None
        self.error = None
        self.eff_wave = None
        self.eff_band = None
        self.calibration_flux = None

        # set the effective wavelength
        self.__set_effective_wavelength()

        # load the data
        self.load_data(**dict(delimiter=delimiter, comment=comment))

    def get_data(self):
        """
        Returns all data stored within the object
        """
        # empty dictionary for the data
        data = {}

        # append data
        data['eff_wave'] = np.ones(self.hjd.size) * self.eff_wave
        data['eff_band'] = np.ones(self.hjd.size) * self.eff_band
        data['error'] = self.error
        data['hjd'] = self.hjd
        data['magnitude'] = self.magnitude

        return data

    def get_effective_wavelength(self):
        """
        Returns the effective wavelength in nanometers.
        """

        return self.eff_wave

    def get_filename(self):
        """
        Returns the filenam.
        :return:
        """

        return self.__filename

    def get_passband(self):
        """
        Returns the passband.
        """
        return self.__passband

    def load_data(self, **kwargs):
        """
        Loads the photometric observations.
        """

        # load columns designation
        columns = copy.deepcopy(self.__columns)

        # assign global error if local unavailable
        if self.__global_error is not None and 'errors' in columns:
            columns.remove('error')

        # load the data
        data = load_ascii(self.__filename, columns, **kwargs)

        for i, c in enumerate(columns):
            # in case we wanted to skip on a column
            if c is None:
                continue
            else:
                setattr(self, c, data[:, i])

        # assign global error
        if self.__global_error is not None:
            setattr(self, 'error', self.__global_error * np.ones(np.shape(data)[0]))

    def __check_columns(self):
        """
        Checks columns definition.
        """
        for c in self.__columns:
            if c not in ['hjd', 'magnitude', 'error', None]:
                raise ValueError('%c is not among allowed columns names (hjd, magnitude, error).' % c)

    def __set_effective_wavelength(self):
        """
        Assigns an effective wavelength to the data.
        """
        idx = filters.index(self.__passband)
        self.eff_wave = eff_wave[idx] * units.meter
        self.eff_band = eff_band[idx] * units.meter
        self.calibration_flux = calibration_flux[idx]



class SEDData(object):
    """Spectral-energy distribution (SED) data"""

    def __init__(self, filename, columns=['hjd', 'eff_wave', 'eff_band', 'flux', 'error', 'dataset'], global_error=None):
        self.__filename = filename
        self.__columns = [c.lower() for c in columns]
        self.__global_error = global_error

        self.hjd = None
        self.eff_wave = None
        self.eff_band = None
        self.flux = None
        self.error = None
        self.dataset = None

        self.load_data(delimiter=None, comment='#')

    def load_data(self, **kwargs):
        columns = copy.deepcopy(self.__columns)

        if self.__global_error is not None and 'errors' in columns:
            columns.remove('error')

        data = load_ascii(self.__filename, columns, **kwargs)

        for i, c in enumerate(columns):
            if c is None:
                continue
            else:
                setattr(self, c, data[:, i])

        if self.__global_error is not None:
            setattr(self, 'error', self.__global_error * np.ones(np.shape(data)[0]))

    def get_data(self):
        data = {}

        data['hjd'] = self.hjd
        data['eff_wave'] = self.eff_wave
        data['eff_band'] = self.eff_band
        data['flux'] = self.flux
        data['error'] = self.error
        data['dataset'] = self.dataset

        return data

    def get_filename(self):
        return self.__filename


class SPEData(SEDData):
    """Spectral (normalized) data"""


