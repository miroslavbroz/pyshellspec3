#!/usr/bin/env python3

import numpy as np
from astropy import constants as cst
from astropy import units

from .auxilliary import filling_factor_to_rpoint
from .auxilliary import l1_position
from .auxilliary import potential_point
from .auxilliary import polar_radius
from .auxilliary import semiamplitude
from .auxilliary import radial_velocity
from .definitions import object_definitions
from .parameter import Parameter


class BaseObject(object):
    """
    """
    def __init__(self, name, parameters=[], **kwargs):
        """
        Initialize the base class.
        """

        # pass command line arguments
        self.__name = name
        self.__parameters = parameters

        # if the name is within object_definitions, 
        # then the class attempts to define parameters 
        # from ''object_definitions''
        if self.__name in list(object_definitions.keys()) and len(parameters) == 0:
            for pdict in object_definitions[self.__name]:
                self.__parameters.append(Parameter(**pdict))
                
        # save defined parameters
        self.__defined_parameters = [p.get_property('name') for p in self.__parameters]
                
        # (optional) sets initial values of each parameter
        self.__set_initial_values(**kwargs)

    def check_boundaries(self):
        """
        Checks that no parametr lies outside boundaries
        :return:
        """
        for p in self.__parameters:
            p.check_boundaries()

    def get_name(self):
        """
        Returns name of the object.
        """
        
        return self.__name
        
    def get_parameter(self, parname, attr='value'):
        """
        Returns an attribute of a parameter.
        :param parname:
        :param attr:
        """
        
        for p in self.__parameters:
            if p.get_property('name') == parname:
                return p.get_property(attr)

    def keys(self):
        """
        Returns a list of names of all defined parameters.
        :return:
        """

        return self.__defined_parameters

    def resolve_constraints(self):
        """
        Resolves constraints separately for each object.
        :return:
        """
        return
        
    def set_parameter(self, parname, **kwargs):
        """
        Sets properties of parameters --- all that are defined in kwargs.
        :param parname:
        :param kwargs:
        """
    
        # make it case insensitive
        parname = parname.lower()
        
        # check that parameter is defined
        if parname not in self.__defined_parameters:
            raise ValueError('Parameter %s is not defined for the object %s.' % (parname, self.__name))

        # assign the parameters
        for i in range(0, len(self.__parameters)):
            if self.__parameters[i].get_property('name') == parname:
                for attr in list(kwargs.keys()):
                    self.__parameters[i].set_property(attr, kwargs[attr])

    def __getitem__(self, key):
        """
        Retrieve parameter value through get_parameter method.
        :param key:
        :return:
        """
        return self.get_parameter(key)

    def __setitem__(self, key, value):
        """
        Sets properties through set_parameter method.
        :param key:
        :param value:
        :return:
        """
        self.set_parameter(key, value=value)
        
    def __str__(self):
        """
        String representation of the class.
        """
        string = ''
        string += "==================================================================================================\n"
        string += "Object: %s\n" % self.__name
        for p in self.__parameters:
            string += str(p) + '\n'
        string += "==================================================================================================\n"
            
        return string
    
    def __set_initial_values(self, **kwargs):
        """
        Sets parameters from an initial set of values.
        """
        
        for k in list(kwargs.keys()):
            self.set_parameter(k, value=kwargs[k])


class CentralObject(BaseObject):
    """
    Class wrapping the central object of shellspec.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        :return:
        """

        # initialize the parent class
        super(CentralObject, self).__init__(name='central_object', parameters=[], **kwargs)


class Companion(BaseObject):
    """
    Class wrapping the companion of shellspec.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        :return:
        """

        # initialize the parent class
        super(Companion, self).__init__(name='companion', parameters=[], **kwargs)


class Disk(BaseObject):
    """
    Class wrapping the disk of shellspec.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        :return:
        """

        # initialize the parent class
        super(Disk, self).__init__(name='disk', parameters=[], **kwargs)


class Nebula(BaseObject):
    """
    Class wrapping the nebula of shellspec.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        :return:
        """

        # initialize the parent class
        super(Nebula, self).__init__(name='nebula', parameters=[], **kwargs)


class Distance(BaseObject):
    """
    Class wrapping distance of the system.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        :return:
        """

        # initialize the parent class
        super(Distance, self).__init__(name='distance', parameters=[], **kwargs)


class Environmental(BaseObject):
    """
    Class wrapping the companion of shellspec.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        :return:
        """

        # initialize the parent class
        super(Environmental, self).__init__(name='environmental', parameters=[], **kwargs)

class Envelope(BaseObject):
    """
    Class wrapping the envelope of shellspec.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        :return:
        """

        # initialize the parent class
        super(Envelope, self).__init__(name='envelope', parameters=[], **kwargs)


class Grid(BaseObject):
    """
    Class wrapping properties of the grid. Mainly
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        :return:
        """

        # initialize the parent class
        super(Grid, self).__init__(name='grid', parameters=[], **kwargs)

    def get_pixel_size(self, grid='los'):
        """
        Returns number of pixels along each axis.
        :param grid:
        :return:
        """
        if grid.lower() == 'los':

            # get the lin-of-sight grid size
            x1 = self['rmdx1'].value
            x2 = self['rmdx2'].value
            y1 = self['rmdy1'].value
            y2 = self['rmdy2'].value
            z1 = self['rmdz1'].value
            z2 = self['rmdz2'].value
            stepxy = self['stepx'].value
            stepz = self['stepz'].value

            # print x1, x2, y1, y2, z1, z2, stepxy, stepz

            # compute the pixel size
            # npx = int((x2 - x1) / stepxy) + 1
            # npy = int((y2 - y1) / stepxy) + 1
            # npz = int((z2 - z1) / stepz) + 1

            npx = int(np.around((x2 - x1) / stepxy))
            npy = int(np.around((y2 - y1) / stepxy))
            npz = int(np.around((z2 - z1) / stepz))

            # print npx, npy

            if npx % 2 == 0:
                npx += 1
            if npy  % 2 == 0:
                npy += 1
            if npz % 2 == 0:
                npz += 1
            # print npx, npy, npz

        # raise an error if anything else is querried
        else:
            raise ValueError('Only pixel size of the line-of-sight (grid=\'los\') grid can be querried now.')

        return npx, npy, npz


class Orbit(BaseObject):
    """
    Class wrapping the orbital elements of the object
    if it is a binary.
    """
    def __init__(self, **kwargs):
        """
        A wrapper of the orbit
        """
        
        # initialize the parent class
        # BaseObject.__init__(self, 'orbit', **kwargs)
        super(Orbit, self).__init__(name='orbit', parameters=[], **kwargs)

    def filling_factor_to_rpole(self, ffact, component='primary'):
        """
        Computes polar radius for a given filling factor
        and component.
        :param ffact: the filling factor
        :param component: component of a binary
        :return: rpole polar radius in Rsol
        """
        # get orbital properties
        sma = self['sma']
        q = self['q'].value
        l1 = self.get_l1_position()

        # compute rpoint
        rpoint = filling_factor_to_rpoint(ffact, l1, component=component)

        # compute Kopal potential
        pot = potential_point(rpoint, q, component=component)

        # compute polar radius
        rpole = polar_radius(pot, q, component=component)

        # convert to absolute units
        return rpole * sma

    def get_barycentre(self):
        """
        Returns position of the barycentre with respect
        to the central object in the body-frozen coordinates.
        :return:
        """
        # get the mass ratio and the semimajor axis
        q = self.get_parameter('q')
        sma = self.get_parameter('sma')

        # compute the barycentre position
        bar_pos = q * sma / (1 + q)

        return bar_pos

    def get_ephemeris(self):
        """
        Returns ephemeris, wrapped within units 
        class.
        """
        t0 = self.get_parameter('t0')
        period = self.get_parameter('period')
        qeph = self.get_parameter('qeph')
        
        return t0, period, qeph

    def get_l1_position(self):
        """
        Returns the position of Lagrange1 point.
        :return:
        """
        q = self['q'].value
        return l1_position(q)

    def get_mass(self, component='primary'):
        """
        Returns mass of a component in kg
        :param component:
        :return:
        """

        # orbital properties
        G = cst.G.value
        sma = self['sma'].to('m').value
        period = self['period'].to('s').value
        q = self['q'].value

        # total systemic mass
        M = (sma ** 3 * 4 * np.pi ** 2 / (G * period ** 2)) * units.kg

        if component.lower() == 'primary':
            return M / (q + 1.)
        elif component.lower() == 'secondary':
            return M * q / (q + 1.)
        else:
            raise ValueError('Component can be either \'primary\' or \'secondary\'.')

    def get_semiamplitude(self, component='primary'):
        """
        Computes radial velocity semiamplitude.
        :return:
        """
        # extract all necessary orbital properties
        sma = self['sma'].to('m').value
        ecc = self['ecc'].value
        q = self['q'].value
        incl = self['dinc'].to('rad').value
        period = self['period'].to('s').value

        K = semiamplitude(ecc, q, sma, incl, period, component)

        return K.to('km / s').value

    def get_radial_velocity(self, t, component='primary'):
        """
        Computes radial velocity at given epoch t.
        :param t:
        :return:
        """
        # extract all necessary orbital properties
        t  = t * units.d
        t = t.to('s').value
        sma = self['sma'].to('m').value
        ecc = self['ecc'].value
        q = self['q'].value
        incl = self['dinc'].to('rad').value
        omega = self['omega'].to('rad').value
        t0 = self['t0'].to('s').value
        period = self['period'].to('s').value
        qeph = self['qeph'].value
        dperiod = 2 * qeph / period

        RV = radial_velocity(t, ecc, q, sma, incl, omega, period, dperiod, t0, component)

        return RV.to('km / s').value

    def resolve_constraints(self):
        """If a sin(i) is set, compute semimajor axis sma for given inclination dinc."""

        asini = self.get_parameter('asini')
        if asini > 0.0:
            self['sma'] = asini / np.sin(self.get_parameter('dinc').to('rad'))
        print("asini = ", asini)
        print("sini = ", self.get_parameter('dinc'))
        print("sma = ", self.get_parameter('sma'))

class Spot(BaseObject):
    """
    Wraps the object spot of Shellspec.
    """
    def __init__(self, **kwargs):
        """
        A wrapper of the spot
        """

        # initialize the parent class
        super(Spot, self).__init__(name='spot', parameters=[], **kwargs)

    def resolve_constraints(self):
        """
        If the polar coordiantes of the spot
        are fitted. This function propagates
        any changes to the xsp and ysp
        :return:
        """
        # propagates polar coordinates into Cartesian
        # if they are fitted
        # if self.get_parameter('rpolsp', 'fitted') or self.get_parameter('pangsp', 'fitted'):
        if self.get_parameter('rpolsp', 'value') > 0.0:
            self.set_position(r=self.get_parameter('rpolsp'), angle=self.get_parameter('pangsp'))

        if abs(self.get_parameter('vpolsp', 'value')) > 0.0:
            self.set_velocity(v=self.get_parameter('vpolsp'), angle=self.get_parameter('pangsp'))

    def set_position(self, x=0.0, y=0.0, r=0.0, angle=0.0, degrees=True):
        """
        Sets the position of the spot
        :param x: Cartesian coordinate in solRad
        :param y: Cartesian coordinate in solRad
        :param r: polar coordinate in solRad
        :param angle: polar coordinate in degrees
        :return:
        """
        # if radius was given transforms polar
        # to Cartesian coordinates
        # print r, angle
        if r > 0.:

            # convert degrees to radians
            if degrees:
                angle = np.radians(angle.value) % (2 * np.pi)

            self['xsp'] = r * np.cos(angle)
            self['ysp'] = r * np.sin(angle)
        else:
            self['xsp'] = x
            self['ysp'] = y

    def set_velocity(self, vx=0.0, vy=0.0, v=0.0, angle=0.0, degrees=True):
        """
        Sets the net velocity of the spot
        :param vx: Cartesian coordinate in solRad
        :param vy: Cartesian coordinate in solRad
        :param v: polar coordinate in solRad
        :param angle: polar coordinate in degrees
        :return:
        """
        if abs(v) > 0.:
            if degrees:
                angle = np.radians(angle.value) % (2 * np.pi)
            self['vxsp'] = v * np.sin(angle)
            self['vysp'] = -v * np.cos(angle)
        else:
            self['vxsp'] = vx
            self['vysp'] = vy


class Ufo(BaseObject):
    """
    Wraps the object ufo of Shellspec.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        """

        # initialize the parent class
        super(Ufo, self).__init__(name='ufo', parameters=[], **kwargs)


class Jet(BaseObject):
    """
    Wraps the object jet of Shellspec.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        """

        # initialize the parent class
        super(Jet, self).__init__(name='jet', parameters=[], **kwargs)

    def resolve_constraints(self):
        """
        If the polar coordiantes of the jet
        are fitted. This function propagates
        any changes to the xjt and yjt
        :return:
        """
        # propagates polar coordinates into Cartesian
        if self.get_parameter('rpoljt', 'value') > 0.0:
            self.set_position(r=self.get_parameter('rpoljt'), angle=self.get_parameter('pangjt'))

        if abs(self.get_parameter('vpoljt', 'value')) > 0.0:
            self.set_velocity(v=self.get_parameter('vpoljt'), angle=self.get_parameter('pangjt'))

    def set_position(self, x=0.0, y=0.0, r=0.0, angle=0.0, degrees=True):
        """
        Sets the position of the jet
        :param x: Cartesian coordinate in solRad
        :param y: Cartesian coordinate in solRad
        :param r: polar coordinate in solRad
        :param angle: polar coordinate in degrees
        :return:
        """
        if r > 0.:
            if degrees:
                angle = np.radians(angle.value) % (2 * np.pi)
            self['xjt'] = r * np.cos(angle)
            self['yjt'] = r * np.sin(angle)
        else:
            self['xjt'] = x
            self['yjt'] = y

    def set_velocity(self, vx=0.0, vy=0.0, v=0.0, angle=0.0, degrees=True):
        """
        Sets the net velocity of the jet
        :param vx: Cartesian coordinate in solRad
        :param vy: Cartesian coordinate in solRad
        :param v: polar coordinate in solRad
        :param angle: polar coordinate in degrees
        :return:
        """
        if abs(v) > 0.:
            if degrees:
                angle = np.radians(angle.value) % (2 * np.pi)
            self['vxjt'] = v * np.sin(angle)
            self['vyjt'] = -v * np.cos(angle)
        else:
            self['vxjt'] = vx
            self['vyjt'] = vy


class Flow(BaseObject):
    """
    Wraps the object flow of Shellspec.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        """

        # initialize the parent class
        super(Flow, self).__init__(name='flow', parameters=[], **kwargs)

    def resolve_constraints(self):
        """
        If the polar coordiantes of the flow
        are fitted. This function propagates
        any changes to the xfw and yfw
        :return:
        """
        # propagates polar coordinates into Cartesian
        # if they are fitted
        if self.get_parameter('rpolfw', 'value') > 0.0:
            self.set_position(r=self.get_parameter('rpolfw'), angle=self.get_parameter('pangfw'))

        if self.get_parameter('r12fw', 'value') > 0.0:
            self.set_radius(r=self.get_parameter('r12fw'))

        if self.get_parameter('z12fw', 'value') > 0.0:
            self.set_zcoord(z=self.get_parameter('z12fw'))

        if abs(self.get_parameter('v12fw', 'value')) > 0.0:
            self.set_velocity(v=self.get_parameter('v12fw'))

    def set_position(self, x=0.0, y=0.0, r=0.0, angle=0.0, degrees=True):
        """
        Sets the position of the flow
        :param x: Cartesian coordinate in solRad
        :param y: Cartesian coordinate in solRad
        :param r: polar coordinate in solRad
        :param angle: polar coordinate in degrees
        :return:
        """
        # if radius was given transforms polar
        # to Cartesian coordinates
        # print r, angle
        if r > 0.0:

            # convert degrees to radians
            if degrees:
                angle = np.radians(angle.value) % (2 * np.pi)

            self['x1fw'] = r * np.cos(angle)
            self['y1fw'] = r * np.sin(angle)
            self['x2fw'] = r * np.cos(angle)
            self['y2fw'] = r * np.sin(angle)
        else:
            self['x1fw'] = x
            self['y1fw'] = y
            self['x2fw'] = x
            self['y2fw'] = y

    def set_radius(self, r=0.0):
        """Sets radius at both ends of the flow"""
        if r > 0.0:
            self['r1fw'] = r
            self['r2fw'] = r

    def set_zcoord(self, z=0.0):
        """Sets z-coordinate of both ends of the flow"""
        if z > 0.0:
            self['z1fw'] = -z
            self['z2fw'] = z

    def set_velocity(self, v=0.0):
        """Sets velocity at both ends of the flow"""
        if abs(v) > 0.0:
            self['v1fw'] = -v
            self['v2fw'] = v


class Shell(BaseObject):
    """
    Class wrapping the shell of shellspec.
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        :return:
        """

        # initialize the parent class
        super(Shell, self).__init__(name='shell', parameters=[], **kwargs)


