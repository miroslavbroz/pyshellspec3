#!/usr/bin/env python3

import sys
from astropy import units

class Parameter:
    """
    Each physical parameter is wrapped within this class.
    """
    
    def __init__(self, name='dummy', value=0.5, unit=None, vmin=0.0, vmax=1.0, fitted=False, dtype=float, docstring=''):
        """
        Constructs the class
        """
        
        # pass the arguments
        self.name = name
        self.unit = unit

        # some parameters can be mere switches, they
        # should remain type integer
        if dtype == int:
            self.value = value
            self.vmin = vmin
            self.vmax = vmax

        # floats are transcripted into units.Quantities
        else:
            self.value = units.Quantity(value=value, unit=unit)
            self.vmin = units.Quantity(value=vmin, unit=unit)
            self.vmax = units.Quantity(value=vmax, unit=unit)

        self.fitted = fitted
        self.dtype = dtype
        self.docstring = docstring
        
        # check boundaries and type
        self.check_boundaries()
        self.check_type()
        
        # define adjustable parameters
        self.__adjustable = ['value', 'vmin', 'vmax', 'fitted']

    def check_boundaries(self):
        """
        Checks that value lies within boundaries.
        """
        if (self.vmin > self.value) or (self.vmax < self.value):
            raise ValueError('Value %s of paramters %s does not lie within boundaries (%s, %s). ' %
                             (str(self.value), self.name, str(self.vmin), str(self.vmax)))
    def check_type(self):
        """
        Checks that value, and boundaries have correct 
        type.
        """
        for attr in ['value', 'vmin', 'vmax']:
            if not isinstance(getattr(self, attr), units.Quantity):
                if not isinstance(getattr(self, attr), self.dtype):
                    raise TypeError('Attribute %s has incorrect type!' \
                                ' Should be %s, but found %s.' % (attr, str(self.dtype), str(type(getattr(self, attr)))))
            else:
                if not isinstance(getattr(self, attr).value, self.dtype):
                    raise TypeError('Attribute %s has incorrect type!' \
                                ' Should be %s, but found %s.' % (attr, str(self.dtype), str(type(getattr(self, attr).value))))
        
    def get_property(self, attr='value'):
        """
        Returns value of an attribute.
        """
        
        return getattr(self, attr.lower())
    
    def set_property(self, attr, val):
        """

        Sets value of an property.
        :param attr: attribute that will be changed
        :param val: new value of teh attribute
        """
        
        # convert to lower case
        attr = attr.lower()
        
        # only some properties can be changed
        if attr not in self.__adjustable:
            raise ValueError('The only adjustable attributes of the class are: %s' % (str(self.__adjustable)))
        
        elif attr in ['value', 'vmin', 'vmax']:
            # first set attributes type
            self.check_type()

            # attach the attribute
            if self.dtype == int:
                setattr(self, attr, val)
            else:
                setattr(self, attr, units.Quantity(val, unit=self.unit))
            
            # check that we are still within boundaries if value was changed
            if attr == 'value':
                self.check_boundaries()
        else:
            setattr(self, attr, val)

    def __getattr__(self):
        """
        Trial...
        """
        return None

        
    def __str__(self):
        """
        String representation of the class.
        """
        
        string = ''
        for attr in ['name', 'value', 'unit', 'vmin', 'vmax', 'fitted', 'dtype', 'docstring']:
            string += '%s:%s, ' % (attr, str(getattr(self, attr)))
        
        return string
    
