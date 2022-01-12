#!/usr/bin/env python3

from .objects import Distance
from .objects import Environmental
from .objects import Grid

class Model:
    """
    In this container all shellspec objects 
    that will be modelled + the orbit are stored.
    """
    
    def __init__(self, name='', objects=[]):
        """
        """
        
        # available objects and their names
        self.__registered_objects = [obj.get_name() for obj in objects]
        
        # set the object name and pass models
        self.__name = name, 
        self.__objects = {obj.get_name():obj for obj in objects}

        # some shellspec objects are common for every
        # model -- if they are not provided, they will
        # be attached automatically
        if not self.has_object('distance'):
            self.add_object(Distance())
        if not self.has_object('grid'):
            self.add_object(Grid())
        if not self.has_object('enviromental'):
            self.add_object(Environmental())
        
    def add_object(self, obj):
        """
        """
        
        # register the name
        obj_name = obj.get_name()
        
        # register the object
        self.__registered_objects.append(obj_name)
        self.__objects[obj_name] = obj

        
    def has_object(self, objname):
        """
        Checks that a shellspec object was registered.
        """
        
        if objname.lower() in self.__registered_objects:
            return True
        else:
            return False

    def keys(self):
        """
        Returns names of all registered objects
        :return:
        """

        return self.__registered_objects
        
    def __getitem__(self, name):
        """
        Returns an object.
        """
        return self.__objects[name]

    def __str__(self):
        """
        Returns string representation of the class.
        :return:
        """

        string = ''
        for objname in list(self.keys()):
            string += str(self[objname])
        return string
