# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 16:05:51 2016

@author: acollin
"""


class Transformer(object):

    def __init__(self):

        # What attributes does the transformer have?
        self.voltage_1 = None
        self.voltage_2 = None
        self.impedance = None
        self.cost = None
        self.mass = None
        self.height = None
        self.weight = None
        self.length = None

        # What attributes do we place on the transformer?
        self.id_ = None
        self.db_key = None
        self.location = None  # id of collection point to which it is connected

    def __str__(self):

        '''Override print command to display some info'''

        return ('This is the ' + self.__class__.__name__ + ' class.')
