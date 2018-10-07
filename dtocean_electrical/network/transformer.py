# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Adam Collin
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Created on Thu Apr 07 16:05:51 2016

.. moduleauthor:: Adam Collin <adam.collin@ieee.org>
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
