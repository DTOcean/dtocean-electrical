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
Created on Thu Apr 07 13:41:07 2016

.. moduleauthor:: Adam Collin <adam.collin@ieee.org>
"""

# from transformer import Transformer
# from connector import Connector, WetMateConnector, DryMateConnector
# import power_quality
# import switchgear


class CollectionPoint(object):

    '''Class to define all properties of the offshore collection point.

    '''

    def __init__(self, index, loc, db_key, data):

        # set attributes of a collection point have
        self.electrical_type_ = None
        self.operating_environment = data.operating_environment.values[0]
        self.n_inputs = data.input.values[0]
        self.n_output = data.output.values[0]
        self.input_type = data.input_connector.values[0]
        self.output_type = data.output_connector.values[0]
        self.foundation_type = data.foundation.values[0]
        self.mass = data.dry_mass.values[0]
        self.centre_of_gravity = data.gravity_centre.values[0]
        self.wet_frontal_area = data.wet_frontal_area.values[0]
        self.wet_beam_area = data.wet_beam_area.values[0]
        self.dry_frontal_area = data.dry_frontal_area.values[0]
        self.dry_beam_area = data.dry_beam_area.values[0]
        self.length = data.depth.values[0]
        self.width = data.width.values[0]
        self.height = data.height.values[0]
        self.volume = self.length * self.width * self.height
        self.profile = 'rectangular'
        self.surface_roughness = 1e-6
        self.orientation_angle = data.orientation_angle.values[0]
        self.foundation_locations = data.foundation_loc.values[0]
#        self.subsea = self._check_operating_environment(
#            data.operating_environment.values[0])

        # se attributes we place on the collection point
        self.id_ = index
        self.db_key = db_key
        self.location = loc
        self.utm_x = loc[0]
        self.utm_y = loc[1]
        self.input_connectors = data.input_connector.item()
        self.output_connectors = data.output_connector.item()
        self.upstream = None  # these are the connectors and added later
        self.downstream = None  # these are the connectors and added later
        self.marker = None  # the network marker is added later

    def __str__(self):

        '''Override print command to display some info'''

        return ('This is the ' + self.__class__.__name__ + ' class. This ' +
                self.__class__.__name__ + ' has ' + str(self.n_inputs) +
                ' inputs.')


class PassiveHub(CollectionPoint):

    '''Special instance of CollectionPoint class.

    '''

    def __init__(self, index, loc, db_key, data):
        super(PassiveHub, self).__init__(index, loc, db_key, data)
        self.type_ = 'passive hub'
        self.configuration = 'busbar'
        self.subsea = True

        # create the connectors
        # input

        # output


class Substation(CollectionPoint):

    '''Special instance of CollectionPoint class.

    '''

    def __init__(self, index, loc, db_key, data):
        super(Substation, self).__init__(index, loc, db_key, data)
        self.type_ = 'substation'
        self.subsea = True
#        self.transformer = Transformer(10)
        self.configuration = ['busbar', 'transformer', 'disconnector']
