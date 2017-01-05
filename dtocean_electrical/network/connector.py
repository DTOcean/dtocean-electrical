# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 16:05:50 2016

@author: acollin
"""


class Connector(object):

    def __init__(self, index, db_key, marker, loc):

        # what attributes does a connector have?
        self.type_ = None

        # what attributes do we impose on a connector?
        self.id_ = index
        self.db_key = db_key
        self.marker = marker
        self.utm_x = loc[0]
        self.utm_y = loc[1]

    def __str__(self):

        '''Override print command to display some info.

        '''

        return ('This is the ' + self.__class__.__name__ + ' class.')


class WetMateConnector(Connector):

    def __init__(self, index, db_key, marker, loc):
        super(WetMateConnector, self).__init__(index, db_key, marker, loc)
        self.type_ = 'wet-mate'


class DryMateConnector(Connector):

    def __init__(self, index, db_key, marker, loc):
        super(DryMateConnector, self).__init__(index, db_key, marker, loc)
        self.type_ = 'dry-mate'
