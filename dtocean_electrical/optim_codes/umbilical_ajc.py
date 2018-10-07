# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Sam Weller
#    Copyright (C) 2017-2018 Mathew Topper
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
Umbilical cable calculations

.. moduleauthor:: Sam Weller <S.Weller@exeter.ac.uk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import math
import logging
from collections import namedtuple

import numpy as np

# Start logging
module_logger = logging.getLogger(__name__)


class Variables(object):
    
    '''Collect input data for the umbilical design module.

    Args: 
        devices (list) [-]: List of device identification numbers.
        gravity (float) [m/s2]: Acceleration due to gravity.
        compdict (dict) [-]: Representation of db data for umbilical design
            module:
                key = 'item3', value = None,
                key = 'item5', value = [minimum break load,
                                        minimum bend radius]
                key = 'item6', value = [diameter]
                key = 'item7', value = [dry mass, wet mass]
        systype (str) [-]: Device type, from: 'tidefloat', 'tidefixed',
            'wavefloat', 'wavefixed'.
        sysorig () [-]: Seabed connection point of each device as (x, y, z)
            coordinates; key = device number, value = (x, y, z).
        umbconpt (np.ndarray) [m]: Umbilical connection point as (x, y, z)
            coordinates. x and y given with respect to device origin, z with
            respect to MSL.
        sysorienang (float) [deg]: System orientation angle.
        preumb (int) [-]: DB key of the selected umbilical.
        umbsf (float) [-]: Umbilical safety factor.
        subcabconpt (dict) [m]: Subsea cable connection point for each device 
            as (x, y, z) coordinates.
        sysdraft (float) [m]: Device equilibrium draft without mooring/cable
            system.

    Attributes:
        None

    '''

    def __init__(self, devices,
                       gravity,
                       compdict,
                       systype,
                       sysorig,
                       umbconpt,
                       sysorienang,
                       preumb,
                       umbsf,
                       subcabconpt,
                       sysdraft):

        self.devices = devices
        self.gravity = gravity
        self.compdict = compdict
        self.systype = systype
        self.sysorig = sysorig
        self.umbconpt = umbconpt
        self.sysorienang = sysorienang
        self.preumb = preumb
        self.umbsf = umbsf
        self.subcabconpt = subcabconpt
        self.sysdraft = sysdraft

        return


class Umbilical(object):
    
    '''Umbilical geometry specification submodule 
    
        Args:
            variables (object) [-]: Instance of Variables class.

        Attributes:
            selumbtyp: (str) [-]: Selected umbilical type.
            umbgeolw (numpy.ndarray): Not yet fully defined.
            umbleng (float) [m]: Umbilical length.
            totumbcost (float) [euros]: Umbilical total cost.
            umbbomat (dict) [-]: Umbilical bill of materials, keys:
                'umbilical type' [str], 'length' [m], 'cost' [euros]
                'total weight' [kg], 'diameter' [m].
            umbhier (list) [-]: Umbilical hierarchy.

        Functions:
            umbdes: Specifies umbilical geometry.
            umbinst: Calculates installation parameters.
            umbcost: Calculates umbilical capital cost.
            umbbom: Creates umbilical bill of materials.
            umbhier: Creates umbilical hierarchy.

    '''

    def __init__(self, variables):
        self._variables = variables

    def umbdes(self, deviceid, syspos):
        """ This method will be used to look-up umbilical properties """
        
        # Define geometry
        # Umbilical defined by WP3
        preumb_record = self._variables.compdict[self._variables.preumb]
        self.selumbtyp = preumb_record['item3']
        
        if self._variables.systype in ("wavefloat", "tidefloat"):
            # """ Lazy-wave geometry comprises three sections; hang-off,
            # buoyancy and decline """
            compblocks = ['hang off', 'buoyancy', 'decline']
        elif self._variables.systype in ("wavefixed", "tidefixed"):
            compblocks = ['hang off']
            
        prop_cols = ['compid',
                     'size',
                     'length',
                     'dry_mass',
                     'wet_mass',
                     'mbl',
                     'mbr']
        UmpProps = namedtuple("UmpProps", prop_cols)
        
        umbwetmass = [0 for row in range(len(compblocks))]
        umbconpt_rotated = [0 for row in range(0, 2)]
        umbtopconn = [0 for row in range(0, 3)]
        angle_rads = -self._variables.sysorienang * math.pi / 180.0
        pre_rotated = self._variables.umbconpt[:2]
        
        # Rotate the connection point
        umbconpt_rotated[0] = pre_rotated[0] * math.cos(angle_rads) - \
                                      pre_rotated[1] * math.sin(angle_rads)
        umbconpt_rotated[1] = pre_rotated[0] * math.sin(angle_rads) + \
                                      pre_rotated[1] * math.cos(angle_rads)
        
        # Move the connection point to global coordinates
        umbtopconn[0] = round(umbconpt_rotated[0] + syspos[0], 3)
        umbtopconn[1] = round(umbconpt_rotated[1] + syspos[1], 3)
        
        if self._variables.systype in ("wavefloat","tidefloat"):  
            umbtopconn[2] = self._variables.umbconpt[2] - \
                                                    self._variables.sysdraft
            klim = 1
        elif self._variables.systype in ("wavefixed","tidefixed"):
            umbtopconn[2] = self._variables.umbconpt[2]
            klim = 500
            
        # Umbilical length set initially as 1.15 x the shortest distance
        # between the upper and lower connection points
        subcabconpts = self._variables.subcabconpt[deviceid]
        
        umbleng = 1.15 * math.sqrt((umbtopconn[0] - subcabconpts[0]) ** 2.0 
                                 + (umbtopconn[1] - subcabconpts[1]) ** 2.0
                                 + (umbtopconn[2] - subcabconpts[2]) ** 2.0)
        umblengcheck = 'False'
        
        for k in range(0,klim):
            # logmsg = ('Umbilical total length {}').format(umbleng)
            # module_logger.debug(logmsg)
            
            # If any z-coordinate is below the global subsea cable connection
            # point reduce umbilical length by 0.5%
            if (self._variables.systype in ("wavefixed","tidefixed") and
                k > 0 and umblengcheck == 'False'):
                
                umbleng = 0.999 * umbleng
                
            if self._variables.systype in ("wavefloat","tidefloat"):                         
                # Initial lengths of hang-off, buoyancy and decline sections
                # (40%, 20% and 40%
                umbsecleng = [0.4 * umbleng, 0.2 * umbleng, 0.4 * umbleng]
            elif self._variables.systype in ("wavefixed","tidefixed"):
                self.umbleng = umbleng
                umbsecleng = [umbleng]
                
            umbcomptab = {}
                            
            for i, block_name in enumerate(compblocks):
                
                if compblocks[i] == 'buoyancy':
                    umbwetmass[i] = -1.4 * preumb_record['item7'][1]
                else:
                    umbwetmass[i] = preumb_record['item7'][1]
                    
                umbprops = UmpProps(self._variables.preumb, 
                                    preumb_record['item6'][0], 
                                    umbsecleng[i], 
                                    preumb_record['item7'][0], 
                                    umbwetmass[i], 
                                    preumb_record['item5'][0], 
                                    preumb_record['item5'][1])

                umbcomptab[block_name] = umbprops
                
            mlim = 5000
            # """ Catenary tolerance """
            tol = 0.01 
            # """ Distance tolerance """
            disttol = 0.001    
            # """ Number of segements along cable """
            numseg = 50
            flipzumb = [0 for row in range(0,numseg)]
            # """ Segment length """
            ds = umbleng / numseg 
            umbxf = math.sqrt((umbtopconn[0] - subcabconpts[0]) ** 2 
                            + (umbtopconn[1] - subcabconpts[1]) ** 2)
            umbzf = umbtopconn[2] - subcabconpts[2]
            
            Humb = [0 for row in range(numseg)]
            Vumb = [0 for row in range(numseg)] 
            Tumb = [0 for row in range(numseg)] 
            thetaumb = [0 for row in range(numseg)] 
            xumb = [0 for row in range(numseg)] 
            zumb = [0 for row in range(numseg)] 
            leng = [0 for row in range(numseg)]
            
            # Approximate catenary profile used in first instance to estimate
            #top end loads
            if umbxf == 0:            
                lambdacat = 1e6
            elif math.sqrt(umbxf ** 2 + umbzf ** 2) >= umbleng:
                lambdacat = 0.2
            else:     
                lambdacat = math.sqrt(3 * (((umbleng ** 2 
                                             - umbzf ** 2) / umbxf ** 2) - 1))

            c1 = 0.5 * umbcomptab['hang off'].wet_mass * \
                                                    self._variables.gravity
            Hf = max(math.fabs(c1 * umbxf / lambdacat), tol)  
            Vf = c1 * ((umbzf / math.tanh(lambdacat)) + umbleng)     
            Tf = math.sqrt(Hf ** 2.0 + Vf ** 2.0)
            theta_0 = math.atan(Vf / Hf)
            Tumb[0] = Tf
            Humb[0] = Hf
            Vumb[0] = Vf
            thetaumb[0] = theta_0    
            xumb[0] = 0.0
            zumb[0] = 0.0
            
            for m in range(0, mlim):   
                if m >= 1:
                    if math.fabs(errumbzf) > disttol * umbzf:                          
                        if (np.diff(zumb[k-2:k+1:2]) == 0.0
                            and np.diff(zumb[k-3:k:2]) == 0.0):
                            Tffactor = 0.0001
                        else:
                            Tffactor = 0.001                        
                        if errumbzf > 0.0:
                            Tumb[0] = Tumb[0] + Tffactor * Tf
                        elif errumbzf < 0.0:
                            Tumb[0] = Tumb[0] - Tffactor * Tf                        
                        Vumb[0] = math.sqrt(Tumb[0] ** 2.0 - Humb[0] ** 2.0)
                        thetaumb[0] = math.atan(Vumb[0] / Humb[0])
                
                    if math.fabs(errumbxf) > disttol * umbxf:
                        if (np.diff(xumb[k-2:k+1:2]) == 0.0
                            and np.diff(xumb[k-3:k:2]) == 0.0):
                            Hffactor = 0.005
                        else:
                            Hffactor = 0.01
                        if errumbxf > 0.0:
                            Humb[0] = Humb[0] + 0.001 * Hf
                        elif errumbxf < 0.0:
                            Humb[0] = Humb[0] - 0.001 * Hf
                        
                        Vumb[0] = Vumb[0]
                        thetaumb[0] = math.atan(Vumb[0] / Humb[0])
                        if Humb[0] < 0.0:
                            thetaumb[0] = math.pi / 2.0
                            Humb[0] = 0.0
                            
                c2 = umbcomptab['hang off'].wet_mass * \
                                                self._variables.gravity * ds
                c3 = umbcomptab['buoyancy'].wet_mass * \
                                                self._variables.gravity * ds
                c4 = umbcomptab['decline'].wet_mass * \
                                                self._variables.gravity * ds
                
                for k in range(1, numseg):
                    leng[k] = ds * k 
                    if leng[k] <= umbsecleng[0]:
                        # """ Hang-off section """
                        Vumb[k] = Vumb[k-1] - c2
                        Humb[k] = Humb[k-1] 
                        thetaumb[k] = math.atan(Vumb[k] / Humb[k])
                        xumb[k] = xumb[k-1] + ds * math.cos(thetaumb[k-1])
                        zumb[k] = zumb[k-1] + ds * math.sin(thetaumb[k-1])
                    elif (leng[k] > umbsecleng[0] and
                          leng[k] <= sum(umbsecleng[0:2])): 
                        # """ Buoyancy section """
                        Vumb[k] = Vumb[k-1] - c3
                        Humb[k] = Humb[k-1]
                        thetaumb[k] = math.atan(Vumb[k] / Humb[k])
                        xumb[k] = xumb[k-1] + ds * math.cos(thetaumb[k-1])
                        zumb[k] = zumb[k-1] + ds * math.sin(thetaumb[k-1])
                    elif (leng[k] > sum(umbsecleng[0:2]) and
                          leng[k] <= umbleng):
                        # """ Decline section """
                        Vumb[k] = Vumb[k-1] - c4
                        Humb[k] = Humb[k-1]                
                        thetaumb[k] = math.atan(Vumb[k] / Humb[k])                
                        if Vumb[k] < 0.0:
                            thetaumb[k] = 0.0 
                            Vumb[k] = 0.0
                        xumb[k] = xumb[k-1] + ds * math.cos(thetaumb[k-1])                
                        zumb[k] = zumb[k-1] + ds * math.sin(thetaumb[k-1]) 
                    Tumb[k] = math.sqrt(Humb[k] ** 2.0 + Vumb[k] ** 2.0)     
                    # """ Allow cable to embed by up to 0.5m """
                    if (self._variables.systype in ("wavefixed", "tidefixed") 
                        and zumb[k] > umbzf + 0.5):
                        umblengcheck = 'False'                        
                        break
                    else:
                        umblengcheck = 'True'
                    
                errumbxf = umbxf - xumb[-1]
                errumbzf = umbzf - zumb[-1]
                if umblengcheck == 'False':
                    break
                
                if (math.fabs(errumbxf) < disttol * umbxf and
                    math.fabs(errumbzf) < disttol * umbzf):
                    # logmsg = ('Solution found, umbilical length '
                    #           '{}').format(umbleng)
                    # module_logger.debug(logmsg)  
                    break  
            
            if umblengcheck == 'True':               
               break    
            
        # """ Maximum tension """
        self.Tumbmax = max(Tumb)
        # """ Radius of curvature along umbilical (starting at device end) """
        dzdx = np.diff(zumb) / np.diff(xumb) 
        d2zdx2 = np.diff(dzdx) / np.diff(xumb[:-1])
        umbradcurv = [abs(number) for number in
                                  (((1 + dzdx[0:-1] ** 2.0) ** 1.5) / d2zdx2)]
        self.umbradmin = min(umbradcurv)
            
        # """ Umbilical load angle from y-axis """
        deltaumbx = subcabconpts[0] - umbtopconn[0]
        deltaumby = subcabconpts[1] - umbtopconn[1]
        posloadang = math.atan(deltaumbx / deltaumby)
        
        if (deltaumbx > 0.0 and deltaumby > 0.0):
            Humbloadang = posloadang
        elif (deltaumbx > 0.0 and deltaumby < 0.0):
            Humbloadang = posloadang + 90.0 * math.pi / 180.0
        elif (deltaumbx < 0.0 and deltaumby < 0.0):
            Humbloadang = posloadang + 180.0 * math.pi / 180.0
        elif (deltaumbx < 0.0 and deltaumby > 0.0):
            Humbloadang = posloadang + 270.0 * math.pi / 180.0
        elif (deltaumbx == 0.0 and deltaumby > 0.0):
            Humbloadang = 0.0
        elif (deltaumbx > 0.0 and deltaumby == 0.0):
            Humbloadang = 90.0 * math.pi / 180.0
        elif (deltaumbx == 0.0 and deltaumby < 0.0):
            Humbloadang = 180.0 * math.pi / 180.0
        elif (deltaumbx < 0.0 and deltaumby == 0.0):
            Humbloadang = 270.0 * math.pi / 180.0                                  
        
        HumbloadX = Humb[0] * math.sin(Humbloadang)
        HumbloadY = Humb[0] * math.cos(Humbloadang)
        Vumbload = Vumb[0] 
            
        if (self.Tumbmax > self._variables.umbsf * preumb_record['item5'][0] or
            self.umbradmin < preumb_record['item5'][1]):
            umbtencheck = 'False'
        else:
            umbtencheck = 'True'
        
        for zind, zval in enumerate(zumb):
            flipzumb[zind] = umbtopconn[2] - zval
        
        # logmsg = ('X coord {}').format(xumb)
        # module_logger.debug(logmsg) 
        # logmsg = ('Z coord {}').format(flipzumb)
        # module_logger.debug(logmsg) 
        
        return umbleng, xumb, flipzumb
