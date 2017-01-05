"""

Input variables are listed in WP4 Input list.xlsx

Module author: Sam Weller

"""
# Start logging
import logging
module_logger = logging.getLogger(__name__)

import pandas as pd
import math
import numpy as np


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
        self.selumbtyp = self._variables.preumb        
        """ Define geometry """        
        """ Umbilical defined by WP3 """
        self.selumbtyp = self._variables.compdict[self._variables.preumb]['item3']
        if self._variables.systype in ("wavefloat","tidefloat"):
            """ Lazy-wave geometry comprises three sections; hang-off, buoyancy and decline """
            compblocks = ['hang off', 'buoyancy', 'decline']
        elif self._variables.systype in ("wavefixed","tidefixed"):
            compblocks = ['hang off']
        umbprops = [0 for row in range(len(compblocks))]
        umbwetmass = [0 for row in range(len(compblocks))]
        umbtopconn = [0 for row in range(0, 3)]
        """ Global position of umbilical top connection """
        umbtopconn[0] = round((self._variables.umbconpt[0] + syspos[0]) * math.cos(
                                    - self._variables.sysorienang * math.pi / 180.0)
                                    - (self._variables.umbconpt[1] + syspos[1])
                                    * math.sin(- self._variables.sysorienang 
                                    * math.pi / 180.0),3)
        umbtopconn[1] = round((self._variables.umbconpt[0] + syspos[0]) * math.sin(
                                    - self._variables.sysorienang * math.pi / 180.0)
                                    + (self._variables.umbconpt[1] + syspos[1])
                                    * math.cos(- self._variables.sysorienang 
                                    * math.pi / 180.0),3)
        if self._variables.systype in ("wavefloat","tidefloat"):  
            umbtopconn[2] = self._variables.umbconpt[2] - (syspos[2] - self._variables.sysdraft)  
            klim = 1
        elif self._variables.systype in ("wavefixed","tidefixed"):
            umbtopconn[2] = self._variables.umbconpt[2]
            klim = 500
            
        """ Umbilical length set initially as 1.15 x the shortest distance between the upper and 
            lower connection points """
        umbleng = 1.15 * math.sqrt((umbtopconn[0] - self._variables.subcabconpt[deviceid][0]) ** 2.0 
                                 + (umbtopconn[1] - self._variables.subcabconpt[deviceid][1]) ** 2.0
                                 + (umbtopconn[2] - self._variables.subcabconpt[deviceid][2]) ** 2.0)
        umblengcheck = 'False'
        
        for k in range(0,klim):
            # logmsg = ('Umbilical total length {}').format(umbleng)
            # module_logger.debug(logmsg)
            """ If any z-coordinate is below the global subsea cable connection point reduce 
                umbilical length by 0.5% """
            if (self._variables.systype in ("wavefixed","tidefixed") and k > 0 and umblengcheck == 'False'):
                umbleng = 0.999 * umbleng
                
            if self._variables.systype in ("wavefloat","tidefloat"):                         
                """ Initial lengths of hang-off, buoyancy and decline sections (40%, 20% and 40%) """
                umbsecleng = [0.4 * umbleng, 0.2 * umbleng, 0.4 * umbleng]
            elif self._variables.systype in ("wavefixed","tidefixed"):
                self.umbleng = umbleng
                umbsecleng = [umbleng]
            
            for i in range(0,len(compblocks)):
                if compblocks[i] == 'buoyancy':
                    umbwetmass[i] = -1.4 * self._variables.compdict[self._variables.preumb]['item7'][1]
                else:
                    umbwetmass[i] = self._variables.compdict[self._variables.preumb]['item7'][1]
                umbprops[i] = [self._variables.preumb, 
                            self._variables.compdict[self._variables.preumb]['item6'][0], 
                            umbsecleng[i], 
                            self._variables.compdict[self._variables.preumb]['item7'][0], 
                            umbwetmass[i], 
                            self._variables.compdict[self._variables.preumb]['item5'][0], 
                            self._variables.compdict[self._variables.preumb]['item5'][1]]                    
            """ Set up umbilical table """      
            colheads = ['compid', 'size', 'length', 'dry mass', 
                        'wet mass', 'mbl', 'mbr']     
            self.umbcomptab = pd.DataFrame(umbprops, 
                                            index=compblocks, 
                                            columns=colheads)                                            
            mlim = 5000
            """ Catenary tolerance """
            tol = 0.01 
            """ Distance tolerance """
            disttol = 0.001    
            """ Number of segements along cable """
            numseg = 50
            flipzumb = [0 for row in range(0,numseg)]
            """ Segment length """
            ds = umbleng / numseg 
            umbxf = math.sqrt((umbtopconn[0] - self._variables.subcabconpt[deviceid][0]) ** 2 
                            + (umbtopconn[1] - self._variables.subcabconpt[deviceid][1]) ** 2)
            umbzf = umbtopconn[2] - self._variables.subcabconpt[deviceid][2]
            
            Humb = [0 for row in range(numseg)]
            Vumb = [0 for row in range(numseg)] 
            Tumb = [0 for row in range(numseg)] 
            thetaumb = [0 for row in range(numseg)] 
            xumb = [0 for row in range(numseg)] 
            zumb = [0 for row in range(numseg)] 
            leng = [0 for row in range(numseg)]
            """ Approximate catenary profile used in first instance to estimate top end loads """
            if umbxf == 0:            
                lambdacat = 1e6
            elif math.sqrt(umbxf ** 2 + umbzf ** 2) >= umbleng:
                lambdacat = 0.2
            else:     
                lambdacat = math.sqrt(3 * (((umbleng ** 2 
                            - umbzf ** 2) / umbxf ** 2) - 1))                             
            Hf = max(math.fabs(0.5 * self.umbcomptab.ix['hang off','wet mass'] * self._variables.gravity * umbxf
                / lambdacat),tol)  
            Vf = 0.5 * self.umbcomptab.ix['hang off','wet mass'] * self._variables.gravity * ((umbzf
                / math.tanh(lambdacat)) + umbleng)     
            Tf = math.sqrt(Hf ** 2.0 + Vf ** 2.0)
            theta_0 = math.atan(Vf / Hf)
            Tumb[0] = Tf
            Humb[0] = Hf
            Vumb[0] = Vf
            thetaumb[0] = theta_0    
            xumb[0] = 0.0
            zumb[0] = 0.0
            for m in range(0,mlim):   
                if m >= 1:
                    if math.fabs(errumbzf) > disttol * umbzf:                          
                        if (np.diff(zumb[k-2:k+1:2]) == 0.0 and np.diff(zumb[k-3:k:2]) == 0.0):
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
                        if (np.diff(xumb[k-2:k+1:2]) == 0.0 and np.diff(xumb[k-3:k:2]) == 0.0):
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
                for k in range(1,numseg):
                    leng[k] = ds * k 
                    if leng[k] <= umbsecleng[0]:
                        """ Hang-off section """
                        Vumb[k] = Vumb[k-1] - self.umbcomptab.ix['hang off','wet mass'] * self._variables.gravity  * ds
                        Humb[k] = Humb[k-1] 
                        thetaumb[k] = math.atan(Vumb[k] / Humb[k])
                        xumb[k] = xumb[k-1] + ds * math.cos(thetaumb[k-1])
                        zumb[k] = zumb[k-1] + ds * math.sin(thetaumb[k-1])
                    elif (leng[k] > umbsecleng[0] and leng[k] <= sum(umbsecleng[0:2])): 
                        """ Buoyancy section """
                        Vumb[k] = Vumb[k-1] - self.umbcomptab.ix['buoyancy','wet mass'] * self._variables.gravity * ds
                        Humb[k] = Humb[k-1]
                        thetaumb[k] = math.atan(Vumb[k] / Humb[k])
                        xumb[k] = xumb[k-1] + ds * math.cos(thetaumb[k-1])
                        zumb[k] = zumb[k-1] + ds * math.sin(thetaumb[k-1])
                    elif (leng[k] > sum(umbsecleng[0:2]) and leng[k] <= umbleng):
                        """ Decline section """
                        Vumb[k] = Vumb[k-1] - self.umbcomptab.ix['decline','wet mass'] * self._variables.gravity * ds
                        Humb[k] = Humb[k-1]                
                        thetaumb[k] = math.atan(Vumb[k] / Humb[k])                
                        if Vumb[k] < 0.0:
                            thetaumb[k] = 0.0 
                            Vumb[k] = 0.0
                        xumb[k] = xumb[k-1] + ds * math.cos(thetaumb[k-1])                
                        zumb[k] = zumb[k-1] + ds * math.sin(thetaumb[k-1]) 
                    Tumb[k] = math.sqrt(Humb[k] ** 2.0 + Vumb[k] ** 2.0)     
                    """ Allow cable to embed by up to 0.5m """
                    if (self._variables.systype in ("wavefixed","tidefixed") 
                        and zumb[k] > umbzf + 0.5):
                        umblengcheck = 'False'                        
                        break
                    else:
                        umblengcheck = 'True'
                    
                errumbxf = umbxf - xumb[-1]
                errumbzf = umbzf - zumb[-1]
                if umblengcheck == 'False':
                    break
                
                if (math.fabs(errumbxf) < disttol * umbxf  and math.fabs(errumbzf) < disttol * umbzf):
                    # logmsg = ('Solution found, umbilical length {}').format(umbleng)
                    # module_logger.debug(logmsg)  
                    break  
            
            if umblengcheck == 'True':               
               break    
            
        """ Maximum tension """
        self.Tumbmax = max(Tumb)
        """ Radius of curvature along umbilical (starting at device end) """
        dzdx = np.diff(zumb) / np.diff(xumb) 
        d2zdx2 = np.diff(dzdx) / np.diff(xumb[:-1])
        umbradcurv = [abs(number) for number in (((1 + dzdx[0:-1] ** 2.0) ** 1.5) / d2zdx2)]
        self.umbradmin = min(umbradcurv)
            
        """ Umbilical load angle from y-axis """
        deltaumbx = self._variables.subcabconpt[deviceid][0] - umbtopconn[0]
        deltaumby = self._variables.subcabconpt[deviceid][1] - umbtopconn[1]
        if (deltaumbx > 0.0 and deltaumby > 0.0):
            Humbloadang = math.atan(deltaumbx / deltaumby) 
        elif (deltaumbx > 0.0 and deltaumby < 0.0):
            Humbloadang = math.atan(deltaumbx / deltaumby) + 90.0 * math.pi / 180.0
        elif (deltaumbx < 0.0 and deltaumby < 0.0):
            Humbloadang = math.atan(deltaumbx / deltaumby) + 180.0 * math.pi / 180.0
        elif (deltaumbx < 0.0 and deltaumby > 0.0):
            Humbloadang = math.atan(deltaumbx / deltaumby) + 270.0 * math.pi / 180.0
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
            
        if (self.Tumbmax > self._variables.umbsf * self._variables.compdict[self._variables.preumb]['item5'][0] 
            or self.umbradmin < self._variables.compdict[self._variables.preumb]['item5'][1]):
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
