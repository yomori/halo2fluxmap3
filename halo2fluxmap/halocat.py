import os
import datetime

import numpy              as np
import pylab              as pl
import scipy              as sp
import cosmolopy.distance as cosdist

from astropy.io import fits
from cosmolopy  import fidcosmo
from globals    import *

class HaloLightCone():
    """
    @brief Class describing a halo light cone
    """
    def __init__(self, **kwargs):
        pass

    def copy(self):
        """
        @brief Creates a copy of the halo light cone
        """
        return copy.copy(self)

def ReadPkscLightCone(i):
    """
    @brief Reads peakpatch output in binary format into an object.
    """    
    halos = HaloLightCone()

    filename = params.folder+params.inputfile
    Nreads   = params.Nreads

    pkfile   = open(filename,"rb")

    halos.Non        = np.fromfile(pkfile, dtype=np.int32, count=1).astype('int')[0]
    halos.RTHmax     = np.fromfile(pkfile, dtype=np.float32, count=1).astype('float')[0]
    halos.redshiftin = np.fromfile(pkfile, dtype=np.float32, count=1).astype('float')[0]

    start = int(np.ceil(float(halos.Non)/Nreads)*i)
    end   = int(min( halos.Non , np.ceil(float(halos.Non)/Nreads)*(i+1) ))

    Nread = end-start

    outnum  = 10 #7 floats per halo
    npkdata = Nread*outnum

    pkfile.seek( (3+start*outnum)*4 ,os.SEEK_SET)  #12 byte header plus array
    halos.data = np.fromfile(pkfile, dtype=np.float32, count=npkdata)
    halos.data = np.reshape(halos.data,(Nread,outnum)).astype('float32')
    pkfile.close()

    rho       = 2.775e11*params.omegam*params.h**2
    halos.data[:,6] = (4*(np.pi)/3)*(halos.data[:,6]**3)*rho # convert 6 from RTH to M

    halos.data = np.column_stack((
        halos.data[:,0],
        halos.data[:,1],
        halos.data[:,2],
        halos.data[:,6]))            

    return halos

def ReadLightConeHeader():

    filename = params.folder+params.inputfile

    pkfile = open(filename,"rb")
    Non    = np.fromfile(pkfile, dtype=np.int32, count=1).astype('int')[0]
    
    return Non
