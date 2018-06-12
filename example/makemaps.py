#!/usr/bin/env python
import numpy             as np
import healpy            as hp
import matplotlib.pyplot as plt
import halo2fluxmap      as h2fm
import mpi4py.rc
import datetime
import os
import psutil
from scipy.ndimage.filters import gaussian_filter

h2fm.initialize('params.txt')

h2fm.report('Reading, Shuffling, Culling Catalogue',2)

# Get catalog of centrals
cen = h2fm.read_catalog()

# Get number of satellites for each halo
ns  = h2fm.cen2ns(cen)

# Populate halos with satellites
sat = h2fm.populate_centrals(cen,ns)

# Write time
h2fm.write_time("HOD completed", h2fm.params.rank)

#loop over frequencies
for inu_map in range(len(h2fm.params.freq_list)):

    # Put halos in map
    h2fm.params.nu_obs     = h2fm.params.freq_list[inu_map] * 1.e9
    h2fm.params.nu_obs_GHz = h2fm.params.freq_list[inu_map]
        
    # Give flux to gals and put in map
    intensity, flux_cen, flux_sat = h2fm.halos2map(cen,ns,sat) 
        
    # Write map
    ns_str   = str(h2fm.params.nside)
    nu_str   = str(int(h2fm.params.freq_list[inu_map]))    
    base = 'maps/'+h2fm.params.outbase 
    base += ('_ns'+ns_str+'_nu'+nu_str)

    h2fm.writemap(base,intensity)

# Write time
h2fm.write_time("Halo projection completed",h2fm.params.rank)

