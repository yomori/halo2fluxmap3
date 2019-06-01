#instrument.py
from astropy.io import fits
import numpyas np
import healpy as hp

def get_transmission(exp,freq):
	''' Given experiment and channel, return frequency/transmission'''
    if exp=='planck':
        d = fits.open('data/HFI_RIMO_R3.00.fits')
        f = d['BANDPASS_F%d'%freq].data['WAVENUMBER']*3e8*1e-7
        T = d['BANDPASS_F%d'%freq].data['TRANSMISSION']
    return f,T
