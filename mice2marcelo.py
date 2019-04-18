from numpy.random import *
from astropy.io import fits
import numpy as np
from cosmotools import cosmo
import camb
from astropy.io import fits
from scipy import interpolate
from camb import model, initialpower
from cosmotools import utils
import sys


pars = camb.CAMBparams()
pars.set_accuracy(AccuracyBoost=3,lSampleBoost=1, lAccuracyBoost=1)
pars.set_cosmology(H0=68.81, ombh2=0.0221589393480, omch2=0.116902635602, mnu=0.06, nnu=3.046, omk=0, tau=0.08,YHe=0.245352)
pars.set_dark_energy(),
pars.InitPower.set_params(ns=0.9676,As=2.260574e-09,pivot_scalar=0.05) # some test As which we fix later
pars.set_matter_power(redshifts=[0.], kmax=500)
pars.NonLinear = model.NonLinear_both
results = camb.get_results(pars)


N      = 10000 # number of objects
chimin = 1     # minimum distance [Mpc]
chimax = 5e3   # maximum distance [Mpc]
latmin = -5.   # minimum latitute [deg]
latmax =  5.   # maximum latitute [deg]
lngmin = -5    # minimum longitute [deg]
lngmax =  5    # maximum longitute [deg]
Mmin   = 11.   # minimum log10 halo mass [Msun]
Mmax   = 15.   # maximum log10 halo mass [Msun]

fileout = open('./catalogs/halocatalog_mice-0.pksc','w')
totpx=[]
totpy=[]
totpz=[]
totRTH=[]
for patch in range(0,1):
	fname='/project2/chihway/sims/MICEv2/halos/halocat.fits'
	y       = fits.open(fname)
	print('loading %s'%fname)
	ids  =y[1].data['unique_halo_id']
	N=len(ids)
	print "reading %d halos"%N
	zs  = y[1].data['z_cgal']
	ra  = y[1].data['ra_gal']
	dec = y[1].data['dec_gal']
	Mh  = y[1].data['lmhalo']
	Rh  = y[1].data['RVIR']

	good_all = np.where(((ids < 1e8) & (zs < 0.34)) | ((ids >= 1e8) & (ids <= 1e9) & (zs > 0.34) & (zs < 0.9)) | ((ids > 1e9) & (zs > 0.9)))[0]
	#good_all = np.where( (ids < 1e8) & (zs < 0.34) )[0]

	print "restricting to %d halos"%len(good_all)

	zs  = y[1].data['z_cgal'][good_all]
	ra  = y[1].data['ra_gal'][good_all]
	dec = y[1].data['dec_gal'][good_all]
	Mh  = y[1].data['lmhalo'][good_all]
	Rh  = y[1].data['RVIR'][good_all]
	px  = y[1].data['zhalo'][good_all]
	py  = y[1].data['yhalo'][good_all]
	pz  = y[1].data['zhalo'][good_all]

	#zs = np.ones_like(zs)*0.5
	tht,phi = utils.rd2tp(ra,dec)
	chis = results.comoving_radial_distance(zs)

	# mean matter density
	omegam = 0.3
	h      = 0.7
	rho = 2.775e11 * omegam * h**2

	mu   = np.cos(tht)

	# generate random positions and halo masses
	chi = chis#uniform(chimin, chimax, N)
	mu  = mu#uniform(mumin,   mumax, N)
	phi = phi#uniform(phimin, phimax, N) 
	M   = Mh#5e13*np.ones_like(zs)#uniform(Mmin,     Mmax, len(zs))
	z   = chi * mu
	#r   = chi * np.sqrt(1.-mu**2)
	#x   = r   * np.cos(phi)
	#y   = r   * np.sin(phi)
	x    = px
	y    = py
	z    = pz
	totpx=np.append(totpx,x)
	totpy=np.append(totpy,y)
	totpz=np.append(totpz,z)

	# convert from log10 M to M
	#M = 10**M

	# convert from M to Lagrangian radius RTH
	RTH = (3*M/4./np.pi/rho)**(1./3.)
	totRTH=np.append(totRTH,RTH)
	# write data to pksc file oriented such that longitude = latitude = 0 is along z-axis

# zeros for velocities and lagrangian positions
ph = np.zeros(len(totRTH))

np.asarray([len(totRTH)]).astype(np.int32).tofile(fileout)
np.asarray([totRTH.max()]).astype(np.float32).tofile(fileout)
np.asarray([0]).astype(np.float32).tofile(fileout)
np.asarray([totpy,totpz,totpx,ph,ph,ph,totRTH,ph,ph,ph]).transpose().astype(np.float32).tofile(fileout)

fileout.close()


