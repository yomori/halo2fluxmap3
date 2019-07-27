import pandas as pd
import pickle
#import cibmapping
import numpy as np
import camb
from camb import model, initialpower
from tabulate import tabulate
#import universemachine as um
import healpy as hp
import sys
from astropy.modeling.blackbody import blackbody_nu as bb
from astropy.io import fits

def greybody(zmid):
    zmid = 2.0
    nu   = np.arange(40001)
    nu   = nu*1e9
    zc      = 2.; sigc  = 2. ; alpha = 2.; beta  = 2.; nup = 4955e9
    kb      = 1.38064852e-23; Tdust = 34.0; h = 6.62607004e-34
    nu_rest = nu*(1+zmid)

    fnu  = bb(nu, Tdust)*nu**beta
    fnup = bb(nup,Tdust)*nu**beta*(nu/nup)**-alpha

    F             = np.zeros_like(nu)
    F[nu<nup]     = fnu[nu<nup]
    F[nu>=nup]    = fnup[nu>=nup]
    F[nu<300e9]   = 0
    F[nu>37500e9] = 0
    return nu,F

def checkslicehit(chilow,chihigh,xx,yy,zz):
    # doing pre-selection so that we're not loading non-intersecting blocks 
    bvx=np.array([0,boxL,boxL,   0,   0,boxL,boxL,   0])
    bvy=np.array([0,   0,boxL,boxL,   0,   0,boxL,boxL])
    bvz=np.array([0,   0,   0,   0,boxL,boxL,boxL,boxL])

    sx  = (bvx - origin[0] + boxL * xx);
    sy  = (bvy - origin[1] + boxL * yy);
    sz  = (bvz - origin[2] + boxL * zz);
    r   = np.sqrt(sx*sx + sy*sy + sz*sz)

    if ( (np.all(chilow*0.8>r)) | (np.all(chihigh*1.2<r)) ):
        return False
    else:
        return True

def sfr2irlum(Mstar,SFR):
    #using modified model from 1611.04517
    KIR    = 1.49e-10
    KUV    = 1.71e-10
    IRX0   = 1.32
    alpha  = 1.5 
    IRX    = alpha*np.log10(Mstar/10**10.35)+IRX0
    irlum  = SFR/(KIR+KUV*10**(-IRX))
    return irlum

def getnearestsnap(alist,zmid):
    zsnap  = 1/alist-1.
    return alist[np.argmin(np.abs(zsnap-zmid))]

nsideout   = int(sys.argv[1]) # nside of the output CIB map
shellnum   = int(sys.argv[2]) # Shell index
shellwidth = int(sys.argv[3]) # Width of shell in Mpc/h

alist   = np.loadtxt('/project2/chihway/sims/MDPL2/universemachine/outputs4/alist')
zlist   = 1/alist-1.
origin  = [0,0,0]
boxL    = 1000
#-------- Running camb to get comoving distances -----------
h          = 0.6777
pars      = camb.CAMBparams()
pars.InitPower.set_params(ns=0.9611,As=np.exp(3.0663)*1e-10,pivot_scalar=0.05)# have to fiddle with this to match sigma8
pars.set_cosmology(H0=67.77, ombh2=0.022161 ,omch2=0.11889, tau=0.0952,num_massive_neutrinos=0,mnu=0,nnu=0)
pars.set_for_lmax(2000, lens_potential_accuracy=3)
pars.set_matter_power(redshifts=[0.], kmax=200.0)
pars.NonLinearModel.set_params(halofit_version='takahashi')
camb.set_feedback_level(level=100)
results   = camb.get_results(pars)
print(results.get_sigma8())

chilow = shellwidth*(shellnum+0)
chiupp = shellwidth*(shellnum+1)
chimid = 0.5*(chilow+chiupp)
ntiles = int(np.ceil(chiupp/boxL))
print("tiling [%dx%dx%d]"%(2*ntiles,2*ntiles,2*ntiles))
zmid   = results.redshift_at_comoving_radial_distance(chimid/h)
print('Generating map for halos in the range [%3.f - %.3f Mpc/h]'%(chilow,chiupp))

nearestsnap = getnearestsnap(alist,zmid)
print('The scalefactor closest to the middle of the shell is [%.6f]'%(nearestsnap))

ret      = np.zeros(hp.nside2npix(nsideout))

#--------------Loading the binary data file------------------------
dtype = np.dtype(dtype=[('id', 'i8'),('descid','i8'),('upid','i8'),
                        ('flags', 'i4'), ('uparent_dist', 'f4'),
                        ('pos', 'f4', (6)), ('vmp', 'f4'), ('lvmp', 'f4'),
                        ('mp', 'f4'), ('m', 'f4'), ('v', 'f4'), ('r', 'f4'),
                        ('rank1', 'f4'), ('rank2', 'f4'), ('ra', 'f4'),
                        ('rarank', 'f4'), ('A_UV', 'f4'), ('sm', 'f4'), 
                        ('icl', 'f4'), ('sfr', 'f4'), ('obs_sm', 'f4'), 
                        ('obs_sfr', 'f4'), ('obs_uv', 'f4'), ('empty', 'f4')],
                 align=True)

halos    = np.fromfile('/project2/chihway/sims/MDPL2/universemachine/outputs4/sfr_catalog_%.6f.bin'%nearestsnap, dtype=dtype)
obssm    = halos['obs_sm']
idx      = np.where(obssm>0)[0]   #only select halos with >0 stellar mass
obssm    = halos['obs_sm'][idx]
obssfr   = halos['obs_sfr'][idx]
IRlum    = sfr2irlum(obssm,obssfr)
del obssm,obssfr
px       = halos['pos'][:,0][idx]
py       = halos['pos'][:,1][idx]
pz       = halos['pos'][:,2][idx]

#------------------------------------------------------------------
d    = fits.open('/project2/chihway/yuuki/repo/halo2fluxmap3/data/HFI_RIMO_R3.00.fits')
fint  = np.arange(40001)
trans = np.zeros((40001,5))
c=0
for freq in (143,217,353,545,857):
    f = d['BANDPASS_F%d'%freq].data['WAVENUMBER']*3e8*1e-7
    T = d['BANDPASS_F%d'%freq].data['TRANSMISSION']
    y = np.interp(fint,f,T)
    trans[:,c]=y
    c+=1

#------------------------------------------------------------------
#for xx in range(-ntiles,ntiles):
for xx in range(0,ntiles):
    for yy in range(-ntiles,ntiles):
        for zz in range(-ntiles,ntiles):

            print("%d %d %d"%(xx,yy,zz))

            slicehit = checkslicehit(chilow,chiupp,xx,yy,zz)             # Check if box intersects with shell

            if slicehit==True:
                #totslicehit+=1
                print('slicehit')
                nu,B = greybody(zmid) # frequency and the greybody spectrum

                for i in range(0,1):
                    sx  = px - origin[0] + boxL * xx                      # positions in the tesselated space [Mpc/h]
                    sy  = py - origin[1] + boxL * yy
                    sz  = pz - origin[2] + boxL * zz
                    del px,py,pz
                    r   = np.sqrt(sx*sx + sy*sy + sz*sz)                    # comoving radial distance [Mpc/h]
                    zi  = results.redshift_at_comoving_radial_distance(r/h) # interpolated distance from position
                    idx = np.where((r>chilow) & (r<chiupp))[0]              # only select halos that are within the shell

                    if idx.size!=0:
                        print(len(idx))
                        pix     = hp.vec2pix(nsideout,sx[idx]/r[idx],sy[idx]/r[idx],sz[idx]/r[idx])
                        q       = np.arange(np.amax(pix)+1)
                        for freq in range(0,5): #143,217,353,545,857
                            Fnu     = np.sum(np.outer(B/np.sum(B)*trans[:,freq],IRlum[idx]/4/np.pi/(r[idx]**2)/(1+zi[idx])),axis=0)
                            ret[q]  = ret[q] + np.bincount(pix, weights=Fnu) # Convert luminosity to flux

hp.write_map('/project2/chihway/sims/MDPL2/universemachine/cibmaps/cibmap_bolo_%d_%d.fits'%(chilow,chiupp),ret,overwrite=True)

#nu,B = greybody(zmid) # frequency and the greybody spectrum
#d    = fits.open('/project2/chihway/yuuki/repo/halo2fluxmap3/data/HFI_RIMO_R3.00.fits')
#fint = np.arange(40001)

#for freq in (143,217,353,545,857):
#    f = d['BANDPASS_F%d'%freq].data['WAVENUMBER']*3e8*1e-7
#    T = d['BANDPASS_F%d'%freq].data['TRANSMISSION']
#    y = np.interp(fint,f*3,T)
#    f = np.sum(B*y)/np.sum(B) # fraction of total bolometric luminosity we receive
#    hp.write_map('/project2/chihway/sims/MDPL2/universemachine/cibmaps/cibmap_%d_%d_%d.fits'%(freq,chilow,chiupp),ret*f,overwrite=True)



