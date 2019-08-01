#import pandas as pd
import pickle
#import cibmapping
import numpy as np
import camb
from camb import model, initialpower
#from tabulate import tabulate
#import universemachine as um
import healpy as hp
import sys
from astropy.modeling.blackbody import blackbody_nu as bb
from astropy.io import fits
import numexpr as ne

def tp2rd(tht,phi):
        """
        Convert tht,phi -> ra,dec
        """
        ra  = phi/np.pi*180.0
        dec = -1*(tht/np.pi*180.0-90.0)
        return ra,dec

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


def getnearestsnap(alist,zmid):
    zsnap  = 1/alist-1.
    return alist[np.argmin(np.abs(zsnap-zmid))]

#nsideout   = int(sys.argv[1]) # nside of the output CIB map
shellnum   = int(sys.argv[1]) # Shell index
shellwidth = int(sys.argv[2]) # Width of shell in Mpc/h

alist   = np.loadtxt('/scratch/users/yomori/mdpl2/halos/alist')
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

#ret      = np.zeros(hp.nside2npix(nsideout))

#--------------Loading the binary data file------------------------
d  = np.load('/scratch/users/yomori/mdpl2/halos/stripped_%.6f.npy'%nearestsnap)
m  = d[:,6]
idx=np.where(m>1e12)[0]
px = d[idx,0]
py = d[idx,1]
pz = d[idx,2]
m  = d[idx,6]
print("using %d halos"%len(idx))
del d
#------------------------------------------------------------------
#for xx in range(-ntiles,ntiles):
#for xx in range(0,ntiles):
#    for yy in range(-ntiles,ntiles):
#        for zz in range(-ntiles,ntiles):
totra  = np.array([])
totdec = np.array([])
totz   = np.array([])
totm   = np.array([])
for xx in range(-ntiles,ntiles):
    for yy in range(-ntiles,ntiles):
        for zz in range(-ntiles,ntiles):

            print("%d %d %d"%(xx,yy,zz))

            slicehit = checkslicehit(chilow,chiupp,xx,yy,zz)             # Check if box intersects with shell

            if slicehit==True:
                print('slicehit')

                for i in range(0,1):
                    sx  = ne.evaluate("px -%d + boxL * xx"%origin[0]) # dramatically faster to evaluate it this way
                    sy  = ne.evaluate("py -%d + boxL * yy"%origin[1])
                    sz  = ne.evaluate("pz -%d + boxL * zz"%origin[2])
                    r   = ne.evaluate("sqrt(sx*sx + sy*sy + sz*sz)")
                    #sx  = px - origin[0] + boxL * xx                      # positions in the tesselated space [Mpc/h]
                    #sy  = py - origin[1] + boxL * yy
                    #sz  = pz - origin[2] + boxL * zz
                    #del px,py,pz
                    #r   = np.sqrt(sx*sx + sy*sy + sz*sz)                    # comoving radial distance [Mpc/h]
                    zi  = results.redshift_at_comoving_radial_distance(r/h) # interpolated distance from position
                    idx = np.where((r>chilow) & (r<chiupp))[0]              # only select halos that are within the shell

                    if idx.size!=0:
                        tht,phi = hp.vec2ang(np.c_[sx[idx]/r[idx],sy[idx]/r[idx],sz[idx]/r[idx]])
                        ra,dec  = tp2rd(tht,phi)
	                totra  = np.append(totra,ra)
	                totdec = np.append(totdec,dec)	
                        totz   = np.append(totz,zi[idx])
                        totm   = np.append(totm,m[idx])

np.save('/scratch/users/yomori/mdpl2/halos/haloslc_%d.npy'%shellnum,np.c_[totra,totdec,totz,totm])
#nu,B = greybody(zmid) # frequency and the greybody spectrum
#d    = fits.open('/project2/chihway/yuuki/repo/halo2fluxmap3/data/HFI_RIMO_R3.00.fits')
#fint = np.arange(40001)

#for freq in (143,217,353,545,857):
#    f = d['BANDPASS_F%d'%freq].data['WAVENUMBER']*3e8*1e-7
#    T = d['BANDPASS_F%d'%freq].data['TRANSMISSION']
#    y = np.interp(fint,f*3,T)
#    f = np.sum(B*y)/np.sum(B) # fraction of total bolometric luminosity we receive
#    hp.write_map('/project2/chihway/sims/MDPL2/universemachine/cibmaps/cibmap_%d_%d_%d.fits'%(freq,chilow,chiupp),ret*f,overwrite=True)



