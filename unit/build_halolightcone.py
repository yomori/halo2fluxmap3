import pandas as pd 
import pickle
#import cibmapping
import numpy as np
import camb 
from camb import model, initialpower
from tabulate import tabulate
#import universemachine as um

#shellnum   = int(sys.argv[1])
chilow = 3550#float(sys.argv[1])#shellwidth*(shellnum+0)
chiupp = 3700#float(sys.argv[2])#shellwidth*(shellnum+1)

nsideout   = 2048
origin     = [0,0,0]
shellwidth = 50       # in Mpc/h 
boxL       = 1000     # in Mpc/h
maxchi     = 9000     # in Mpc/h --> 9000~z=6
minchi     = 5000     # in Mpc/h --> 5000~z=2

alist      = '/project2/chihway/sims/UNIT/UNITSIMS_GADGET/fixedAmp_002/hlist/alist'
dir_hlist  = '/project2/chihway/sims/UNIT/UNITSIMS_GADGET/fixedAmp_002/hlist/' #Path where the hlist files are
#dir_hlist  = '/project2/chihway/sims/UNIT/UNITSIMS_GADGET/fixedAmp_001/ROCKSTAR_HALOS/unitsims.ft.uam.es/DATABASE/UNITSIMS_GADGET/ROCKSTAR_HALOS/fixedAmp_001/hlist/' #Path where the hlist files are
dir_out    = '/project2/chihway/sims/UNIT/UNITSIMS_GADGET/fixedAmp_002/cib/'
chunk_size = 5000000 #size of chunks relies on your available memory #---> 323 iterations for 1 block

#--------------------------------------------------

h          = 0.6774
pars      = camb.CAMBparams()
pars.InitPower.set_params(ns=0.9667,As=2.142e-9,pivot_scalar=0.05)
pars.set_cosmology(H0=67.74, ombh2=0.02230 , omch2=0.1188, tau=0.066 )
pars.set_for_lmax(2000, lens_potential_accuracy=1)
results   = camb.get_results(pars)

zmid      = results.redshift_at_comoving_radial_distance(0.5*(chilow+chiupp)/h)
#chi    = results.comoving_radial_distance(2)*h

def getnearestsnap(dir_hlist,zmid):
    import glob 
    files = glob.glob(dir_hlist+"*.bz2")
    asnap = [np.float(line.rstrip('\n')[-16:-9]) for line in files]
    zsnap = 1./np.asarray(asnap)-1.
    return files[np.argmin(np.abs(zsnap-zmid))]

hlist_file = getnearestsnap(dir_hlist,zmid)

# num of lines = 175357379
'''
reader     = pd.read_csv('hlist_0.09140.list',\
                         #hlist_file,\
	                     chunksize = chunk_size,\
	                     engine    = 'c',\
	                     delim_whitespace=True,\
	                     skiprows  = 65,\
	                     usecols   = [0,17,18,19,60,62],\
	                     names     = ['scale','px','py','pz','Mpeak','Vpeak']\
	                     #nrows     = 500000
	                     )    
'''
ret        = np.zeros(hp.nside2npix(nsideout))

ntiles = int(np.ceil(maxchi/boxL))

def checkslicehit(chilow,chihigh,xx,yy,zz):
    # doing pre-selection so that we're not loading non-intersecting blocks 
    bvx=np.array([0,boxL,boxL,0,0,boxL,boxL,0])
    bvy=np.array([0,0,boxL,boxL,0,0,boxL,boxL])
    bvz=np.array([0,0,0,0,boxL,boxL,boxL,boxL])
    sx  = (bvx - origin[0] + boxL * xx);
    sy  = (bvy - origin[1] + boxL * yy);
    sz  = (bvz - origin[2] + boxL * zz);
    r   = np.sqrt(sx*sx + sy*sy + sz*sz)
    if ( (np.all(chilow<r)) | (np.all(chihigh<r)) ):
        return False
    else:
        return True

reader     = pd.read_csv('hlist_0.33030.list_lite',\
                         #hlist_file,\
                         #chunksize = chunk_size,\
                         #dtype     = {'scale': np.float64, 'px': np.float64, 'py': np.float64, 'pz': np.float64, 'Mpeak': np.float64, 'Vpeak': np.float64},\
                         #engine    = 'c',\
                         #delim_whitespace=True,\
                         #skiprows  = 1,\
                         #usecols   = [0,1,2,3,4,5],\
                         #names     = ['scale','px','py','pz','Mpeak','Vpeak'],\
                         #low_memory=False\
                         #nrows     = 500000
                         )   

ppx   = reader['px'].values
ppy   = reader['py'].values
ppz   = reader['pz'].values
a     = reader['scale'].values
Mpeak = reader['Mpeak'].values
Vpeak = reader['Vpeak'].values

z     = 1./a-1.
del a

ret        = np.zeros(hp.nside2npix(nsideout))

totslicehit=0
for xx in range(-ntiles,ntiles+1):
    for yy in range(-ntiles,ntiles+1):
        for zz in range(-ntiles,ntiles+1):

            print("%d %d %d"%(xx,yy,zz))
        
            slicehit = checkslicehit(chilow,chiupp,xx,yy,zz)
            #slicehit = True
            if slicehit==True:
                totslicehit+=1
                print('slicehit')

                for chunk in reader:
                    sx  = ppx - origin[0] + boxL * xx
                    sy  = ppy - origin[1] + boxL * yy
                    sz  = ppz - origin[2] + boxL * zz
                    r   = np.sqrt(sx*sx + sy*sy + sz*sz)
                    idx = np.where((r>chilow) & (r<chiupp))[0]

                    pix     = hp.vec2pix(nsideout,sx[idx]/r[idx],sy[idx]/r[idx],sz[idx]/r[idx])
                    IRflux  = um.halo2irflux(z,Mpeak,vMpeak)
                    q       = np.arange(np.amax(pix)+1)
                    ret[q]  = ret[q] + np.bincount(pix, weights=IRflux)

print(totslicehit)
hp.write_map(dir_out+'/cibflux_%.5f.fits'%zmid,ret,overwrite=True)




