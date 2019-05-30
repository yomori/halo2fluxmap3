import numpy as np
import camb 
from camb import model, initialpower
from tabulate import tabulate

h          = 0.6774

zmax       = 10.0   # maximum redshift to make the halo lightcone
boxL       = 1000   # length of simulation box in Mpc/h
dL         = 50     # thickness of shell to use in Mpc/h
nside      = 8192   # output healpix map resolution
list_hlist = '/project2/chihway/sims/UNIT/UNITSIMS_GADGET/fixedAmp_001/ROCKSTAR_HALOS/unitsims.ft.uam.es/DATABASE/UNITSIMS_GADGET/ROCKSTAR_HALOS/fixedAmp_001/hlist/list_hlist'
LightConeOrigin = np.array([0,0,0]) # where to place the origin of lightcone

#-------------------------------------------------
# Doing comoving distance calculations
pars      = camb.CAMBparams()
pars.set_cosmology(H0=67.74, ombh2=0.02230 , omch2=0.1188, tau=0.066 )
pars.InitPower.set_params(ns=0.9667,As=2.142e-9,pivot_scalar=0.05)
pars.set_for_lmax(2000, lens_potential_accuracy=1)
results   = camb.get_results(pars)

chimax    = results.comoving_radial_distance(zmax)*h
Ntess     = int(np.ceil(chimax/boxL))
Nslices   = int(np.ceil(chimax/dL))

chi_edges = (np.arange(Nslices+1)*dL)*h
#-------------------------------------------------

print(tabulate([['Maximum chi [Mpc/h]'    , chimax],
                ['Box size [Mpc/h]'       , boxL],
                ['Number of tesselations' , Ntess],
                ['Number of slices'       , Nslices],
                ['Lightcone origin'       , LightConeOrigin]],
                headers=['Parameter', 'Value']))

print("---------------------------------------")

hlists = [line.rstrip('\n') for line in open(list_hlist)]
Nsnaps = len(hlists)
alist  = np.zeros(Nsnaps)
zlist  = np.zeros(Nsnaps)
for i in range(0,Nsnaps):
    tmp = hlists[i]
    alist[i] = float(tmp[tmp.index("_")+1:-5])
    zlist[i] = float(1/(alist[i])-1.)
    print("%s a=%.5f z=%.5f"%(hlists[i],alist[i],zlist[i]))

halo_phi = np.array([])
halo_tht = np.array([])

#if(np.max(zlist<zmax)):
#    raise Exception("Halo catalog didn't reach requested zmax")

for si in range(0,Nslices):

    chilow  = chi_edges[si]   # lower end of zbin
    chihigh = chi_edges[si+1] # upper end of zbin

    print("processing zslice %d/%d, chi range [%.2f/%.2f]"%(si+1,Nslices,chilow,chihigh))

    px,py,pz,Mh = np.loadtxt(hlists[i],unpack=True,skiprows=100)

    ret         = np.zeros(hp.nside2npix(nside))

    if chilow<chimax:

        print("Adding halos")

        for ti in range(-Ntess,Ntess+1):
            for tj in range(-Ntess,Ntess+1):
                for tk in range(-Ntess,Ntess+1):
                    vec      = np.zeros((Mh.shape[0],3))
                    vec[:,0] = (px - LightConeOrigin[0] + boxL * ti);
                    vec[:,1] = (py - LightConeOrigin[1] + boxL * tj);
                    vec[:,2] = (pz - LightConeOrigin[2] + boxL * tk);
                    r        = np.sqrt(vec[:,0]*vec[:,0] + vec[:,1]*vec[:,1] + vec[:,2]*vec[:,2])
                    vec[:,0]  /= r
                    vec[:,1]  /= r
                    vec[:,2]  /= r
                    tht,phi  = hp.vec2ang(vec)
                    #ra,dec   = tp2rd(tht,phi)

                    halo_tht = np.concatenate([halo_tht,tht])
                    halo_phi = np.concatenate([halo_phi,phi])
        pix = hp.ang2pix(nside,halo_tht,halo_phi)
        ret[pix]=1 