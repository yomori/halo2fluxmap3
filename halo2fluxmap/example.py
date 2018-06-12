#!/usr/bin/env python
from halo2fluxmap import *
import default_params as params
import mpi4py.rc
mpi4py.rc.initialize = False
from mpi4py import MPI

params = getparams('params.py')

if MPI.Is_initialized():
    params.comm     = MPI.COMM_WORLD
    params.rank     = params.comm.Get_rank()
    params.size     = params.comm.Get_size()
    params.parallel = True
else:
    params.rank     = 0
    params.size     = 1
    params.parallel = False
    
print params.rank, params.size

fmt       = '%H:%M:%S on %m/%d/%Y'
timestamp = datetime.datetime.now().strftime(fmt)

if(params.rank==0):
    bar = 72*'-'
    print ''
    print bar
    print 'hod2map running on',params.size,'processor(s)'
    print 'Time:      '+timestamp
    print 'Directory: '+os.getcwd()
    print bar
    print ''

params.proc     = psutil.Process(os.getpid())
params.overhead = params.proc.memory_info().rss
params.maxmem   = 0

params.justify = 25*' '

Non    = get_catalog_Non()
Nreads = max(params.Nreads,params.size/2)

report('Reading, Shuffling, Culling Catalogue',2)
#To load balance without reading in entire catalogue on each processor
for i in range(Nreads):  
    ceni, Mmax = get_catalog(i,Nreads)
    ceni       = cull_catalog(ceni)
    ceni       = shuffle_catalog(ceni)
    if i==0:
        cen = distribute_catalog(ceni)
    else:
        cen  = np.concatenate((cen,distribute_catalog(ceni)))

    report('  '+str(i+1)+' of '+str(Nreads)+' chunks read and culled, Ncen = '+str(np.shape(cen)[0]),2)

del ceni

#Center map on most massive halo
#cen = center_catalog(cen,Mmax)

#Get number of satellites for each halo
report('cen2ns starting',2)
ns  = cen2ns(cen)
report('cen2ns done',2)

#loop over redshifts
for izslic in range(params.num_redshift):
    dz_slice = (params.max_redshift-params.min_redshift) / params.num_redshift
    zmin     = params.min_redshift + izslic * dz_slice
    zmax     = params.min_redshift + (izslic+1) * dz_slice

    if(params.rank==0): print '\n\tStarting z slice %f < z < %f \n' % (zmin,zmax) 
    #max z cut
    dm = [ r2z(np.sqrt(cen[:,0]**2+cen[:,1]**2+cen[:,2]**2)) < zmax]
    cen_i = cen[dm]
    ns_i  = ns[dm]
    #min z cut
    dm = [ r2z(np.sqrt(cen_i[:,0]**2+cen_i[:,1]**2+cen_i[:,2]**2)) > zmin ]
    cen_i = cen_i[dm]
    ns_i  = ns_i[dm]
    
    #Populate remaining halos with satellites
    report('populate_centrals starting',2)
    sat    = populate_centrals(cen_i,ns_i)
    report('populate_centrals done',2)

    #Write timing
    if params.num_redshift==1: write_time("HOD completed", params.rank)

    #loop over frequencies
    for inu_map in range(len(params.freq_list)):

        #Put halos in map
        params.nu_obs     = params.freq_list[inu_map] * 1.e9
        params.nu_obs_GHz = params.freq_list[inu_map]
        
        #Give flux to gals and put in map
        cib_tot, flux_cen, flux_sat = halos2map(cen_i,ns_i,sat) 
        
        if params.parallel: params.comm.Barrier()

        #SAVE MAP TO DISK
        ns_str   = str(params.nside)
        nu_str   = str(int(params.freq_list[inu_map]))
        zmin_str = str("%.2f" % zmin)
        zmax_str = str("%.2f" % zmax)

        #creating output file names
        dirout = "maps/" 
        if params.numdens==0: base = dirout+"cib_"
        if params.numdens==1: base = dirout+"opt_"

        if params.flat==0: base += "fullsky_"
        if params.flat==1: base += "flatsky_"

        if params.numdens==0: base += ('ns'+ns_str+'_zmin'+zmin_str+
                                       '_zmax'+zmax_str+'_nu'+nu_str)

        if params.numdens==0: base += '_ns'+ns_str
        if params.numdens==1: base += 'ns'+ns_str+'_zmin'+zmin_str+'_zmax'+zmax_str

        if(params.rank==0): 

            if params.flat==0:
                file_tot = base+'_tot_hp.fits'
                hp.write_map(file_tot,cib_tot)

            elif params.flat==1:
                file_tot = base+'_tot_fs.map'
                fov      = np.array([params.fov*np.pi/180]).astype('float32')
                npix     = np.array([params.nside]).astype('int32')
                cib_tot  = cib_tot.astype('float32')

                outfile = open(file_tot,'wb') 
                npix.tofile(outfile)
                npix.tofile(outfile)
                fov.tofile(outfile)
                fov.tofile(outfile)
                cib_tot.tofile(outfile)

            print params.justify,'Wrote Map nu = '+nu_str+'GHz'

    #Write timing
write_time("Halo projection completed", params.rank)


