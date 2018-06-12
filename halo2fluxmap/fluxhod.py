from mods               import *  
from hod                import *
from globals import *

def printbounds(a):
	print a.min(),a.max()

def sigma_cen(m):

	arg  = np.log10(m)
	arg -= np.log10(params.shang_Mpeak)
	arg  = arg**2
	arg /= 2*params.shang_sigmaM
	arg  = np.exp(-arg)
	arg *= m
	arg *= 1/np.sqrt(2*np.pi*params.shang_sigmaM)

	return arg

def sigma_sat(m):

	# Make table of <L_sat> vs. M_halo
	x = np.linspace(np.log(params.shang_Mmin),np.log(m.max()),100)
	L_mean     = np.zeros(len(x))
	for i in range(len(x)):
		L_mean[i],err = quad(integrand_L,np.log(params.shang_Mmin),
				     x[i],args=x[i])
	f_L = interpolate.interp1d(x,L_mean)

	return f_L(np.log(m))

def jiang_shmf(m,M_halo):	
	gamma_1    = 0.13
	alpha_1    = -0.83
	gamma_2    = 1.33
	alpha_2    = -0.02
	beta_2     = 5.67
	zeta       = 1.19 	
	
	dndm = (((gamma_1*((m/M_halo)**alpha_1))+
		 (gamma_2*((m/M_halo)**alpha_2)))*
		(np.exp(-(beta_2)*((m/M_halo)**zeta))))

	return dndm

# <L_sat> interpolation values
def integrand_L(lm,lM_halo):
	m      = np.exp(lm)
	M_halo = np.exp(lM_halo)
	
	dns_dm = jiang_shmf(m,M_halo)
	dns_dm_sigma = sigma_cen(m) * dns_dm
	
	return dns_dm_sigma
	
def f2t(intensity):
	x = h*params.nu_obs/k/2.726
	T = ((np.exp(x)-1)**2/np.exp(x)*intensity).astype('float32')

	return T.astype('float32')

def nu2theta(nu):
	xnu     = h*nu/k/params.shang_Td
	Thetanu = (1/(nu*params.shang_I0)*
		   xnu**(4.+params.shang_beta)/(np.exp(xnu)-1.))
	return Thetanu

def LF(M,x,y,z,gtype):

	r  = x**2 
	r += y**2
	r += z**2
	r  = r**0.5
	r *= Rf	

	z  = r2z(r)
	r  = (1+z)*params.nu_obs      #observed frequency in Hz

	if (params.LM=="Planck2013"):
		if(gtype == 'cen'): L = sigma_cen(M)
		if(gtype == 'sat'): L = sigma_sat(M)
	if (params.LM=="Planck2015"):              #const * M_500/1e14M_sun
		L = M*np.sqrt(200./500) / 1.e14 #sqrt(200/500) to convert M_200 to M_500

	L *= nu2theta(r)
	L *= (1+z)**params.shang_eta

	return L

def l2f(L,x,y,z):	

	r  = x**2 
	r += y**2
	r += z**2
	r  = r**0.5
	r *= Rf	

	return L / r**2 / (1+r2z(r)) * Rf**2

def dimensionless(cen,sat,i):

	cen[:,3]   *= Mf**i

	cen[:,0:3] *= Rf**i
	sat[:,0:3] *= Rf**i
	
	return cen, sat

def halos2map(cen,ns,sat):
	import time

	if params.flat==0: 
		fi = 12.*Nsf**2/4/np.pi               # dimensionless
	elif params.flat==1:
		fi = Nsf**2 / np.radians(params.fov)**2  # dimensionless


	t2=time.time()

	N = np.shape(sat)[0]

	# Convert to dimensionless units
	cen, sat = dimensionless(cen,sat,-1)

	C = np.zeros((np.shape(cen)[0],1),dtype='float32')
	S = np.zeros((np.shape(sat)[0],2),dtype='float32')
	
	C[:,0] = cen[:,3] # Mass
	S[:,:] = flux2map.cen2sat(np.column_stack((C[:,0],ns)),ns) # Parent Mass and Nsat

	#CALCULATE FLUX FROM CENTRALS AND SATELLITES
	if params.numdens==0:
		if (params.LM=="Planck2013"):
			C[:,0] = LF(C[:,0],cen[:,0],cen[:,1],cen[:,2],'cen')
			S[:,0] = LF(S[:,0],sat[:,0],sat[:,1],sat[:,2],'sat') / S[:,1]
		if (params.LM=="Planck2015"):
		#divide total luminosity evenly between centrals and satellites
			C[:,0] = LF(C[:,0],cen[:,0],cen[:,1],cen[:,2],'cen') / (ns+1)
			S[:,0] = LF(S[:,0],sat[:,0],sat[:,1],sat[:,2],'sat') / (S[:,1]+1)
		
		C[:,0] = l2f(C[:,0],cen[:,0],cen[:,1],cen[:,2])
		S[:,0] = l2f(S[:,0],sat[:,0],sat[:,1],sat[:,2])

	elif params.numdens==1:
		C[:,0] = 1.0
		S[:,0] = 1.0

	cen_fluxes = C[:,0]
	sat_fluxes = S[:,0]

	total_fluxl=np.zeros(1)
	total_flux=np.zeros(1)

	if params.flat==1:				
		thetaxc = np.abs(np.arctan(cen[:,0]/cen[:,2]))*2
		thetayc = np.abs(np.arctan(cen[:,1]/cen[:,2]))*2	
		dmc = [(thetaxc < np.radians(params.fov)) & (thetayc < np.radians(params.fov))
		       & (cen[:,2]>0)]

		thetaxs = np.abs(np.arctan(sat[:,0]/sat[:,2]))*2
		thetays = np.abs(np.arctan(sat[:,1]/sat[:,2]))*2	
		dms = [(thetaxs < np.radians(params.fov)) & (thetays < np.radians(params.fov))
		       & (sat[:,2]>0)]

		total_fluxl = cen_fluxes[dmc].sum() + sat_fluxes[dms].sum() # fiducial units
	else:
		total_fluxl = cen_fluxes.sum() + sat_fluxes.sum() # fiducial units
		
	total_fluxl = total_fluxl * ff*fl*fs       # erg/sec/cm^2/Hz

	#MAKE MAP
	#To save memory no longer makes T_cen,T_sat,T_tot maps
	if params.flat==0:
		omega_map = 4.*np.pi		
		T_tot  = flux2map.makemap(
			cen[:,0],cen[:,1],cen[:,2],C[:,0],params.nside)
		T_tot += flux2map.makemap(
			sat[:,0],sat[:,1],sat[:,2],S[:,0],params.nside)
	elif params.flat==1:
		omega_map = np.radians(params.fov)**2
		T_tot  = flux2map.makemapflat(
			cen[:,0],cen[:,1],cen[:,2],C[:,0],params.nside,params.fov)
		T_tot += flux2map.makemapflat(
			sat[:,0],sat[:,1],sat[:,2],S[:,0],params.nside,params.fov)

	# Convert back to physical units
	cen, sat = dimensionless(cen,sat,1)
	# Convert from dimensionless temperature to muK	
	T_tot  = T_tot.astype("float32")

	if params.numdens==0:
		T_tot *= fi*ff*fl*fs

	# Sum over all processes
	T = np.zeros(np.shape(T_tot)[0],dtype='float32')

        if params.parallel:
                from mpi4py import MPI 
                params.comm.Reduce([T_tot,MPI.FLOAT],[T,MPI.FLOAT],
                                   op = MPI.SUM,root=0)
                params.comm.Reduce(total_fluxl,total_flux,op=MPI.SUM,root=0)
        else:
                T          = T_tot
                total_flux = total_fluxl

	if(params.rank==0):
		# normalize final map values to specified frequency
		if params.nu_obs_GHz  == params.norm_freq:
			params.norm_value = np.mean(T)

		T = T/params.norm_value
		T = T*params.Inu_norm

                # get total flux of map by multiplying mean intensity by map solid angle
		map_flux = (T).mean()*omega_map

                # normalize sources to have the same total flux as map
                tot_source_flux = cen_fluxes.sum() + sat_fluxes.sum()
                cen_fluxes *= map_flux / tot_source_flux
                sat_fluxes *= map_flux / tot_source_flux
                tot_source_flux = cen_fluxes.sum() + sat_fluxes.sum()

                report('',2)
                report('Finished making map nu:     '+str(params.nu_obs_GHz)+' GHz',2)
                report('     mean source intensity: '+str((cen_fluxes.sum() + sat_fluxes.sum()) / omega_map),2)
                report('        mean map intensity: '+str(T.mean()),2)
                report('    mean central intensity: '+str(cen_fluxes.sum() / omega_map),2)
                report('  mean satellite intensity: '+str(sat_fluxes.sum() / omega_map),2)
                report('           min, max of map: '+str(np.min(T))+', '+str(np.max(T)),2)

		nonzeropix = np.shape(T[T>0])[0]
		
                if params.iwantcat==1:   #save centrals and satellites
			np.savez(params.inputfile+"_censat_"+str(params.nu_obs_GHz),x_cen=cen[:,0],y_cen=cen[:,1],z_cen=cen[:,2],
				 M_cen=cen[:,3],ns_cen=ns,x_sat=sat[:,0],y_sat=sat[:,1],z_sat=sat[:,2],
				 M_sat=S[:,0],cen_fluxes=cen_fluxes,sat_fluxes=sat_fluxes)                        

		t1=time.time()
		dt = (t1-t2)/60.
		ng = np.shape(sat)[0]+np.shape(cen)[0]
		if(params.rank==0 and params.verbose>0 and 
		   params.nu_obs_GHz == params.norm_freq):
			report('Time to project galaxies: '+str(dt)+' minutes',2)
			report('Number of nonzero pixels:   '+str(nonzeropix),-1)
		return T, cen_fluxes, sat_fluxes
	else:
		return 0, 0, 0 
