# uses universe machine relations to go from Halo mass to IR luminosiry

def vMpeak2SFR(z,vMpeak):
	"""
	takes in vMpeak and estimates the SFR of star forming galaxies
	The specific forms are from Appendix H of 1806.07893
	"""
	a   = 1/(z+1)

	V   = 10**(2.060+(-0.274)*(1-a)+(0.696)*np.log(1+z)+(-0.099*z))
	vs  = vMpeak/V

    alpha = -6.783+(-0.373)*(1-a)+(2.483)*np.log(1+z)+(-0.278)*z
    beta  = -0.670+(-0.673)*(1-a)+0.705*z
    gamma = 10**(-1.313+4.315*(1-a)+(-0.969)*z)
    delta = 0.078

	SFR   = eps *( (vs**alpha+vs**beta)**-1+gamma*np.exp(-np.log10(v**2)/2/delta**2)  )

	return SFR

def Mpeak2Mstar(z,Mpeak):
	if par==obsAllAllExcl:
		eps0    = -1.505
		alpha0  =  1.998
		gamma0  = -0.738
		epsa    =  0.697
		alphaa  = -1.394 
		gammaa  = -1.697
		epslna  =  0.002
		alphalna= -1.175
		gammaz  = -0.573
		epsz    =  0.124
		alphaz  =  0.166
		M0      = 11.979
		beta0   =  0.512
		Ma      =  2.293
		betaa   = -0.181
		Mlna    =  2.393
		betaz   = -0.160
		Mz      = -0.380
		delta0  =  0.382
		
	gamma = 10**( gamma0+gammaa*(a-1)+gammaz*z )
	delta = delta0
	beta  = beta0+betaa*(a-1)+betaz*z
	alpha = alpha0+alpha1*(a-1)-alphalna*np.log(a)+alphaz*z
	eps   = eps0+epsa*(a-1)-epslna*eps*np.log(a)+epsz*z
	M1    = 10**( M0+Ma*(a-1)-Mlna*np.log(a)+Mz*z ) # in Msol

    x     = np.log10(Mpeak/M1)

    Mstar = M1*10**( eps-np.log10(10**(-alpha*x)+10**(-beta*x)) +gamma*np.exp(-0.5*(x/delta)**2) ) 

    return Mstar

def sfr2LIR(Mstar,SFR):
	rUV = 4.07*np.log(Mstar)-39.32
	LIR = SFR/1e-10*10**(0.4*rUV)/(1+10**(0.4*rUV))
	return LIR

def halo2LIR(px,py,pz,z,Mpeak,vMpeak):
	#parallelize this
	Mstar = Mpeak2Mstar(z,Mpeak)
	SFR   = vMpeak2SFR(z,vMpeak)
	LIR   = fr2LIR(Mstar,SFR)

