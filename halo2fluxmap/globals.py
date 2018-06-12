import default_parameters as params 
import numpy as np
                         # constants in cgs 
h     = 6.62606957e-27   # erg.s
c     = 3e+10            # cm/s
k     = 1.3806488e-16    # erg/K
Msun  = 2e33             # g
Mpc   = 3.086e24         # cm

                         # fiducial values in cgs 
#nuf   = params.nu_obs   # 150 GHz
Rf    = 1e3              # Fiducial radius in Mpc
Rfg   = Rf*Mpc           # Fiducial radius in cm
Mf    = 1                # Fiducial mass value in Msun
Mfg   = Mf*Msun          # Fiducial mass value in g
Nsf   = params.nside     # Nside value

                         # conversion factors from 
                         # dimensionless to cgs units
#ft = c**2/(2*k*nuf**2)   # cm^2/K/erg
ff = 1/Rfg**2/4/np.pi    # 1/cm^2
fl = params.shang_L0        # erg/sec/Hz/g
fs = Mfg                 # g

                         # final conversion factor from 
                         # dimensionless temperature to muK
#T0 = ft*fi*ff*fl*fs*1e6 # muK = cm^2/K/erg*1/cm^2*erg/sec/Hz/g*g*1e6

# convert model parameters to fiducial units
params.shang_Mpeak /= Mf

