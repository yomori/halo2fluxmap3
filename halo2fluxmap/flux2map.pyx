import healpy as hp
import numpy as np
cimport numpy as np

def makemap(np.ndarray x,np.ndarray y, np.ndarray z,np.ndarray F, nside):
    
    pixind = hp.vec2pix(nside,x,y,z) #to get same orientation as pks2map
    hmap = np.zeros(hp.nside2npix(nside))
    for i in range(np.shape(F)[0]):
        hmap[pixind[i]] += F[i]

    return hmap

def makemapflat(np.ndarray x,np.ndarray y, np.ndarray z,np.ndarray F, nside,fov):
        fov = fov*np.pi/180           
        dp  = fov/nside  
        #+fov/2 to ensure numbers in pixind are positive

        dm = [(z>0)]   
        z = z[dm]
        y = y[dm]         
        x = x[dm] 
        F = F[dm]

        r = np.sqrt(z**2+x**2+y**2)
        y = y/r
        x = x/r
        y = np.arcsin(y) #now theta
        x = np.arcsin(x) #now phi
#        y = np.arctan(y/z)# + fov/2 #now theta
#        x = np.arctan(x/z)# + fov/2 #now phi

        dm = [(abs(y)<fov/2)]   
        y = y[dm]         
        x = x[dm] 
        F = F[dm]
        dm = [(abs(x)<fov/2)] 
        y = y[dm] 
        x = x[dm]
        F = F[dm]

        pixind = np.floor((x+fov/2)/dp) + nside*np.floor((y+fov/2)/dp) 
        pixind = pixind.astype(int)
        
        hmap = np.zeros(nside**2)
        for i in range(np.shape(F)[0]):
                hmap[pixind[i]]+=F[i] 
        return hmap
    
def cen2sat(np.ndarray cen, np.ndarray n):
    
    cdef int N_sat = np.sum(n)
    cdef int N_cen = np.shape(cen)[0]
    cdef int N_prp = np.shape(cen)[1]
    
    sat = np.zeros((N_sat,N_prp),dtype='float32')
    
    cdef int count = 0
    cdef int i
    for i in range(N_cen):
        sat[count:count+n[i],:] = cen[i,:]
        count += n[i]	
        
    return sat
    
