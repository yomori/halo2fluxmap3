import matplotlib.pyplot as plt
import numpy as np
from   scipy.ndimage.filters import gaussian_filter

mapfile = open('maps/example_ns4096_nu545.map')

npix = np.fromfile(mapfile,count=2,dtype=np.int32)
fov  = np.fromfile(mapfile,count=2,dtype=np.float32)
flux = np.fromfile(mapfile,count=npix[0]*npix[1],dtype=np.float32)

spix = 1.0

flux = np.reshape(flux,(npix[0],npix[1]))
flux = gaussian_filter(flux,spix)

vmin = np.log10(flux[flux>0].min())

plt.imshow(np.log10(flux),vmin=vmin,cmap='binary')

plt.show()



