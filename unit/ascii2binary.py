import pandas as pd
import sys
import numpy as np
a = float(sys.argv[1])

dir_hlist='/project2/chihway/sims/MDPL2/hlists/'

reader     = pd.read_csv(dir_hlist+'hlist_%.5f.list'%a,\
                         #hlist_file,\
                         #chunksize = 5000000,\
                         engine    = 'c',\
                         delim_whitespace=True,\
                         skiprows  = 65,\
                         usecols   = [0,17,18,19,20,21,22,39],\
                         names     = ['scale','px','py','pz','vx','vy','vz','M200c'],\
                         #nrows     = 5000000
                         )

np.save('/project2/chihway/sims/MDPL2/hlists/stripped_%.6f.npy'%a,\
        np.c_[(reader['px'].values).astype(np.float32),
              (reader['py'].values).astype(np.float32),
              (reader['pz'].values).astype(np.float32),
              (reader['vx'].values).astype(np.float32),
              (reader['vy'].values).astype(np.float32),
              (reader['vz'].values).astype(np.float32),
              (reader['M200c'].values).astype(np.float32)
             ]
        ) 
#i=0
#for chunk in reader:
#    if i==0:
#        chunk.to_csv(dir_hlist+'hlist_%.5f.list_lite'%a,mode='w')
#    else:
#        chunk.to_csv(dir_hlist+'hlist_%.5f.list_lite'%a,mode='a',header=False)
#    i+=1
