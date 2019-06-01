import pandas as pd

reader     = pd.read_csv('/project2/chihway/sims/UNIT/UNITSIMS_GADGET/fixedAmp_002/hlist/hlist_0.33030.list',\
                         #hlist_file,\
                         chunksize = 5000000,\
                         engine    = 'c',\
                         delim_whitespace=True,\
                         skiprows  = 65,\
                         usecols   = [0,17,18,19,60,62],\
                         names     = ['scale','px','py','pz','Mpeak','Vpeak'],\
                         #nrows     = 5000000
                         )
i=0
for chunk in reader:
    if i==0:
        chunk.to_csv('/project2/chihway/sims/UNIT/UNITSIMS_GADGET/fixedAmp_002/hlist/hlist_0.33030.list_lite',mode='w')
    else:
        chunk.to_csv('/project2/chihway/sims/UNIT/UNITSIMS_GADGET/fixedAmp_002/hlist/hlist_0.33030.list_lite',mode='a',header=False)
    i+=1