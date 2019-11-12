''' read a previous simulation, grab the final pressure result to use as a pressure condition'''

import h5py

f = h5py.File('pflotran.h5','r')
f2 = h5py.File('ic_pressure.h5','w')

p = f[f.keys()[-2]]['Liquid Pressure [Pa]'].value

f2['Cell Ids']=range(0,len(p))


f.close()
f2.close()
