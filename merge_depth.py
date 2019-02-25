import numpy as np
import netCDF4
import sys

ncfile = sys.argv[1]
nc = netCDF4.Dataset(ncfile, 'a', format='NETCDF3_64BIT')

hraw = nc.variables['hraw'][:]
h = hraw[3,:,:]

h = np.where(hraw[1,:,:]==0, h, hraw[1,:,:])

nc.variables['hraw'][4,:,:] = h

nc.close()
