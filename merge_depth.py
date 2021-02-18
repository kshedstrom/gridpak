import numpy as np
import netCDF4
import sys

#ncfile = sys.argv[1]
ncfile = 'grid_Arctic_5.nc'
nc = netCDF4.Dataset(ncfile, 'a', format='NETCDF3_64BIT')

hraw = nc.variables['hraw'][:]
h = hraw[1,:,:]

h = np.where(hraw[5,:,:]==0, h, hraw[5,:,:])   #ARDEM
h = np.where(hraw[10,:,:]==0, h, hraw[10,:,:]) #IBCAO4

nc.variables['hraw'][11,:,:] = h

nc.close()
