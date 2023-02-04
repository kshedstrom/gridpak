import numpy as np
import netCDF4
import sys

#ncfile = sys.argv[1]
ncfile = 'grid_Arctic_6.nc'
nc = netCDF4.Dataset(ncfile, 'a')

mask_rho = nc.variables['mask_rho'][:]
mask_2 = nc.variables['mask_2'][:]

mask_rho[884:973,384:500] = mask_2[884:973,384:500]

nc.variables['mask_rho'][:] = mask_rho

nc.close()
