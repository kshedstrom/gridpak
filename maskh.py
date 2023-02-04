import numpy as np
import netCDF4
import sys

ncfile = sys.argv[1]
nc = netCDF4.Dataset(ncfile, 'a')

h = nc.variables['h'][:]
mask = nc.variables['mask_rho'][:]

spval = 1.e30
spval_arr = spval*np.ones(h.shape)

h_mask = np.where(mask==0, spval_arr, h)

#nc.createVariable('h_mask', 'f8', ('eta_rho', 'xi_rho'), fill_value=spval)
nc.variables['h_mask'][:] = h_mask

nc.close()
