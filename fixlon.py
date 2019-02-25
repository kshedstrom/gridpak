import numpy as np
import netCDF4
import sys

ncfile = sys.argv[1]
nc = netCDF4.Dataset(ncfile, 'a', format='NETCDF3_64BIT')

lon_rho = nc.variables['lon_rho'][:]
lon_psi = nc.variables['lon_psi'][:]
lon_u = nc.variables['lon_u'][:]
lon_v = nc.variables['lon_v'][:]

lon_rho = np.where(lon_rho>=0, lon_rho, lon_rho+360.0)
lon_psi = np.where(lon_psi>=0, lon_psi, lon_psi+360.0)
lon_u = np.where(lon_u>=0, lon_u, lon_u+360.0)
lon_v = np.where(lon_v>=0, lon_v, lon_v+360.0)

nc.variables['lon_rho'][:] = lon_rho
nc.variables['lon_psi'][:] = lon_psi
nc.variables['lon_u'][:] = lon_u
nc.variables['lon_v'][:] = lon_v

nc.close()
