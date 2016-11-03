import pycnal
import pycnal_toolbox
from mpl_toolkits.basemap import Basemap, shiftgrid
import netCDF4

# BE SURE TO RUN sphere FIRST!!! (sets spherical to true)

grd = pycnal.grid.get_ROMS_grid('BENGUELA')
grd = pycnal.grid.get_ROMS_grid('S_AFRICA')

# dir(grd), dir(grd.hgrid)
#
# Resolution values are 'c', 'l', 'i', 'h', and 'f'.
# S_Africa
#m = Basemap(projection='lcc',lat_1=-35, lat_2=-10, lon_0=25, lat_0=-20, width=7 000000, height=7000000, resolution='h')
# Benguela
m = Basemap(projection='lcc',lat_1=-35, lat_2=-10, lon_0=10, lat_0=-13, width=4000000, height=5000000, resolution='h')

coast = pycnal.utility.get_coast_from_map(m)
pycnal.grid.edit_mask_mesh_ij(grd.hgrid, coast=coast)

#pyroms.grid.edit_mask_mesh(grd.hgrid, proj=m)

pycnal.grid.write_ROMS_grid(grd, filename='grid_py.nc')
# ncks -v mask_rho,mask_u,mask_v,mask_psi grid_py.nc grid_SW_Africa.nc
# ncks -v mask_rho,mask_u,mask_v,mask_psi grid_py.nc grid_S_Africa.nc


