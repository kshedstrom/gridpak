import pyroms
import pyroms_toolbox
from mpl_toolkits.basemap import Basemap, shiftgrid
import netCDF4

grd = pyroms.grid.get_ROMS_grid('CHUKCHI5')

# dir(grd), dir(grd.hgrid)
#
# Resolution values are 'c', 'l', 'i', 'h', and 'f'.
map = Basemap(projection='npstere', boundinglat=60, lon_0=240, resolution='h')

pyroms.grid.edit_mask_mesh(grd.hgrid, proj=map)

pyroms.grid.write_ROMS_grid(grd, filename='grid_py.nc')
# ncks -A -v mask_rho,mask_u,mask_v,mask_psi grid_py.nc grid_Chukchi_5.nc
# ncks -A -v mask_rho,mask_u,mask_v,mask_psi grid_py.nc grid_Chukchi_4.nc
# ncks -A -v mask_rho,mask_u,mask_v,mask_psi grid_py.nc grid_Chukchi_3.nc


