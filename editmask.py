import pyroms
import pyroms_toolbox
from mpl_toolkits.basemap import Basemap, shiftgrid
import netCDF4

grd = pyroms.grid.get_ROMS_grid('ARCTIC6')

# dir(grd), dir(grd.hgrid)
#
# Resolution values are 'c', 'l', 'i', 'h', and 'f'.
map = Basemap(projection='npstere', boundinglat=45, lon_0=160, resolution='f')

pyroms.grid.edit_mask_mesh(grd.hgrid, proj=map)

pyroms.grid.write_ROMS_grid(grd, filename='grid_py.nc')
# ncks -A -v mask_rho,mask_u,mask_v,mask_psi grid_py.nc grid_Arctic_6.nc


