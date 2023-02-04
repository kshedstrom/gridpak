import netCDF4
import subprocess
import pyroms
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import datetime

grd = pyroms.grid.get_ROMS_grid('ARCTIC6')
clat = grd.hgrid.lat_rho
clon = grd.hgrid.lon_rho

m = Basemap(projection='npstere', boundinglat=45, lon_0=160, resolution='h')
x, y = m(clon, clat)
cmap = plt.cm.get_cmap("magma")

print("Plotting Size")
m.drawcoastlines()
m.drawmapboundary()
#m.drawmapboundary(fill_color='#99ffff')
m.fillcontinents(color='0.3',lake_color='0.3')
#   m.fillcontinents(color='coral',lake_color='aqua')

dx = grd.hgrid.dx
dy = grd.hgrid.dy
size = np.sqrt(dx**2 + dy**2)/np.sqrt(2.0)
cs = m.contourf(x, y, size, cmap=cmap, extend='both')
# labels = [left,right,top,bottom]
parallels = np.arange(10.,90.,10.)
m.drawparallels(parallels,labels=[0, 1, 1, 0])
meridians = np.arange(10.,361.,10.)
m.drawmeridians(meridians,labels=[1, 0, 0, 1])

cbar = plt.colorbar(cs, orientation='vertical', pad=0.2, shrink=0.6)
#cbar.ax.tick_params(labelsize=15)

plt.savefig('Arctic6_size.png')
plt.close()
