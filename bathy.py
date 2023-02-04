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

colors = np.array([ [0.65, 0.60, 0.60], \
                      [0.65, 0.75, 0.75], \
                      [0.65, 0.90, 0.90], \
                      [0.60, 1.00, 1.00], \
                      [0.50, 1.00, 1.00], \
                      [0.40, 1.00, 1.00], \
                      [0.25, 1.00, 1.00], \
                      [0.10, 1.00, 1.00], \
                      [0.05, 0.95, 1.00], \
                      [0.00, 0.85, 1.00], \
                      [0.00, 0.75, 1.00], \
                      [0.00, 0.65, 1.00], \
                      [0.00, 0.50, 1.00], \
                      [0.00, 0.40, 0.85], \
                      [0.00, 0.30, 0.75], \
                      [0.00, 0.20, 0.60], \
                      [0.00, 0.00, 0.50]], 'f')
cmap = mpl.colors.ListedColormap(colors)

levels = [ 0, 10, 20, 30, 40, 50, 60, 100, 200, 400, 600, 800, \
           1000, 1500, 2000, 2500, 3000, 4000, 5000 ]

fig = plt.figure(figsize=(9,8))
#ax1 = fig.add_subplot(111)

print("Plotting Depth")
m.drawcoastlines()
m.drawmapboundary()
#m.drawmapboundary(fill_color='#99ffff')
m.fillcontinents(color='0.3',lake_color='0.3')
#   m.fillcontinents(color='coral',lake_color='aqua')

cs = m.contourf(x, y, grd.vgrid.h, levels=levels, cmap=cmap, extend='both')
# labels = [left,right,top,bottom]
parallels = np.arange(10.,90.,10.)
m.drawparallels(parallels,labels=[0, 1, 1, 0])
meridians = np.arange(10.,361.,10.)
m.drawmeridians(meridians,labels=[1, 0, 0, 1])

cbar = plt.colorbar(cs, orientation='vertical', pad=0.2, shrink=0.6)
#cbar.ax.tick_params(labelsize=15)

plt.savefig('bathy_Arctic6.png')
plt.close()
