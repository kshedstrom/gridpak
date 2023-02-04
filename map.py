from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import sys

ncfile = 'grid_Arctic_6.nc'
nc = netCDF4.Dataset(ncfile, 'a')

# setup Lambert Conformal basemap.
m = Basemap(projection='npstere', boundinglat=45, lon_0=160, resolution='h')

# draw a boundary around the map, fill the background.
# this background will end up being the ocean color, since
# the continents will be drawn on top.
m.drawmapboundary(fill_color='aqua')
# fill continents, set lake color same as ocean color.
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
# label parallels on right and top
# meridians on bottom and left
parallels = np.arange(40.,81,10.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(20.,361.,20.)
m.drawmeridians(meridians,labels=[True,False,False,True])
# plot blue dot somewhere and label it as such.
lon, lat = -150., 60.0
# convert to map projection coords.
# Note that lon,lat can be scalars, lists or numpy arrays.
xpt,ypt = m(lon,lat)
# convert back to lat/lon
lonpt, latpt = m(xpt,ypt,inverse=True)
m.plot(xpt,ypt,'bo')  # plot a blue dot there
# put some text next to the dot, offset a little bit
# (the offset is in map projection coordinates)
plt.text(xpt+100000,ypt+100000,'Spot (%5.1fW,%3.1fN)' % (lonpt,latpt))
plt.text(xpt-400000,ypt+400000,'Spot (%5.2f,%5.2f)' % (xpt*1.e-6,ypt*1.e-6))
plt.show()
