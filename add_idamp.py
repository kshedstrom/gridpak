import numpy as np
import sys
from subprocess import call
import xarray as xr


output_file = 'Arctic6_idamp.nc'

compress = True
if compress:
   fileformat, zlib, complevel, area_dtype = 'NETCDF4', True, 1, 'f4'
else:
   fileformat, zlib, complevel, area_dtype = 'NETCDF3_64BIT_OFFSET', None, None, 'd'

NJ = 1392
NI = 1080
width = 10
edge_val = 1.e-6

idamp = np.zeros((NJ,NI))
idamp_u = np.zeros((NJ,NI+1))
idamp_v = np.zeros((NJ+1,NI))

for i in range(width,-1,-1):
    idamp[:, i] = (width-i) * edge_val/width
    idamp[:, -i] = (width-i) * edge_val/width
    idamp[i, :] = (width-i) * edge_val/width
    idamp[-i, :] = (width-i) * edge_val/width
    idamp_u[:, i] = (width-i) * edge_val/width
    idamp_u[:, -i] = (width-i) * edge_val/width
    idamp_u[i, :] = (width-i) * edge_val/width
    idamp_u[-i, :] = (width-i) * edge_val/width
    idamp_v[:, i] = (width-i) * edge_val/width
    idamp_v[:, -i] = (width-i) * edge_val/width
    idamp_v[i, :] = (width-i) * edge_val/width
    idamp_v[-i, :] = (width-i) * edge_val/width

ida = xr.DataArray(idamp, dims=['yh', 'xh'])
ida.attrs["units"] = "s-1"
ida.attrs["cell_methods"] = "time: point"

idu = xr.DataArray(idamp_u, dims=['yh', 'xq'])
idu.attrs["units"] = "s-1"
idu.attrs["cell_methods"] = "time: point"

idv = xr.DataArray(idamp_v, dims=['yq', 'xh'])
idv.attrs["units"] = "s-1"
idv.attrs["cell_methods"] = "time: point"

base = xr.Dataset(data_vars = {"Idamp" : ida,
                               "Idamp_u" : idu,
			       "Idamp_v" : idv})
base.attrs["grid_type"] = "regional MOM6"
base.attrs["grid"] = "Arctic6"

comp = dict(zlib=zlib, complevel=complevel)
encoding = {var: comp for var in base.data_vars}
base.to_netcdf(
    output_file,
    format=fileformat,
    encoding=encoding
)
