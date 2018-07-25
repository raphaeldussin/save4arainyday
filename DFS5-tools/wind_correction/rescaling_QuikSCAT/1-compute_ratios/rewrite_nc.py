import PyRaf
import numpy

lon = PyRaf.readfull('./ratio_QuikSCAT_on_ERAi_2000-2008.nc','lon')
lat = PyRaf.readfull('./ratio_QuikSCAT_on_ERAi_2000-2008.nc','lat')
lon2,lat2 = numpy.meshgrid(lon,lat)
ratio = PyRaf.readfull('./ratio_QuikSCAT_on_ERAi_2000-2008.nc','U10')

spval = ratio.min()
ratio[numpy.where(ratio == spval)] = 0.

PyRaf.write_2d_nemo_file('./nemo_ratio_QuikSCAT_on_ERAi_2000-2008.nc',lon2,lat2,0.,ratio,'ratio')
