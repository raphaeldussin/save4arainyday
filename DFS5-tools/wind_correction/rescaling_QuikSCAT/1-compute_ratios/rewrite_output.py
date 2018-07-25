import PyRaf
import numpy

lon = PyRaf.readfull('./output.nc','lon')
lat = PyRaf.readfull('./output.nc','lat')
#lon2,lat2 = numpy.meshgrid(lon,lat)
ratio = PyRaf.readfull('./output.nc','ratio')

spval = ratio.min()
ratio[numpy.where(ratio == spval)] = 0.
# in order to remove drown values in the arctic
ratio[numpy.where(lat >= 69. )] = 0.


PyRaf.write_2d_nemo_file('./extended_ouput.nc',lon,lat,0.,ratio,'ratio')
