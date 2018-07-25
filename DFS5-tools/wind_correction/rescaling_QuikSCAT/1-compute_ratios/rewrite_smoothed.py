import PyRaf
import numpy as npy

lon = PyRaf.readfull('./extended_ouput.ncB007','nav_lon')
lat = PyRaf.readfull('./extended_ouput.ncB007','nav_lat')

ratio = PyRaf.readfull('./extended_ouput.ncB007','ratio')

lonout   = npy.concatenate((lon[:,512:622]-360.,lon[:,110:512]),axis=1)
latout   = npy.concatenate((lat[:,512:622]     ,lat[:,110:512]),axis=1)
ratioout = npy.concatenate((ratio[:,512:622]   ,ratio[:,110:512]),axis=1)

lon1d = lonout[0,:].squeeze()
lat1d = latout[:,0].squeeze()

spval = 0.
ratioout[npy.where(ratioout == spval)] = 1.

## create a version with 1 below 60S and north of 60N
ratiotmp = ratioout.copy()

# au sud de 50S
#ratiotmp[198:255,:] = 1.

# au sud de 60S
ratiotmp[213:,:] = 1.
# au nord de 60N
ratiotmp[:43,:] = 1.

## deal with accepatable min and max
zmax = 1.15
zmin = 0.92

ratiotmp[npy.where(ratiotmp >= zmax)] = zmax
ratiotmp[npy.where(ratiotmp <= zmin)] = zmin

## merge the 2 version
ratiofin = ratiotmp.copy()

ratiofin[212,:] = 0.1 * ratioout[212,:] + 0.9 * 1.
ratiofin[211,:] = 0.2 * ratioout[211,:] + 0.8 * 1.
ratiofin[210,:] = 0.3 * ratioout[210,:] + 0.7 * 1.
ratiofin[209,:] = 0.4 * ratioout[209,:] + 0.6 * 1.
ratiofin[208,:] = 0.5 * ratioout[208,:] + 0.5 * 1.
ratiofin[207,:] = 0.6 * ratioout[207,:] + 0.4 * 1.
ratiofin[206,:] = 0.7 * ratioout[206,:] + 0.3 * 1.
ratiofin[205,:] = 0.8 * ratioout[205,:] + 0.2 * 1.
ratiofin[204,:] = 0.9 * ratioout[204,:] + 0.1 * 1.

ratiofin[43,:] = 0.1 * ratioout[43,:] + 0.9 * 1.
ratiofin[44,:] = 0.2 * ratioout[44,:] + 0.8 * 1.
ratiofin[45,:] = 0.3 * ratioout[45,:] + 0.7 * 1.
ratiofin[46,:] = 0.4 * ratioout[46,:] + 0.6 * 1.
ratiofin[47,:] = 0.5 * ratioout[47,:] + 0.5 * 1.
ratiofin[48,:] = 0.6 * ratioout[48,:] + 0.4 * 1.
ratiofin[49,:] = 0.7 * ratioout[49,:] + 0.3 * 1.
ratiofin[50,:] = 0.8 * ratioout[50,:] + 0.2 * 1.
ratiofin[51,:] = 0.9 * ratioout[51,:] + 0.1 * 1.

PyRaf.write_2d_reg_file('./ratio_windmodule_u.nc',lon1d,lat1d,0.,ratiofin,'u10')
PyRaf.write_2d_reg_file('./ratio_windmodule_v.nc',lon1d,lat1d,0.,ratiofin,'v10')
