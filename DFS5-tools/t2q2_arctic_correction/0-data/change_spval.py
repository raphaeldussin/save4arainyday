
import PyRaf
import numpy as np

monnc = '../TRUCS_LB/ifrac_monthly_SSMI-ERA40_1979-1998.nc'

lon   = PyRaf.readfull( monnc , 'lon' )
lat   = PyRaf.readfull( monnc , 'lat' )
time  = PyRaf.readfull( monnc , 'time' )
ifrac = PyRaf.readfull( monnc , 'ifrac' )

ifracout = ifrac.copy()

ifracout[np.where(ifracout <= -9998.)] = 0.

PyRaf.write_2dpt_file('./temp.nc',lon,lat,time,ifracout,'ifrac')
