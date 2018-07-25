import PyRaf

lon = PyRaf.readfull('./new_ratio_windmodule_u.nc','lon0')
lat = PyRaf.readfull('./new_ratio_windmodule_u.nc','lat0')
alp = PyRaf.readfull('./new_ratio_windmodule_u.nc','u10')

u10 = PyRaf.readfull('./mean_u10.nc','u10')

out = alp * u10

PyRaf.write_2d_reg_file('./correction_u10.nc', lon, lat, 0., out, 'u10')


lon = PyRaf.readfull('./new_ratio_windmodule_v.nc','lon0')
lat = PyRaf.readfull('./new_ratio_windmodule_v.nc','lat0')
alp = PyRaf.readfull('./new_ratio_windmodule_v.nc','v10')

v10 = PyRaf.readfull('./mean_v10.nc','v10')

out2 = alp * v10

PyRaf.write_2d_reg_file('./correction_v10.nc', lon, lat, 0., out2, 'v10')
