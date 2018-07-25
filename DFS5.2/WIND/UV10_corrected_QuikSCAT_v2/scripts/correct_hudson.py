import PyRaf

lon = PyRaf.readfull('./ratio_windmodule_u.nc','lon0')
lat = PyRaf.readfull('./ratio_windmodule_u.nc','lat0')
u10 = PyRaf.readfull('./ratio_windmodule_u.nc','u10')

u10out = u10.copy()

# mask hudson bay
u10out[40:55,370:400] = 1.

# remove 1 for new computation
u10out = u10out -1.

PyRaf.write_2d_reg_file('./new_ratio_windmodule_u.nc', lon, lat, 0., u10out, 'u10')

lon = PyRaf.readfull('./ratio_windmodule_v.nc','lon0')
lat = PyRaf.readfull('./ratio_windmodule_v.nc','lat0')
v10 = PyRaf.readfull('./ratio_windmodule_v.nc','v10')

v10out = v10.copy()

# mask hudson bay
v10out[40:55,370:400] = 1.

# remove 1 for new computation
v10out = v10out -1.

PyRaf.write_2d_reg_file('./new_ratio_windmodule_v.nc', lon, lat, 0., v10out, 'v10')

