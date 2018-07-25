import PyRaf
import PyRaf_plots
import numpy as npy
import matplotlib.pylab as plt

dir_ERAint='/fsnet/data/meom/workdir/dussin/TMPDIR_TQ2_CORR_ANTARCTIC/monthly/'
dir_POLES='/home/users/dussin/DFS5/DFS5.2/TQ2/t2q2_arctic_correction/0-data/'
dir_MASK='/fsnet/data/meom/workdir/dussin/MESHMASKS/'

f_ERAint = 'drowned_t2_ERAinterim_corr_antarctic_y1979-1998_monthly.nc'
f_POLES  = 't2_POLES-ERAi_1979-1998.nc'
f_mask = 'lsm_erainterim.nc'

t2_ERAint = PyRaf.readfull(dir_ERAint + f_ERAint,'t2')
t2_POLES  = PyRaf.readfull(dir_POLES + f_POLES,'t2')

lon = PyRaf.readfull(dir_ERAint + f_ERAint,'lon0')
lat = PyRaf.readfull(dir_ERAint + f_ERAint,'lat0')
lsm = PyRaf.readfull(dir_MASK + f_mask,'lsm')

print t2_ERAint.shape
print t2_POLES.shape

for kt in npy.arange(12):
	#plt.figure()
        contours = npy.arange(-10,10,0.1)
	limits = [-10,10,-10,10]
	tab = (t2_ERAint[kt,:,:] - t2_POLES[kt,:,:]) * lsm
	tab = PyRaf.mask_value(tab,0.)
	tab = PyRaf.mask_value_up(tab,50.)
	fname='diff_t2_ERAint-POLES_m' + str(kt+1) + '.png'
	name='t2 ERAinterim - POLES month ' + str(kt + 1)
	PyRaf_plots.reg_global_diffplot(lon,lat,tab,contours,limits,name=name,filename=fname)
	#plt.contourf(t2_ERAint[kt,:,:] - t2_POLES[kt,:,:],contours) ; plt.colorbar()
	#plt.savefig('diff_t2_ERAint-POLES_m' + str(kt) + '.png')
