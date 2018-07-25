import PyRaf
import PyRaf_plots
import numpy as npy
import matplotlib.pylab as plt

dir_offset='/home/users/dussin/DFS5/DFS5.2/TQ2/t2q2_arctic_correction/2-offset/'
dir_MASK='/fsnet/data/meom/workdir/dussin/MESHMASKS/'

f_offset = 'offset_poles-eraint_1979-1998.nc.sm4'
f_mask = 'lsm_erainterim.nc'

offset = PyRaf.readfull(dir_offset + f_offset,'offst')

lon = PyRaf.readfull(dir_MASK + f_mask,'lon')
lat = PyRaf.readfull(dir_MASK + f_mask,'lat')
lsm = PyRaf.readfull(dir_MASK + f_mask,'lsm')

for kt in npy.arange(12):
        contours = npy.arange(-10,10,0.1)
	limits = [-10,10,-10,10]
	tab = (- offset[kt,:,:] ) * lsm
	tab = PyRaf.mask_value(tab,0.)
	tab = PyRaf.mask_value_up(tab,50.)
	fname='offset_sm4_m' + str(kt+1) + '.png'
	name='-1 * Offset for month ' + str(kt + 1)
	PyRaf_plots.reg_global_diffplot(lon,lat,tab,contours,limits,name=name,filename=fname)
