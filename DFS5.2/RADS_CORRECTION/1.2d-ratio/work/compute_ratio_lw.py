import PyRaf
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
import numpy as npy
import scipy.ndimage as sn
from matplotlib import mpl
import os

lwdfs = PyRaf.readfull('./radlw_DFS4.3-512x256_1y_1984-2006.nc','radlw')
lwera = PyRaf.readfull('./drowned_radlw_ERAinterim_1y_1984-2006.nc','radlw')

lon0 = PyRaf.readfull('./lsm_erainterim.nc','lon0')
lat0 = PyRaf.readfull('./lsm_erainterim.nc','lat0')
lsm  = PyRaf.readfull('./lsm_erainterim.nc','lsm')

diff_lw = ( lwera - lwdfs ) * lsm

# we keep only differences greater than 10 W/m2
diff_lw = PyRaf.mask_value(diff_lw,0.)
diff_lw = PyRaf.mask_value_up(diff_lw,-2.5)

alpha_lw = ( lwdfs / lwera ) * lsm
#alpha_lw[ npy.where( alpha_lw > 1 ) ] = 1.

# we apply the diff mask on the ratio
ratio = npy.ma.array(data=alpha_lw,mask=diff_lw.mask)
ratio[ npy.where(ratio.mask) ] = 1.

ratio[:112,:62] = 1.
ratio[:80,:90] = 1.
ratio[:96,:81] = 1.
ratio[:45,:] = 1.
ratio[:55,375:405] = 1.
ratio[210:,:] = 1.
ratio[125:130,40:50] = 1.

#ratio = ratio * lsm

# This is the first part
PyRaf.write_2d_reg_file('./undrowned_ratio_radlw.nc', lon0, lat0, 0., ratio, 'radlw')


os.system('~molines/sosie_2012/bin/mask_drown_field.x -D -i undrowned_ratio_radlw.nc -v radlw -m lsm_erainterim.nc')
os.system('mv fout.nc drowned_ratio_radlw.nc')

ratio = PyRaf.readfull('./drowned_ratio_radlw.nc','radlw')

# ratio becomes 1 on land
#ratio[ npy.where( ratio == 0 ) ] = 1.
# we extend to prevent discontinuity at lon = 0
ratio_extended    = npy.concatenate((ratio,ratio),axis=1)
# smooth with a gaussian filter
smoothed_extended = sn.gaussian_filter(ratio_extended, sigma=3)
smoothed = smoothed_extended[:,:512]

## small adjustements
smoothed[86:93,73:78] = 1.

# in order to remove residus from the filtering
smoothed[npy.where( smoothed <= 1.01) ] = 1.

smoothed = smoothed * lsm

PyRaf.write_2d_reg_file('./smoothed_ratio_radlw.nc', lon0, lat0, 0., smoothed, 'radlw')

## figures
lon, lat = npy.meshgrid(lon0,lat0)

lonplt = npy.concatenate((lon-360,lon))
latplt = npy.concatenate((lat,lat))

tab = smoothed
#tabplt = npy.concatenate((tab,tab),axis=1)
tab = PyRaf.mask_value(tab,1.)

norm = mpl.colors.Normalize(vmin=1., vmax=1.08)
contours = npy.arange(1.00,1.08,0.01)

plt.figure(figsize=[12.,10.])
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='grey',lake_color='white')
# draw parallels and meridians.
m.drawparallels(npy.arange(-90.,91.,30.))
m.drawmeridians(npy.arange(0.,360.,30.))
m.drawmapboundary(fill_color='white') 

m.contourf(lon,lat,tab.squeeze(),contours,norm=norm) ; plt.colorbar(shrink=0.45,ticks=[1.0,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08])

# x axis
locs, labels = plt.xticks()
newlocs   = npy.array([0,30,60,90,120,150,180,210,240,270,300,330,360],'f')
newlabels = npy.array([0,30,60,90,120,150,180,210,240,270,300,330,0],'i')
plt.xticks(newlocs,newlabels)
plt.xlabel('Longitude',fontsize=14)
#
# y axis
locsY,labelsy = plt.yticks()
newlocsy   = npy.array([-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80],'f')
newlabelsy = npy.array([-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80],'i')
plt.yticks(newlocsy,newlabelsy)
plt.ylabel('Latitude',fontsize=14)

plt.title('Multiplicative ratio applied on ERAinterim longwave radiation\n (values equal to unity are masked)',fontsize=16)
plt.savefig('./smoothed_ratio_lw.png')
plt.close()


#plt.figure()
#plt.contourf(npy.arange(0,512),-npy.arange(0,256),tab); plt.colorbar()
##plt.contourf(npy.arange(0,1024),-npy.arange(0,256),tabplt); plt.colorbar()
#plt.grid()
##plt.xticks((npy.arange(0,532,20),npy.arange(0,532,20)))
##plt.yticks((npy.arange(0,276,20),npy.arange(0,276,20)))
#plt.show()
