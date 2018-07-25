import PyRaf
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
import numpy as npy
import scipy.ndimage as sn
from matplotlib import mpl
import os

swdfs = PyRaf.readfull('./radsw_DFS4.3-512x256_1y_1984-2006.nc','radsw')
swera = PyRaf.readfull('./drowned_radsw_ERAinterim_1y_1984-2006.nc','radsw')

lon0 = PyRaf.readfull('./lsm_erainterim.nc','lon0')
lat0 = PyRaf.readfull('./lsm_erainterim.nc','lat0')
lsm  = PyRaf.readfull('./lsm_erainterim.nc','lsm')

diff_sw = ( swera - swdfs ) * lsm

# we keep only differences greater than 10 W/m2
diff_sw = PyRaf.mask_value(diff_sw,0.)
diff_sw = PyRaf.mask_value_down(diff_sw,10.)

alpha_sw = ( swdfs / swera ) * lsm
alpha_sw[ npy.where( alpha_sw > 1 ) ] = 1.

# we apply the diff mask on the ratio
ratio = npy.ma.array(data=alpha_sw,mask=diff_sw.mask)
ratio[ npy.where(ratio.mask) ] = 1.

ratio[:112,:62] = 1.
ratio[:80,:90] = 1.
ratio[:96,:81] = 1.
ratio[:36,:] = 1.
ratio[:55,350:395] = 1.
ratio[210:,60:440] = 1.

#ratio = ratio * lsm

# This is the first part
PyRaf.write_2d_reg_file('./undrowned_ratio_radsw.nc', lon0, lat0, 0., ratio, 'radsw')


os.system('~molines/sosie_2012/bin/mask_drown_field.x -D -i undrowned_ratio_radsw.nc -v radsw -m lsm_erainterim.nc')
os.system('mv fout.nc drowned_ratio_radsw.nc')

ratio = PyRaf.readfull('./drowned_ratio_radsw.nc','radsw')

# ratio becomes 1 on land
#ratio[ npy.where( ratio == 0 ) ] = 1.
# we extend to prevent discontinuity at lon = 0
ratio_extended    = npy.concatenate((ratio,ratio),axis=1)
# smooth with a gaussian filter
smoothed_extended = sn.gaussian_filter(ratio_extended, sigma=3)
smoothed = smoothed_extended[:,:512]

## small adjustements
smoothed[204:223,330:350] = 1.
smoothed[180:200,234:262] = 1.
smoothed[218:225,82:93] = 1.
smoothed[34:50,375:395] = 1.
smoothed[94:106,50:61] = 1.
smoothed[11:34,0:32] = 1.
smoothed[:85,:30] = 1.
smoothed[:85,504:] = 1.
smoothed[113:123,427:440] = 1.

# in order to remove residus from the filtering
smoothed[npy.where( smoothed >= 0.99) ] = 1.

smoothed = smoothed * lsm

PyRaf.write_2d_reg_file('./smoothed_ratio_radsw.nc', lon0, lat0, 0., smoothed, 'radsw')

## figures
lon, lat = npy.meshgrid(lon0,lat0)

lonplt = npy.concatenate((lon-360,lon))
latplt = npy.concatenate((lat,lat))

tab = smoothed
#tabplt = npy.concatenate((tab,tab),axis=1)
tab = PyRaf.mask_value(tab,1.)

norm = mpl.colors.Normalize(vmin=0.75, vmax=1.)
contours = npy.arange(0.75,1.00,0.01)

plt.figure(figsize=[12.,10.])
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='grey',lake_color='white')
# draw parallels and meridians.
m.drawparallels(npy.arange(-90.,91.,30.))
m.drawmeridians(npy.arange(0.,360.,30.))
m.drawmapboundary(fill_color='white') 

m.contourf(lon,lat,tab.squeeze(),contours,norm=norm) ; plt.colorbar(shrink=0.45,ticks=[0.75,0.8,0.85,0.9,0.95,1.0,1.05])

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

plt.title('Multiplicative ratio applied on ERAinterim shortwave radiation\n (values equal to unity are masked)',fontsize=16)
plt.savefig('./smoothed_ratio_sw.png')
plt.close()


#plt.figure()
#plt.contourf(npy.arange(0,512),-npy.arange(0,256),tab); plt.colorbar()
##plt.contourf(npy.arange(0,1024),-npy.arange(0,256),tabplt); plt.colorbar()
#plt.grid()
##plt.xticks((npy.arange(0,532,20),npy.arange(0,532,20)))
##plt.yticks((npy.arange(0,276,20),npy.arange(0,276,20)))
#plt.show()
