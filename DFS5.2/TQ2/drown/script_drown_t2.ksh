#!/bin/ksh

MYTMPDIR=/fsnet/data/meom/workdir/dussin/TMPDIR_TQ2_CORR_ARCTIC/t2_2poles_correct
flonlat=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/u10/drowned_u10_DFS5.2_y1979.nc

HERE=$( pwd )

ybeg=2002
yend=2012

LYEAR=$( seq $ybeg $yend )

cd $MYTMPDIR
if [ ! -d drowned ] ; then 

   mkdir drowned 
   cp $HERE/lsm_erainterim.nc ./drowned/. 
   cp $HERE/mask_drown_field.x ./drowned/.

fi

cd drowned

for year in $LYEAR ; do

   for var in t2 ; do

     FILEIN=${var}_DFS5.2_y${year}.nc

     ln -s $MYTMPDIR/$FILEIN .

    ./mask_drown_field.x -i t2_DFS5.2_y${year}.nc -v t2 -D -m lsm_erainterim.nc -o drowned_t2_DFS5.2_y${year}.nc

    ncrename -h -d x,lon0 -d y,lat0 drowned_t2_DFS5.2_y${year}.nc
    ncks -h -F -A -v lon0,lat0 $flonlat -o drowned_t2_DFS5.2_y${year}.nc
    mv drowned_t2_DFS5.2_y${year}.nc tmp.nc
    ncks -h -F -v lon0,lat0,time,$var tmp.nc -o drowned_t2_DFS5.2_y${year}.nc
    rm tmp.nc

    done

done
