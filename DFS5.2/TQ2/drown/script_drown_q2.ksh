#!/bin/ksh

MYTMPDIR=/fsnet/data/meom/workdir/dussin/TMPDIR_TQ2_CORR_ARCTIC/q2_2poles_correct
flonlat=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/u10/drowned_u10_DFS5.2_y1979.nc

HERE=$( pwd )

ybeg=2009
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

   for var in q2 ; do

     FILEIN=${var}_DFS5.2_${year}.nc

     ln -s $MYTMPDIR/$FILEIN .

    ./mask_drown_field.x -i q2_DFS5.2_${year}.nc -v q2 -D -m lsm_erainterim.nc -o drowned_q2_DFS5.2_y${year}.nc

    ncrename -h -d x,lon0 -d y,lat0 drowned_q2_DFS5.2_y${year}.nc
    ncks -h -F -A -v lon0,lat0 $flonlat -o drowned_q2_DFS5.2_y${year}.nc
    mv drowned_q2_DFS5.2_y${year}.nc tmp.nc
    ncks -h -F -v lon0,lat0,time,$var tmp.nc -o drowned_q2_DFS5.2_y${year}.nc
    rm tmp.nc

    done

done
