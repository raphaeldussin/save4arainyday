#!/bin/ksh

T2dir=/fsnet/data/meom/workdir/dussin/TMPDIR_TQ2_CORR_ARCTIC/t2_2poles_correct/drowned
T2dir2=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/t2/

Q2dir=/fsnet/data/meom/workdir/dussin/TMPDIR_TQ2_CORR_ARCTIC/q2_2poles_correct/drowned
Q2dir2=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/q2/

for year in $( seq 1979 2012 ) ; do

    sum1=$( md5sum $T2dir/drowned_t2_DFS5.2_y$year.nc   | awk '{ print $1 }' )
    sum2=$( md5sum $T2dir2/drowned_t2_DFS5.2_y$year.nc  | awk '{ print $1 }' )

    if [ $sum1 = $sum2 ] ; then
       echo drowned_t2_DFS5.2_y$year.nc : copy is OK
       echo $sum1 $sum2
    fi

    sum3=$( md5sum $Q2dir/drowned_q2_DFS5.2_y$year.nc   | awk '{ print $1 }' )
    sum4=$( md5sum $Q2dir2/drowned_q2_DFS5.2_y$year.nc  | awk '{ print $1 }' )

    if [ $sum3 = $sum4 ] ; then
       echo drowned_q2_DFS5.2_y$year.nc : copy is OK
       echo $sum3 $sum4
    fi
   
done 
