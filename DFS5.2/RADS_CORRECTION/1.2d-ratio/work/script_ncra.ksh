#!/bin/ksh

HERE=$(pwd)

list=''
for year in $( seq 1984 2006 ) ; do

    list="$list drowned_radlw_ERAinterim_y${year}.nc"

done

cd /data1/dussin/ERAinterim

ncra -F $list -o $HERE/drowned_radlw_ERAinterim_1y_1984-2006.nc
