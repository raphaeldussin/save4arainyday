#!/bin/ksh

HERE=$(pwd)

list=''
for year in $( seq 1984 2006 ) ; do

    list="$list drowned_radsw_ERAinterim_y${year}.nc"

done

cd /data1/dussin/ERAinterim_1979-2010/drowned/all

ncra -F $list -o $HERE/drowned_radsw_ERAinterim_1y_1984-2006.nc
