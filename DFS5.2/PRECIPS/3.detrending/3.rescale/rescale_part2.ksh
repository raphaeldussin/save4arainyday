#!/bin/ksh

fyear=1992
lyear=2003

SRC=/home/users/dussin/TOOLS/DFS5-tools
DIRIN=../2.detrend

DATASET=ERAinterim_corr_astorto_lowlats

moyperiod=mean_precip_${DATASET}_detrended_y${fyear}-${lyear}.nc
moyglobal=precip_${DATASET}_1y_1979-2012.nc

for year in $( seq $fyear $lyear ) ; do

    filein=precip_${DATASET}_detrended_y${year}.nc
    echo $filein
    $SRC/rescale_precip $DIRIN/$filein $moyglobal $DIRIN/$moyperiod 
    mv precip_rescaled.nc precip_${DATASET}_dtrd_y${year}.nc

done
