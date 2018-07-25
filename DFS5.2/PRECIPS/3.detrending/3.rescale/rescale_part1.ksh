#!/bin/ksh

DATASET=ERAinterim_corr_astorto_lowlats
## create full period mean
cd ../original_blend

ncra drowned_precip_ERAinterim_corr_astorto_lowlats_y????.nc -o ../3.rescale/precip_${DATASET}_1y_1979-2012.nc

cd ../3.rescale

fyear=1979
lyear=1991

SRC=/home/users/dussin/TOOLS/DFS5-tools
DIRIN=../2.detrend


moyperiod=mean_precip_${DATASET}_detrended_y${fyear}-${lyear}.nc
moyglobal=precip_${DATASET}_1y_1979-2012.nc

for year in $( seq $fyear $lyear ) ; do

    filein=precip_${DATASET}_detrended_y${year}.nc
    echo $filein
    $SRC/rescale_precip $DIRIN/$filein $moyglobal $DIRIN/$moyperiod 
    mv precip_rescaled.nc precip_${DATASET}_dtrd_y${year}.nc

done
