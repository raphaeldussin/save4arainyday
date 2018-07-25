#!/bin/ksh

SRCDIR=/home/users/dussin/TOOLS/DFS5-tools/precip_rads_offline_corrections
ERADIR=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/ERA_INTERIM/drowned
# for 1989 to 2009
#ERADIR=/data1/dussin/ERAinterim # no those files are bugged

for year in $( seq 2011 2012 ) ; do

    ln -s $ERADIR/drowned_radlw_ERAinterim_y${year}.nc .
    $SRCDIR/correc_radlw_coef drowned_radlw_ERAinterim_y${year}.nc smoothed_ratio_radlw.nc
    mv radlw_corrected.nc radlw_DFS5.1_y${year}.nc
    rm drowned_radlw_ERAinterim_y${year}.nc

done
    
