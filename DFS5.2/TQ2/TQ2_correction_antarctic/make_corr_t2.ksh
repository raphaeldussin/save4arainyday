#!/bin/ksh

MYTMPDIR=/fsnet/data/meom/workdir/dussin/TMPDIR_TQ2_CORR_ANTARCTIC
MYINPUTS=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/ERA_INTERIM/drowned

HERE=$( pwd )

cd $MYTMPDIR

for year in $(seq 1979 2012 ) ; do

    $HERE/correct_t2 t2 $MYINPUTS/drowned_t2_ERAinterim_y${year}.nc
    mv output.nc drowned_t2_ERAinterim_corr_antarctic_y${year}.nc

done
