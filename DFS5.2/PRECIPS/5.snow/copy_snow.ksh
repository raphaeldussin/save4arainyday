#!/bin/ksh

## Snow is the same in ERAinterim and DFS5.2

INDIR=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/ERA_INTERIM/drowned
OUTDIR=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/snow

fyear=1980
lyear=2012

for year in $( seq $fyear $lyear ) ; do

    sudo -u viking rsync -av $INDIR/drowned_snow_ERAinterim_y${year}.nc $OUTDIR/drowned_snow_DFS5.2_y${year}.nc

done
