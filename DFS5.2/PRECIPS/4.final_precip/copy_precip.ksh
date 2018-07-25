#!/bin/ksh

## Snow is the same in DFS5.2.0 and DFS5.2

INDIR=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/working_versions/DFS5.2.0/precip
OUTDIR=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/precip

fyear=1980
lyear=2012

for year in $( seq $fyear $lyear ) ; do

    sudo -u viking rsync -av $INDIR/drowned_precip_DFS5.2.0_y${year}.nc $OUTDIR/drowned_precip_DFS5.2_y${year}.nc

done
