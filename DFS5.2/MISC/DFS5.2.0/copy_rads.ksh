#!/bin/ksh

FROMDIR=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.1.1
TODIR=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/working_versions/DFS5.2.0

for year in $( seq 1979 2010 ) ; do

    rsync -av $FROMDIR/radsw/drowned_radsw_DFS5.1.1_y${year}.nc $TODIR/radsw/drowned_radsw_DFS5.2.0_y${year}.nc
    rsync -av $FROMDIR/radlw/drowned_radlw_DFS5.1.1_y${year}.nc $TODIR/radlw/drowned_radlw_DFS5.2.0_y${year}.nc

done
