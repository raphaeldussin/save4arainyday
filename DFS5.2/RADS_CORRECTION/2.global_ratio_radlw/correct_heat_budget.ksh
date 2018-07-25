#!/bin/ksh

## To correct the heat budget of DFS5.2.0 (-3.15 W/m2), we apply a ratio
## of 1.0088 on the longwave radiation

DIRIN=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/working_versions/DFS5.2.0/radlw
DIROUT=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/radlw

fyear=1980
lyear=2012

for year in  $( seq $fyear $lyear ) ; do

    rsync -av $DIRIN/drowned_radlw_DFS5.2.0_y${year}.nc .
    ncap2 -s " radlw = radlw * float(1.0088) " ./drowned_radlw_DFS5.2.0_y${year}.nc -o ./drowned_radlw_DFS5.2_y${year}.nc
    sudo -u viking rsync -av ./drowned_radlw_DFS5.2_y${year}.nc $DIROUT/drowned_radlw_DFS5.2_y${year}.nc
    rm drowned_radlw_DFS5.2*

done
