#!/bin/ksh

INDIR=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2
SRCDIR=/home/users/dussin/TOOLS/FORCING_TOOLS/CLIMATO

for var in precip snow radlw radsw ; do

    ln -s $INDIR/$var/drowned_${var}_DFS5.2_y* .
    #ncea $INDIR/$var/drowned_${var}_DFS5.2_y????.nc -o ./drowned_${var}_DFS5.2_y1958-1978.nc
    $SRCDIR/mkclimato_daily drowned_${var}_DFS5.2_y????.nc 
    rm drowned_${var}_DFS5.2_y*nc
    mv climato.nc ./drowned_${var}_DFS5.2_y1958-1978.nc
    mv climato-f.nc ./drowned_${var}_DFS5.2_y1958-1978_filtered.nc

done
