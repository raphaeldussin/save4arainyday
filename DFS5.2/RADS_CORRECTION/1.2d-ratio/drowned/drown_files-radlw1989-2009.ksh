#!/bin/ksh

SRCDIR=~molines/sosie_2012/bin/

for year in $( seq 1989 2009 ) ; do

#    $SRCDIR/mask_drown_field.x -D -i ../nodrown/radsw_DFS5.1_y${year}.nc -v radsw -x lon0 -y lat0 \
#                               -m ../nodrown/lsm_erainterim.nc -o drowned_radsw_DFS5.1_y${year}.nc

    $SRCDIR/mask_drown_field.x -D -i ../nodrown/radlw_DFS5.1_y${year}.nc -v radlw -x lon0 -y lat0 \
                               -m ../nodrown/lsm_erainterim.nc -o drowned_radlw_DFS5.1_y${year}.nc
done


