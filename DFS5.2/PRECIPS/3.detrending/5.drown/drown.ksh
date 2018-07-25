#!/bin/ksh

SOSIE=/home/users/dussin/BIN/SOSIE

for year in $(seq 1979 2010) ; do

    cat namelist.precip | sed -e "s/<YEAR>/$year/g" > namelist.precip.$year
    $SOSIE/sosie.x -f namelist.precip.$year
    mv precip_DFS5.1.2-drown_y${year}.nc drowned_precip_DFS5.1.2_y${year}.nc

done
