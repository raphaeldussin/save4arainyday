#!/bin/ksh

dirin=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/ERA_INTERIM/drowned/
dirout=/fsnet/data/meom/workdir/dussin/TMPDIR_PRECIPS_DFS5.2/ERAinterim_ORCA025

HERE=$( pwd )

ybeg=2011
yend=2012

LYEAR=$( seq $ybeg $yend )

for year in $LYEAR ; do

    ln -s $dirin/drowned_precip_ERAinterim_y$year.nc

    cat namelist.precip.skel | sed -e "s;<DIRIN>;$dirin;g" -e "s;<DIROUT>;$dirout;g" \
                                   -e "s;<YEAR>;$year;g" > namelist.precip.$year

    $HOME/BIN/SOSIE/sosie.x -f namelist.precip.$year

done
