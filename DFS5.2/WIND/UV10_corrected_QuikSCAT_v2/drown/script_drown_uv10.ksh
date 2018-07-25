#!/bin/ksh

MYTMPDIR=/fsnet/data/meom/workdir/dussin/TMPDIR_DFS5.2
HERE=$( pwd )

ybeg=1979
yend=1979

LYEAR=$( seq $ybeg $yend )

cd $MYTMPDIR
if [ ! -d drowned ] ; then mkdir drowned ; fi
cd drowned

for year in $LYEAR ; do

   for var in u10 v10 ; do

     FILEIN=${var}_DFS5.2_y${year}.nc

     ln -s $MYTMPDIR/$FILEIN .

     cat ${HERE}/namelist.$var.skel | sed -e "s/<FILE>/$FILEIN/g" \
                                -e "s/<YEAR>/$year/g" > namelist.$var.$year

    ~/BIN/SOSIE/sosie.x -f namelist.$var.$year

    mv ${var}_ERAinterim-ERAinterim_${year}.nc drowned_${var}_DFS5.2_y${year}.nc

    done

done
