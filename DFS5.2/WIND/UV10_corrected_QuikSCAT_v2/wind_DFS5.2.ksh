#!/bin/ksh

INDIR=/fsnet/data/meom/workdir/dussin/TMPDIR_DFS5.2/drowned
OUTDIR=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2

for year in $(seq 1979 2012) ; do

    for var in u10 v10 ; do

    filein=drowned_${var}_DFS5.2_y${year}.nc
    fileout=drowned_${var}_DFS5.2_y${year}.nc
    sudo -u viking nice -18 ionice -c 3 rsync -av $INDIR/$filein $OUTDIR/$var/$fileout 

    done

done
