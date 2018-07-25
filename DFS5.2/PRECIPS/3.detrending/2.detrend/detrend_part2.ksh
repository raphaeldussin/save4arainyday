#!/bin/ksh

fyear=1992
lyear=2003

SRC=/home/users/dussin/TOOLS/DFS5-tools
DIRIN=/fsnet/data/meom/workdir/dussin/Precip_detrending_v3/original_blend
DATASET=ERAinterim_corr_astorto_lowlats # blending of Storto 30N/30S and original ERAinterim

trendin=trend_ERAint_Storto_lowlat_${fyear}-${lyear}.nc
ln -s ../1.compute_trends/$trendin .

LIST=''

for year in $( seq $fyear $lyear ) ; do

    filein=drowned_precip_${DATASET}_y${year}.nc
    ln -s $DIRIN/$filein .
    $SRC/detrend_precip_v2 $filein $trendin $year $fyear
    mv precip_detrended.nc precip_${DATASET}_detrended_y${year}.nc
    rm $filein
    LIST="$LIST precip_${DATASET}_detrended_y${year}.nc"

done

ncra $LIST -o mean_precip_${DATASET}_detrended_y$fyear-$lyear.nc

