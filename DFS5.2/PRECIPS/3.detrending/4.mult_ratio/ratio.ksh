#!/bin/ksh

for year in $( seq 1979 2010 ) ; do

    filein=../3.rescale/precip_ERAinterim_corr_astorto_lowlats_dtrd_y$year.nc
    fileout=./precip_DFS5.1.2_y$year.nc
    ncap2 -s "precip = precip * 1.07" $filein $fileout

done
