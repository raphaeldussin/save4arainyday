#!/bin/ksh

fyear=1979
lyear=1979

DIRIN=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/ERA_INTERIM/1979_2010/drowned/all

for year in $( seq $fyear $lyear ) ; do

    filein=drowned_precip_ERAinterim_y${year}.nc
    filesnow=drowned_snow_ERAinterim_y${year}.nc
    fileout=precip_liquid_ERAinterim_y${year}.nc

    cp $DIRIN/$filein .
    cp $DIRIN/$filesnow .
    ncrename -v snow,precip $filesnow

    ncap -O -s " precip = precip * 86400" $filein -o $filein
    ncap -O -s " precip = precip * 86400" $filesnow -o $filesnow

    ncdiff $filein $filesnow -o $fileout
    #rm $filein $filesnow

done

