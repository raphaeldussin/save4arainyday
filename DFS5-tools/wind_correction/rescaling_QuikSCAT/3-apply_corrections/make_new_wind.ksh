#!/bin/ksh

DIRIN=/data1/dussin/ERAinterim_1979-2010/drowned/all

coefu=correction_u10.nc
coefv=correction_v10.nc

for year in $(seq 1980 2010) ; do

    fileu=drowned_u10_ERAinterim_y$year.nc
    filev=drowned_v10_ERAinterim_y$year.nc

    leapyear=0
    yearref=1976
    leaptest=$( echo $year $yearref | awk '{print ($1 - $2)%4}' )
    if (( $leaptest == 0 )) ; then
        echo $year is a leap year
        leapyear=1
    fi

 
    nframes=$(( 2920 + $leapyear * 8 ))
    echo $year has $nframes

    for kt in $( seq 1 $nframes ) ; do
 
        tttt=$( printf "%04d" $kt )

        filetmpu=u10_INTERIM-512x256_y${year}_${tttt}.nc
        filetmpv=v10_INTERIM-512x256_y${year}_${tttt}.nc
        fileoutu=u10_corr_INTERIM-512x256_y${year}_${tttt}.nc
        fileoutv=v10_corr_INTERIM-512x256_y${year}_${tttt}.nc

        ncks -F -d time,$kt,$kt $DIRIN/$fileu -o ./$filetmpu
        ncks -F -d time,$kt,$kt $DIRIN/$filev -o ./$filetmpv

        ncbo --op_typ=+ ./$filetmpu ./$coefu -o ./$fileoutu
        ncbo --op_typ=+ ./$filetmpv ./$coefv -o ./$fileoutv

    done

    ncrcat u10_corr_INTERIM-512x256_y${year}_????.nc -o u10_ERAinterim_corr_QuikSCAT_v2_y${year}.nc
    ncrcat v10_corr_INTERIM-512x256_y${year}_????.nc -o v10_ERAinterim_corr_QuikSCAT_v2_y${year}.nc

    rm -f u10_INTERIM-512x256_y${year}_????.nc
    rm -f v10_INTERIM-512x256_y${year}_????.nc
    rm -f u10_corr_INTERIM-512x256_y${year}_????.nc
    rm -f v10_corr_INTERIM-512x256_y${year}_????.nc

done
