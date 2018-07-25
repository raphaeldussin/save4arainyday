#!/bin/ksh

DIRIN=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/ERA_INTERIM/drowned
MYTMPDIR=/fsnet/data/meom/workdir/dussin/TMPDIR_DFS5.2

coefu=correction_u10.nc
coefv=correction_v10.nc

cp $coefu $MYTMPDIR
cp $coefv $MYTMPDIR

cd $MYTMPDIR

ncrename -v lon,lon0 -v lat,lat0 -d lon,lon0 -d lat,lat0 $coefu
ncrename -v lon,lon0 -v lat,lat0 -d lon,lon0 -d lat,lat0 $coefv

for year in $(seq 1979 1979) ; do

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

    ncrcat u10_corr_INTERIM-512x256_y${year}_????.nc -o u10_DFS5.2_y${year}.nc
    ncrcat v10_corr_INTERIM-512x256_y${year}_????.nc -o v10_DFS5.2_y${year}.nc

    rm -f u10_INTERIM-512x256_y${year}_????.nc
    rm -f v10_INTERIM-512x256_y${year}_????.nc
    rm -f u10_corr_INTERIM-512x256_y${year}_????.nc
    rm -f v10_corr_INTERIM-512x256_y${year}_????.nc

done
