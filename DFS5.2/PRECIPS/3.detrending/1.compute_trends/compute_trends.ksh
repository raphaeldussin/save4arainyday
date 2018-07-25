#!/bin/ksh

FARC=/home/users/dussin/TOOLS/FARC/utils

# compute annual means
for year in $(seq 1979 2012) ; do

    ncra ../original_blend/drowned_precip_ERAinterim_corr_astorto_lowlats_y$year.nc -o ./mean_y$year.nc

done


#####################################################################################################

list=''
for year in $(seq 1979 1991) ; do
    list="$list mean_y$year.nc"
done

    $FARC/trend -var precip -dataset ERAint_STORTO_LOWLAT -fyear 1979 -lyear 1991 -diroutput ./ \
            -time time $list 

    mv trend.nc trend_ERAint_Storto_lowlat_1979-1991.nc



#####################################################################################################

list=''
for year in $(seq 1992 2003) ; do
    list="$list mean_y$year.nc"
done

    $FARC/trend -var precip -dataset ERAint_STORTO_LOWLAT -fyear 1992 -lyear 2004 -diroutput ./ \
            -time time $list 

    mv trend.nc trend_ERAint_Storto_lowlat_1992-2003.nc


#####################################################################################################

list=''
for year in $(seq 2004 2012) ; do
    list="$list mean_y$year.nc"
done

    $FARC/trend -var precip -dataset ERAint_STORTO_LOWLAT -fyear 2005 -lyear 2010 -diroutput ./ \
            -time time $list 

    mv trend.nc trend_ERAint_Storto_lowlat_2004-2012.nc


rm mean*nc

