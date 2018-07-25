#!/bin/ksh

epd=/home/users/dussin/BIN/EPD/epd-7.0-2-rh5-x86_64/bin
cdftools=/home/users/dussin/TOOLS/CDFTOOLS_3.0/bin

### This is the directory where wind module files (from FOTO) are stored
#
DIR1=/media/sdc1/dussin/RUNS_FOTO/PREPARATION_OF_DFS5/ERAinterim_artic_corr/wind_module

## compute annual means for wind module
for year in $( seq 2000 2008 ) ; do

   echo working on year $year
   ncra -F $DIR1/U10_6h_$year.nc -o mean_U10_ERAi_$year.nc

done

## compute the interannual mean
ncra -F mean_U10_ERAi_????.nc -o mean_U10_ERAi_2000-2008.nc

## compute a raw ratio
ncbo --op_typ=/ U10_QuikScat-512x256_2000-2008.nc mean_U10_ERAi_2000-2008.nc -o ratio_QuikSCAT_on_ERAi_2000-2008.nc

## this ratio have to be rewritten in nemo form to use cdftools
$epd/python rewrite_nc.py

## then apply cdf2matlab to extend space
$cdftools/cdf2matlab nemo_ratio_QuikSCAT_on_ERAi_2000-2008.nc ratio 1

## then rewrite
$epd/python rewrite_output.py

## then apply cdfsmooth
$cdftools/cdfsmooth extended_ouput.nc 7 B 2.

## rewrite again but in forcing-style file
$epd/python rewrite_smoothed.py

ncrename -v lon,lon0 -v lat,lat0 -d lon,lon0 -d lat,lat0 ratio_windmodule_u.nc
ncrename -v lon,lon0 -v lat,lat0 -d lon,lon0 -d lat,lat0 ratio_windmodule_v.nc


