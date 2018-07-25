#!/bin/ksh

#listdataset='DFS5.1.1 ERAinterim'
listdataset='DFS4.4 '
#listvars='precip t2 q2 u10 v10 radsw radlw'
listvars='snow'
#listvars='Qnet Qtrb Qlat Qsen Qrad Qlw Qsw U10 Taux Tauy EmP evap '
FARCOUT=/fsnet/data/meom/workdir/dussin/FARC_OUT
PYT=/home/users/dussin/python/epd-7.0-2-rh5-x86_64/bin/python

### Climatos

echo '#####################################################'
echo '### ANNUAL CLIMATO '
echo '#####################################################'

for var in $listvars ; do

    for dataset in $listdataset ; do

        # find the min/max for annual climato
	ffound=$( ls $FARCOUT/$dataset/climatos_data | grep $var | grep $dataset | grep 1y | grep nc | awk '{print $1}' )
        filetmp=$FARCOUT/$dataset/climatos_data/$ffound
        tmp=$( $PYT minmax.py $filetmp $var )
        echo $var $dataset $tmp

    done

done

echo '#####################################################'
echo '### MONTHLY CLIMATO '
echo '#####################################################'

for var in $listvars ; do

    for dataset in $listdataset ; do

        # find the min/max for annual climato
	ffound=$( ls $FARCOUT/$dataset/climatos_data | grep $var | grep $dataset | grep 1m | grep nc | awk '{print $1}' )
        filetmp=$FARCOUT/$dataset/climatos_data/$ffound
        tmp=$( $PYT minmax.py $filetmp $var )
        echo $var $dataset $tmp

    done

done

#### timeseries

echo '#####################################################'
echo '### TIMESERIES GLOBAL'
echo '#####################################################'
for var in $listvars ; do

    for dataset in $listdataset ; do

        # find the min/max for annual climato
	ffound=$( ls $FARCOUT/$dataset/timeseries_data | grep $var | grep $dataset | grep global | grep nc | \
                  grep full_ | awk '{print $1}' )
        filetmp=$FARCOUT/$dataset/timeseries_data/$ffound
        tmp=$( $PYT minmax.py $filetmp ${var}_${dataset} )
        echo $var $dataset $tmp

    done

done

for area in Polar_South Subpolar_South Subtropical_South Tropical_South Equatorial_Band Tropical_North \
                Subtropical_North Subpolar_North Polar_North global ; do

    echo '#####################################################'
    echo "### TIMESERIES $area"
    echo '#####################################################'


    for var in $listvars ; do

    for dataset in $listdataset ; do

        # find the min/max for annual climato
        ffound=$( ls $FARCOUT/$dataset/timeseries_data | grep $var | grep $dataset | grep $area | grep nc | \
                  grep full_ | awk '{print $1}' )
        filetmp=$FARCOUT/$dataset/timeseries_data/$ffound
        tmp=$( $PYT minmax.py $filetmp ${var}_${dataset} )
        echo $var $dataset $tmp

    done

    done

done








