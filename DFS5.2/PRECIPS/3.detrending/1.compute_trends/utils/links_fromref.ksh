#!/bin/ksh

# In some case some fields can have the same values in two
# experiments so that linking can help save space on disk

EXPREF=ERAinterim
EXPNEW=ERAinterim_corr_QuikSCAT

# Where results from FARC are stored
FARC_OUT=/media/sdc1/dussin/FARC_OUT

# List of variables to proceed
LIST="Qrad Qlw Qsw precip"


###################################################################################################
###################################################################################################
## 1. In climato data

for var in $LIST ; do

    listfile=$( ls $FARC_OUT/$EXPREF/climatos_data | grep $var )

    for file in $listfile ; do

        fileout=$( echo $file | sed -e "s/$EXPREF/$EXPNEW/g" )
        ln -s $FARC_OUT/$EXPREF/climatos_data/$file $FARC_OUT/$EXPNEW/climatos_data/$fileout

    done

done

## 2. In timeseries data

for var in $LIST ; do

    listfile=$( ls $FARC_OUT/$EXPREF/timeseries_data | grep $var )

    for file in $listfile ; do

        fileout=$( echo $file | sed -e "s/$EXPREF/$EXPNEW/g" )
        cp $FARC_OUT/$EXPREF/timeseries_data/$file $FARC_OUT/$EXPNEW/timeseries_data/$fileout
        ncrename -v ${var}_${EXPREF},${var}_${EXPNEW} $FARC_OUT/$EXPNEW/timeseries_data/$fileout

    done

done

## 3. In trends data

for var in $LIST ; do

    listfile=$( ls $FARC_OUT/$EXPREF/trends_data | grep $var )

    for file in $listfile ; do

        fileout=$( echo $file | sed -e "s/$EXPREF/$EXPNEW/g" )
        ln -s $FARC_OUT/$EXPREF/trends_data/$file $FARC_OUT/$EXPNEW/trends_data/$fileout

    done

done

