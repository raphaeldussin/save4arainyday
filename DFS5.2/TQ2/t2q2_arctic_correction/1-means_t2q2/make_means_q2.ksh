#!/bin/ksh

DIRIN=/data1/dussin/STORAGE2/ERAinterim/ALL
var=q2

DIROUT=/data2/dussin/WORK/PREPA_GRD100/forcages/CORRECTION_ARCTIC/$var

cd $DIRIN

listfiles=$( ls | grep $var | grep .nc )

for file in $listfiles ; do

   #############################################################
   # date stuff
   #############################################################

    year=$( echo $file | sed -e "s/y/ /" -e "s/.nc/ /" | awk '{ print $NF }' )
    echo processing year $year

    # dealing with leaping years
    leapyear=0 # by default
    set -A dayinmonth 00 31 28 31 30 31 30 31 31 30 31 30 31 

    leaptest=$( echo $year 2000 | awk '{print ($1 - $2)%4}' )
    if (( $leaptest == 0 )) ; then
        echo $year is a leap year
        leapyear=1
        set -A dayinmonth 00 31 29 31 30 31 30 31 31 30 31 30 31 
    fi

   #############################################################
   # doing monthly means
   #############################################################

   firstindex=1 ; lastindex=0 ; nspd=8

   for month in $( seq 1 12 ) ; do

       # update lastindex
       lastindex=$(( $lastindex + ( nspd * ${dayinmonth[${month}]} ) ))

       mm=$( printf "%02d" $month )
       fileout=$( echo $file | sed -e "s/y${year}/y${year}m${mm}/" )
       ncra -F -O -d time,$firstindex,$lastindex $file -o $DIROUT/$fileout
       echo $fileout done

       # update firstindex
       firstindex=$(( $lastindex + 1 ))

   done

   ## concatenate
   files=$( echo $fileout | sed -e "s/12/??/" )
   filecat=$( echo $fileout | sed -e "s/m12/_monthly/" )
   echo ncrcat -F $DIROUT/$files -o $DIROUT/$filecat
   ncrcat -F $DIROUT/$files -o $DIROUT/$filecat
   rm $DIROUT/$files

done
