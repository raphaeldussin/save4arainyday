#!/bin/ksh

for file in $( ls | grep drowned_radlw ) ; do

    fileout=$( echo $file | sed -e "s/DFS5.1/DFS5.1.1/g")
    mv $file $fileout

done

