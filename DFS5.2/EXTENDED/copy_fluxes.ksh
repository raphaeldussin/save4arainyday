#!/bin/ksh

for var in radsw radlw precip snow ; do

    sudo -u viking rsync -av drowned_${var}* /fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/DFS5.2/extension_1958-1978/$var/.

done
