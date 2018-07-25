#!/bin/ksh

MYWORKDIR=/fsnet/data/meom/workdir/dussin/TMPDIR_PRECIPS_DFS5.2/ERAinterim_ORCA025
CODES=/home/users/dussin/TOOLS/DFS5-tools

cd $MYWORKDIR

ln -s /fsnet/data/meom/workdir/dussin-adeplacer1/MESH-MASKS/ORCA025-MJM91_byte_mask.nc ./mask.nc
ln -s /fsnet/data/meom/workdir/dussin-adeplacer1/MESH-MASKS/ORCA025-MJM91_mesh_hgr.nc ./mesh_hgr.nc

ln -s /home/users/dussin/DFS5/DFS5.1/STORTO/weights/PMWC_precip_correction.nc .

ln -s /home/users/dussin/DFS5/DFS5.2/PRECIPS/2.apply_storto/ORCA025-ERAint.map .

for year in $(seq 2012 2012 ) ; do

    $CODES/correc_precip_astorto precip_ERAinterim-ORCA025_y$year.nc PMWC_precip_correction.nc -nointerp

    mv precip_corrected.nc precip_ERAinterim_corr-ORCA025_y$year.nc

    cat /home/users/dussin/DFS5/DFS5.2/PRECIPS/2.apply_storto/namelist.precip.step2.skel | sed -e "s/<YEAR>/$year/g" \
        >> namelist.precip.$year

    $HOME/BIN/SOSIE/sosie.x -f namelist.precip.$year

    ln -s $ISLANDSPATH_DATA_SET/FORCING_ATMOSPHERIQUE/ERA_INTERIM/drowned/drowned_precip_ERAinterim_y${year}.nc .

    $CODES/blend_precip precip_ORCA025-ERAint_${year}.nc drowned_precip_ERAinterim_y${year}.nc
    mv output.nc drowned_precip_ERAinterim_corr_astorto_lowlats_y${year}.nc

done

