#!/bin/ksh

INDIR=/fsnet/data/meom/DATA_SET/FORCING_ATMOSPHERIQUE/ERA_INTERIM/drowned

ncra -F $INDIR/drowned_u10_ERAinterim_y????.nc -o mean_u10.nc

ncra -F $INDIR/drowned_v10_ERAinterim_y????.nc -o mean_v10.nc
