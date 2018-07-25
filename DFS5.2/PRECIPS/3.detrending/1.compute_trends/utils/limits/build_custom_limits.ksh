### FARC : utility to build a limits.custom file from list of variables          ###
### Advise : rename your limits.custom files to avoid overwriting and comment it ###

listvar='t2 q2 u10 v10 radlw radsw precip snow Qnet Qtrb Qlat Qsen Qrad Qlw Qsw U10 Taux Tauy EmP evap'

rm limits.custom

for var in $listvar ; do

    cat limits.template | sed -e "s/<VAR>/$var/g" >> limits.custom

done
