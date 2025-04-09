#!/bin/bash



read_backup_file=$1

fname=`basename -s .sav $read_backup_file`

for i in 31 32 33 34 35 36 38 39 40 41;
do
    root -b -q 'phoebus_default_ecal_hv.C('$i',"'$read_backup_file'")'
    chmod +x 'HV_inputs_hv_crate'$i'_'$fname'.sh'
    chmod +x 'power_on_HV_channels_crate'$i'_'$fname'.sh'
done


    
