#!/bin/bash
for i in 32 64 128 256 512
do
    for j in 0.7 1.0 1.2 1.5 1.7 2.0
    do
        # Note, to use the gen_restart.py tool, first install the 
        # f90nml, numpy, and netCDF4 python modules
        python gen_params.py $i $j
        python gen_restart_p.py
        ../../ROME
        python check.py >> converge.dat
    done
done