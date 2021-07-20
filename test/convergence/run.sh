#!/bin/bash
for i in 32 64 128 256 512 1024
do
    for j in 5.0e-5 2.0e-5 1.0e-5 5.0e-6 2.0e-6 1.0e-6 5.0e-7 2.0e-7
    do
        python gen_params.py $i $j
        python gen_restart.py
        ../../ROME
        python check.py >> converge.dat
    done
done