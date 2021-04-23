#!/bin/bash

program_at=/Users/yasu/genesis/bin/atdyn
export OMP_NUM_THREADS=1
mpirun -np 8 ${program_at} tmp.inp | tee tmd.log

