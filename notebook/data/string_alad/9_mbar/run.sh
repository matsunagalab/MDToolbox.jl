#!/bin/bash

program_at=/Users/yasu/genesis/bin/mbar_analysis
export OMP_NUM_THREADS=4
${program_at} mbar.inp | tee mbar.log

