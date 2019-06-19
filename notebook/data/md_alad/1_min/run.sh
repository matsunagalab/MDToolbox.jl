#!/bin/bash -x

pmemd.cuda -O \
 -p ../0_init/alad_water.prmtop \
 -i run.in \
 -c ../0_init/alad_water.crd \
 -ref ../0_init/alad_water.crd \
 -o run.out \
 -r run.rst \
 -x run.nc


