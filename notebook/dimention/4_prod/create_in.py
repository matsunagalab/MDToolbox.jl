#!/usr/bin/python
# coding: utf-8

import random

############# define functions
def print_lines(f, phi, psi):
    ig = random.randint(0,1000000)
    f.write("production run\n")
    f.write(" &cntrl\n")
    f.write("   ig=%d, \n" % (ig))
    f.write("   irest=1, ntx=5,\n")
    f.write("   nstlim=50000, dt=0.002,\n")
    f.write("   ntc=2, ntf=1, tol=0.000001,\n")
    f.write("   ntb=2, nscm=50000,\n")
    f.write("   ntt=3, gamma_ln=2.0, temp0=300.0,\n")
    f.write("   ntp=1, taup=5.0,\n")
    f.write("   cut=10.0,\n")
    f.write("   ioutfm=1,\n")
    f.write("   ntpr=100, ntwx=100, ntwr=50000, iwrap=1,\n")
    f.write("   nmropt=1,\n")
    f.write(" /\n")
    f.write(" &wt\n")
    f.write("  type='DUMPFREQ', istep1=100,\n")
    f.write(" /\n")
    f.write(" &wt\n")
    f.write("  type='END',\n")
    f.write(" /\n")
    f.write("DISANG=run_%d_%d.disang\n" % (phi, psi))
    f.write("DUMPAVE=run_%d_%d.dat\n" % (phi, psi))
    f.write("\n")

############# main
phi = range(-180, 0, 15)
psi = range(0, 180, 15)

for i in phi:
    for j in psi:
        filename = "run_%d_%d.in" % (i, j)
        print "writing %s..." % (filename)
        f = open(filename, 'w')
        print_lines(f, i, j)
        f.close()

