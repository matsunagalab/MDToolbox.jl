# Input and Output

Typical usages of I/O functions for MD files are as follows. 

## PDB files

```julia
# read pdb data and generate a new variable t whose type is TrjArray
t = mdload("protein.pdb")

# after some editings
mdsave("protein_edited.pdb", t)
```

## AMBER files

AMBER NetCDF trajectory file
```julia
# pdb files can be used for obtaining topology information which is used in atom selections
t = mdload("run.pdb")

# read AMBER NetCDF data and generate a new variable t whose type is TrjArray
# the topology information can be attached as the option `top=t`
t = mdload("run.nc", t=top)

# after some calculations
mdsave("run_edited.nc", t)
```

## CHARMM/NAMD files
```julia
# read PSf data and generate a new variable t whose type is TrjArray
# Although PSF file does not contain coordinates, topology information is used in atom selections
t = mdload("run.psf")
# or if you don't have PSF files, pdb can be used to obtain topology information
t = mdload("run.pdb")

# read dcd file and generate a new variable t whose type is TrjArray
# the topology information can be attached as the option `top=t`
t = mdload("run.dcd", top=t)
# after some calculations
mdsave("run_edit.dcd", t)
```

## GROMACS files

```
not available yet
```
