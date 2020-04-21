package require psfgen
topology ./toppar/top_all22_prot.rtf

segment PRO {pdb alad_woH.pdb}
coordpdb alad_woH.pdb PRO
regenerate angles dihedrals
guesscoord

writepdb sys.pdb
writepsf sys.psf

