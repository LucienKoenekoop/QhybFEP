import sys

"""
Module for setting up scripts for qprep, qdyn and qfep.
"""

def qprep(hybrid_lib, ff_lib, hybrid_prm, hybrid_pdb, center, radius, solvent, cyx):
    """
    Function to set up qprep input file for reading
    """
    with open('tmp/qprep.inp', 'w') as qprep_inp:
        qprep_inp.write(
f"""\
rl {hybrid_lib}
rl {ff_lib}
rprm {hybrid_prm}
rp {hybrid_pdb}
boundary 1 {center} {radius}
solvate {center} {radius} 1 {solvent}
! addbond {cyx[0]} {cyx[1]}
maketop hybrid_top
writetop tmp/hybrid_top.top
wp tmp/hybrid_top.pdb y
rt tmp/hybrid_top.top
mask none
mask not excluded
wp tmp/not_excluded.pdb y
q
"""
        )
        qprep_inp.close()
    return

