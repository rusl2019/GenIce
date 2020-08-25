#!/usr/bin/python
"""
Usage: genice c2te
"""
from logging import getLogger
import numpy as np

def usage():
    logger = getLogger()
    logger.info(__doc__)

desc={"ref": {"C2": "Page 12 of the Supplemenrary Material of P. Teeratchanan and A. Hermann, Computational phase diagrams of noble gas hydrates under pressure, J. Chem. Phys. 143, 154507 (2015); https://doi.org/10.1063/1.4933371",
},
      "usage": usage(),
      "brief": "Filled ice C2 (cubic ice) by Teeratchanan (Hydrogen disordered). (Positions of guests are supplied.)"
      }

def pick_atoms(atoms, names, repeat=(1,1,1)):
    nrep = np.array(repeat)
    for atomname, fracpos in atoms:
        if atomname in names:
            for x in range(repeat[0]):
                for y in range(repeat[1]):
                    for z in range(repeat[2]):
                        yield atomname, (fracpos+np.array([x,y,z]))/nrep



import genice.lattices

class Lattice(genice.lattices.Lattice):
    def __init__(self):
        logger = getLogger()

        # Ref. C2
        atoms="""
    O1 0.0000 0.0000 0.1937
    Ne1 0.0000 0.0000 0.7066
        """

        # Ref. C2
        # space group: I4_1md No. 109 http://img.chem.ucl.ac.uk/sgp/large/109az1.htm
        symops="""
          x,            y,            z
         -x,           -y,            z
         -y,          1/2+x,        1/4+z
          y,          1/2-x,        1/4+z
         -x,            y,            z
          x,           -y,            z
          y,          1/2+x,        1/4+z
         -y,          1/2-x,        1/4+z

          x+1/2,            y+1/2,            z+1/2
         -x+1/2,           -y+1/2,            z+1/2
         -y+1/2,          1/2+x+1/2,        1/4+z+1/2
          y+1/2,          1/2-x+1/2,        1/4+z+1/2
         -x+1/2,            y+1/2,            z+1/2
          x+1/2,           -y+1/2,            z+1/2
          y+1/2,          1/2+x+1/2,        1/4+z+1/2
         -y+1/2,          1/2-x+1/2,        1/4+z+1/2
        """.replace(',', ' ')

        # Ref. C2
        a=4.409 / 10.0 #nm
        c=6.251 / 10.0 #nm

        from genice.cell import cellvectors
        self.cell  = cellvectors(a,a,c)

        # helper routines to make from CIF-like data
        from genice import CIF
        atomd = CIF.atomdic(atoms)
        atoms = CIF.fullatoms(atomd, CIF.symmetry_operators(symops))

        self.cagetype = []
        self.cagepos  = []
        for name, pos in pick_atoms(atoms, ("Ne1",), repeat=(2,2,2)):
            self.cagetype.append(name)
            self.cagepos.append(pos)

        self.waters, self.pairs = CIF.waters_and_pairs(self.cell, atomd, CIF.symmetry_operators(symops), rep=(2,2,2))

        self.density = 18*len(self.waters)/6.022e23 / (np.linalg.det(self.cell)*1e-21)
        self.coord = "relative"
