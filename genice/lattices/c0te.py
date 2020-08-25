#!/usr/bin/python
"""
Usage: genice c0te
"""

from logging import getLogger
import numpy as np

def usage():
    logger = getLogger()
    logger.info(__doc__)

desc={"ref": {"C0": "Page 11 of the Supplemenrary Material of P. Teeratchanan and A. Hermann, Computational phase diagrams of noble gas hydrates under pressure, J. Chem. Phys. 143, 154507 (2015); https://doi.org/10.1063/1.4933371",
},
      "usage": usage(),
      "brief": "Filled ice C0 by Teeratchanan (Hydrogen-disordered.) (Positions of guests are supplied.)"
      }




from genice.cell import cellvectors
import genice.lattices

class Lattice(genice.lattices.Lattice):
    def __init__(self):
        logger = getLogger()

        # Ref. 2atom
        atoms="""
        O1 0.2342 0.4721 0.8019
        O2 0.7648 0.5306 0.2941
        Ne1 -0.0647 0.7868 0.7669
        """

        # Ref.
        # space group: P3_2
        symops="""
          x,            y,            z
         -y,x-y,z+2/3
         -x+y,-x,z+1/3
        """.replace(',', ' ')

        # Ref. 2cell
        a=6.177 / 10.0 #nm
        c=6.054 / 10.0 #nm
        C=120.0

        self.cell  = cellvectors(a,a,c,C=C)

        # helper routines to make from CIF-like data
        from genice import CIF
        atomd = CIF.atomdic(atoms)
        atoms = CIF.fullatoms(atomd, CIF.symmetry_operators(symops))

        self.cagetype = []
        self.cagepos  = []
        for atomname, pos in atoms:
            if atomname == "Ne1":
                self.cagetype.append(atomname)
                self.cagepos.append(pos)

        self.waters, self.pairs = CIF.waters_and_pairs(self.cell, atomd, CIF.symmetry_operators(symops))

        self.density = 18*len(self.waters)/6.022e23 / (np.linalg.det(self.cell)*1e-21)
        self.coord = "relative"
