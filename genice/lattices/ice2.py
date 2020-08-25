#!/usr/bin/python
"""
Usage: genice two
"""
from logging import getLogger
from genice import CIF
import numpy as np
from genice.cell import cellvectors
import genice.lattices


def usage():
    logger = getLogger()
    logger.info(__doc__)

desc={"ref": {"2atom": "Kamb, B., Hamilton, W. C., LaPlaca, S. J. & Prakash, A. Ordered Proton Configuration in Ice II, from Single‐Crystal Neutron Diffraction. J. Chem. Phys. 55, 1934–1945 (2003).",
              "2cell": "Kamb, B.IUCr. Ice. II. A proton-ordered form of ice. Acta Cryst 17, 1437–1449 (1964).",
              "C1": "D. Londono, W. F. Kuhs, J. L. Finney, Nature, 1988, DOI:10.1038/332141a0."},
      "usage": usage(),
      "brief": "Hydrogen-ordered ice II."
      }



class Lattice(genice.lattices.Lattice):
    def __init__(self):
        logger = getLogger()

        # Ref. 2atom
        atoms="""
        O1   0.2716    0.0259   -0.1471
        O2   0.4798    0.7571    0.3389
        D1   0.7284    0.4038    0.4034
        D2   0.1491    0.0412   -0.2023
        D3   0.7420    0.1978    0.3708
        D4   0.4232    0.1954   -0.0164
        """

        # Ref. 2cell
        # space group: R-3
        symops="""
          x,            y,            z
          z,            x,            y
          y,            z,            x

         -x,           -y,           -z
         -z,           -x,           -y
         -y,           -z,           -x
        """.replace(',', ' ')

        # Ref. 2cell
        a=7.78 / 10.0 #nm
        A=113.1

        self.cell  = cellvectors(a,a,a,A,A,A)

        # helper routines to make from CIF-like data
        atomd = CIF.atomdic(atoms)
        sops  = CIF.symmetry_operators(symops)
        self.waters, self.fixed = CIF.waters_and_pairs(self.cell, atomd, sops)

        # set pairs in this way for hydrogen-ordered ices.
        self.pairs = self.fixed

        self.density = 18*len(self.waters)/6.022e23 / (np.linalg.det(self.cell)*1e-21)

        self.coord = "relative"
