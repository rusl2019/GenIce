# coding: utf-8

desc = { "ref": {},
         "brief": "Cell-reshaper.",
         "usage": """
A formatter plugin for GenIce to produce a python lattice plugin.

Usage:
    genice ice5 -f reshape[1,0,0,1,1,0,0,0,1] > reshaped_ice5.py

Options:
    i,j,k,l,m,n,o,p,q  Nine integers separated by commas.

    Cell vectors of the reshaped cell, a', b' and c' is calculated by linear combinations of the original cell vectors, a, b, and c, as

        a' = i a + j b + k c
        b' = l a + m b + n c
        c' = o a + p b + q c

    Water molecules are relocated appropriately.
""" }

import numpy as np
import itertools as it
from math import floor
import re
from logging import getLogger
import genice.formats


MAXCELL=11
torr = 1e-8


def isZero(x):
    return -torr < x < torr



def FlagEquivCells(nv, flags, ijk):
    if nv not in flags:
        for d in it.product((-2,-1,0,1,2), repeat=3):
            vv = tuple(np.dot(d,ijk)+nv)
            flags.add(vv)


def FindEmptyCells(cellmat, ijk, relpositions, labels=None):
    newcell = np.dot(ijk, cellmat)
    newcelli = np.linalg.inv(newcell)
    L1 = np.linalg.norm(newcell[0])
    L2 = np.linalg.norm(newcell[1])
    L3 = np.linalg.norm(newcell[2])
    rL = np.array([L1,L2,L3])
    # print(rL)
    ncell = 0
    flags = set()
    queue = [(0,0,0)]
    s = ""
    while len(queue) > 0:
        nv = queue.pop(0)
        if nv not in flags:
            ncell += 1
            # Place atoms
            for i, xyz in enumerate(relpositions):
                # rel to abs
                xxv = np.dot(xyz + nv, cellmat)
                # inner products with axis vector of the newcell
                pv = np.dot(xxv, newcelli)
                # print(pv,rL)
                pv -= np.floor(pv)
                # print(nv)
                label = ""
                if labels is not None:
                    label = labels[i]
                s += "{3} {0:9.4f} {1:9.4f} {2:9.4f}\n".format(pv[0],pv[1],pv[2], label)
            #print("NV:",nv)
            FlagEquivCells(nv, flags, ijk)
            for xi,yi,zi in it.product((-1,0,1), repeat=3):
                nei = nv[0]+xi,nv[1]+yi,nv[2]+zi
                if nei not in flags:
                    #print("nei", nei)
                    queue.append(nei)
    return ncell, s



class Format(genice.formats.Format):
    ijk = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])


    def __init__(self, **kwargs):
        unknown = dict()
        for k, v in kwargs.items():
            if re.match("^[-+0-9,]+$", k) is not None and v is True:
                self.ijk = np.array([int(x) for x in k.split(",")]).reshape(3,3)
            else:
                unknown[k] = v
        super().__init__(**unknown)


    def hooks(self):
        return {1:self.hook1}


    def hook1(self, lattice):
        logger = getLogger()
        logger.info("Hook1: Output as a python module.")
        # Original cell matrix.
        cellmat = lattice.repcell.mat
        logger.info("  Reshaping the unit cell.")
        logger.info("    i:{0}".format(self.ijk[0]))
        logger.info("    j:{0}".format(self.ijk[1]))
        logger.info("    k:{0}".format(self.ijk[2]))
        # reshaped cell might not be rect.
        newcell = np.dot(self.ijk, cellmat)
        # replication ratio.
        vol = abs(np.linalg.det(self.ijk))
        vol = floor(vol*8192+0.5)/8192
        # Unit orthogonal vectors for the newcell.
        e1 = newcell[0].astype(float)
        e1 /= np.linalg.norm(e1)
        e3 = np.cross(newcell[0], newcell[1]).astype(float)
        e3 /= np.linalg.norm(e3)
        e2 = np.cross(e3,e1)
        R = np.array([e1,e2,e3])
        RI = np.linalg.inv(R)
        # Regularization
        # Let x axis of the newcell be along e1
        # and let y axis of the newcell be on the e1-e2 plane.
        regcell = np.dot(newcell, RI)
        a = (-torr < regcell)
        b = (regcell < torr)
        c = a & b
        regcell[c] = 0.0

        # header
        s = ""
        s += '"""\n'
        s += "\n".join(lattice.doc) + "\n"
        s += "Reshaping the unit cell.\n"
        s += "  i:{0}\n".format(self.ijk[0])
        s += "  j:{0}\n".format(self.ijk[1])
        s += "  k:{0}\n".format(self.ijk[2])
        s += '"""\n'

        s += "bondlen={0}\n".format(lattice.bondlen*10)
        s += "coord='relative'\n"
        if isZero(regcell[1,0]) and isZero(regcell[2,0]) and isZero(regcell[2,1]):
            s += "from genice.cell import cellvectors\n"
            s += "cell = cellvectors(a={0:.8f}, b={1:.8f}, c={2:.8f})\n".format(regcell[0,0]*10,regcell[1,1]*10,regcell[2,2]*10)
        else:
            s += "\nimport numpy as np\n"
            s += "cell=np.array(["
            for d in range(3):
                s += "[{0:.8f}, {1:.8f}, {2:.8f}], ".format(regcell[d,0]*10,regcell[d,1]*10,regcell[d,2]*10)
            s += "])\n"
        # s += "cell='{0} {1} {2}'\n".format(ri,rj,rk)
        s += "density={0}\n".format(lattice.density)
        s += 'waters="""'+"\n"
        logger.info("  Total number of molecules: {0}".format(vol*len(lattice.reppositions)))

        ncell, ss = FindEmptyCells(cellmat, self.ijk, lattice.reppositions)
        assert vol == ncell

        s += ss + '"""' + "\n\n"

        if lattice.cagepos is not None:
            s += 'cages="""'+"\n"
            ncell, ss = FindEmptyCells(cellmat, self.ijk, lattice.repcagepos, labels=lattice.repcagetype)
            s += ss + '"""'+"\n\n"

        print(s,end="")

        logger.info("Hook1: end.")


    # Reshaping matrix (Must be integers)
    # for now they are hardcoded.  It will be given as options for the plugin in the future.
    # ijk = np.array([[1, 1, 1], [1, -1, 0], [1, 1, -2]])
    # ijk = np.array([[2, 0, 1], [0, 1, 0], [-1, 0, 2]])
    # ijk = np.array([[1, 0, 0], [0, 1, 0], [1, 0, 2]])
    # ijk = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
