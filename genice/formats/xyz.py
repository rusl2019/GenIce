# coding: utf-8
"""
XMol file format (.xyz)
"""

import numpy as np

from genice import rigid
from logging import getLogger

import genice.formats
class Format(genice.formats.Format):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {7:self.hook7}


    def hook7(self, ice):
        logger = getLogger()
        logger.info("Hook7: Output in XYZ format.")
        s = ""
        s += "{0}\n".format(len(ice.atoms))
        s += "Generated by GenIce. https://github.com/vitroid/GenIce\n"
        for atom in ice.atoms:
            molorder, resname, atomname, position, order = atom
            s += "{0:5} {1:9.4f} {2:9.4f} {3:9.4f}\n".format(atomname,position[0]*10,position[1]*10,position[2]*10)
        s = '#' + "\n#".join(ice.doc) + "\n" + s
        print(s,end="")
        logger.info("Hook7: end.")
