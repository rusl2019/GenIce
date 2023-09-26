# coding: utf-8

from genice2.molecules import serialize
import genice2.formats
from genice2.decorators import timeit, banner
from collections import defaultdict

desc = {
    "ref": {},
    "brief": "Povray2.",
    "usage": """
Usage: genice2 icename -f povray2

options:
    No options available.
""",
}


def Block(name, content):
    return " " + name + " { " + content + " } "


def Vector(v):
    return "<{0:.3f}, {1:.3f}, {2:.3f}> ".format(*v)


def Juxtapose(v):
    return ",".join(v)


def Atom(atomtype, pos):
    return (
        Block(
            "sphere",
            Juxtapose([Vector(pos), "R{0}".format(atomtype)])
            + Block("material", "MAT{0}".format(atomtype)),
        )
        + "\n"
    )


def Bond(bondtype, pos1, pos2):
    return (
        Block(
            "cylinder",
            Juxtapose([Vector(pos1), Vector(pos2), "R{0}".format(bondtype)])
            + Block("material", "MAT{0}".format(bondtype)),
        )
        + "\n"
    )


def Include(filename):
    return '#include "{0}"\n'.format(filename)


class Format(genice2.formats.Format):
    """
    The atomic positions of the molecules are output in Povray format.
    No options available.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {7: self.Hook7, 6: self.Hook6}

    @timeit
    @banner
    def Hook6(self, ice):
        "Output water molecules in Povray format."
        atoms = []
        for mols in ice.universe:
            atoms += serialize(mols)

        # prepare the reverse dict

        waters = defaultdict(dict)
        for atom in atoms:
            resno, resname, atomname, position, order = atom
            if "O" in atomname:
                waters[order]["O"] = position
            elif "H" in atomname:
                if "H0" not in waters[order]:
                    waters[order]["H0"] = position
                else:
                    waters[order]["H1"] = position
        s = Include("default.inc")
        s += "union {\n"
        for order, water in waters.items():
            O = water["O"]
            s += Atom("O", O)
        for i, j in ice.spacegraph.edges(data=False):
            if i in waters and j in waters:  # edge may connect to the dopant
                O1 = waters[j]["O"]
                O2 = waters[i]["O"]
                d0 = O2 - O1
                rr0 = d0 @ d0
                if rr0 < 0.3**2:
                    s += Bond("HB", O1, O2)
        self.output = s

    @timeit
    @banner
    def Hook7(self, ice):
        "Output guest molecules in Povray format."
        gatoms = []
        for mols in ice.universe[1:]:
            gatoms += serialize(mols)
        cellmat = ice.repcell.mat
        s = ""
        H = []
        O = ""
        for atom in gatoms:
            resno, resname, atomname, position, order = atom
            s += Atom(atomname, position)
        s = "//" + "\n//".join(ice.doc) + "\n" + s
        s += (
            "  translate "
            + Vector(-(cellmat[0, :] + cellmat[1, :] + cellmat[2, :]) / 2)
            + "\n}\n\n"
        )
        self.output += s
