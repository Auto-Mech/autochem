"""
 Build the rotor and torsion objects
"""

import automol
from automol.rotor._par import TorsionParam


class Rotor:
    """ Describes a Rotor
    """

    def __init__(self, zma, names=None, model=None, rotor_dct=None):
        """ constructor
        """
        self.names = names
        self.zma = zma

        torsions = ()
        for name in names:

    # Use input rotor dct to construct all te torsion objects if there
    # else
    # Build the dict of tors  names {name: bnd_idx}
    # Build the torsion objects for each name
    # Group the names nad torsions according to model
    # Rotor is a dct of torsion objects maybe?


class Torsion:
    """ Describes a torsion
    """

    def __init__(self, name, zma, tors_dct=None):
        """ constructor
        """
        self.name = name
        self.zma = zma
        if tors_dct is not None:
            self.symmetry = tors_dct.get(TorsionParam.ZMA, None)
            self.span = tors_dct.get(TorsionParam.ZMA, None)
            self.axis = tors_dct.get(TorsionParam.ZMA, None)
            self.groups = tors_dct.get(TorsionParam.ZMA, None)
        else:
            self.symmetry = symmetry_number(self)
            self.span = span(self)
            self.axis, self.groups = rotational_groups(self)


def torsion_axis():

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(ggra, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))


def torsion_groups(gra, axis=None):

    if axis is None:
        axis = torsion_axis()

    return automol.graph.rotational_groups(ggra, axis)


def torsion_symmetry(gra, axis=None):
    """ Obtain the symmetry number for the torsion
    """

    if axis is None:
        axis = torsion_axis()

    return automol.graph.rotational_symmetry_number(
        ggra, axis, lin_keys=lin_keys)


def span(sym_num):
    """
    """
    return (2.0 * numpy.pi) / sym_num


def _torsion_bond_keys(zma, tors_name):
    """ Generate the bond keys for the torsion

        or just get the torsion names and keys?
        build a dictionary (build full dct for all tors in rotor
                            split dcts into subdcts for groupings (including single torsion)
    """

    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    gra = automol.geom.graph(geo)

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(ggra, lin_keys=lin_keys)

    # names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    for keys in gbnd_keys:
        name = automol.zmat.torsion_coordinate_name(zma, keys)


