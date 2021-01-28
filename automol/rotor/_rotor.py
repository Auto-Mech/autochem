"""
 Build the rotor and torsion objects
"""

import numpy
import automol
from automol.rotor._par import TorsionParam


# add a grouping function for the rotor

class Rotor:
    """ Describes a Rotor
    """

    def __init__(self, zma, names=None, model=None, rotor_dct=None):
        """ constructor
        """
        torsion_axes = all_torsion_axes(zma, names=names)
        print(_torsion_bond_keys(zma))
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
            self.symmetry = tors_dct.get(TorsionParam.SYMMETRY, None)
            self.span = tors_dct.get(TorsionParam.SPAN, None)
            self.axis = tors_dct.get(TorsionParam.AXIS, None)
            self.groups = tors_dct.get(TorsionParam.GROUPS, None)
        else:
            pass
            # self.symmetry = symmetry_number(self)
            # self.span = span(self)
            # self.axis, self.groups = rotational_groups(self)


def build_torsion_axes(zma, names=None):
    """ Generate the bond keys for the torsion

        or just get the torsion names and keys?
        build a dictionary
        (build full dct for all tors in rotor
         split dcts into subdcts for groupings, including single torsion)
    """

    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    gra = automol.geom.graph(geo)

    lin_keys = sorted(gdummy_key_dct.keys())
    tors_keys = automol.graph.rotational_bond_keys(gra, lin_keys=lin_keys)

    # Determine all of the torsion names, take subset if requested
    tors_names = tuple(automol.zmat.torsion_coordinate_name(zma, *keys)
                       for keys in tors_keys)
    if names is not None:
        tors_names = tuple(tors_names for name in tors_names
                           if name in names)

    return dict(zip(tors_names, tors_keys))


def torsion_groups(gra, axis):
    """ Generate torsion groups
    """
    return automol.graph.rotational_groups(gra, *axis)


def torsion_symmetry(gra, axis, lin_keys):
    """ Obtain the symmetry number for the torsion
    """
    return automol.graph.rotational_symmetry_number(
        gra, *axis, lin_keys=lin_keys)


def torsion_span(gra, axis, lin_keys):
    """ Obtain the torsional span
    """
    sym_num = automol.graph.rotational_symmetry_number(
        gra, *axis, lin_keys=lin_keys)
    return (2.0 * numpy.pi) / sym_num
