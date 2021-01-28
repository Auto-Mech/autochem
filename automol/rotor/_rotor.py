"""
 Build the rotor and torsion objects

 Rotor: (tors_obj_1, tors_obj_2, tors_obj_3)
"""

import numpy
import automol
# from automol.rotor._par import TorsionParam


class Torsion:
    """ Describes a torsion, which one or more make up a rotor
    """

    def __init__(self, zma, name, axis, groups, symmetry):
        """ constructor
        """

        self.name = name
        self.zma = zma
        self.symmetry = symmetry
        self.groups = groups
        self.axis = axis
        self.span = Torsion.span(symmetry)
        self.indices = Torsion._indices(zma, name)

        # Attributes defaulted to none
        self.pot = None
        self.grid = None


    def span(symmetry):
        """ Obtain the torsional span
        """
        return (2.0 * numpy.pi) / symmetry


    def _indices(zma, name):
        """ Build indices for the torsion
        """
        mode_idxs = automol.zmat.coord_idxs(zma, name)
        mode_idxs = tuple((idx+1 for idx in mode_idxs))


# constructors
def from_zma(zma, names=None, model=None):
    """ Construct a list-of-lists of torsion objects
    """

    # add a grouping function for the rotor

    # Build a graph that is used to get torsion object info
    gra, lin_keys = _graph(zma)

    # Build the torsion objects
    tors_name_dct = torsion_name_dct(zma, gra, lin_keys, names=names)

    rotors = ()
    for name_grp in name_grps:
        torsions = ()
        for name in name_grp:
            axis = tors_name_dct[name]
            groups = torsion_groups(gra, axis)
            symmetry = torsion_symmetry(gra, axis, lin_keys)
            torsions += (Torsion(zma, name, axis, groups, symmetry),)
        rotors += (torsions,)


# def from_data(tors_dct, names=None):
#     """ Build the rotors objects from existing data
#     """
#     assert set(tors_dct.keys) >= {'sym', 'span', 'group'}


# Build fundamental info needed to build torsions from basis automol structures
def torsion_name_dct(zma, gra, lin_keys, names=None):
    """ Generate the bond keys for the torsion

        or just get the torsion names and keys?
        build a dictionary
        (build full dct for all tors in rotor
         split dcts into subdcts for groupings, including single torsion)
    """

    tors_keys = automol.graph.rotational_bond_keys(gra, lin_keys=lin_keys)

    # Determine all of the torsion names, take subset if requested
    tors_names = tuple(automol.zmat.torsion_coordinate_name(zma, *keys)
                       for keys in tors_keys)
    if names is not None:
        tors_names = tuple(tors_names for name in tors_names
                           if name in names)

    return dict(zip(tors_names, tors_keys))


def torsion_axes(gra, lin_keys):
    """ Build the torsion axes
    """
    return automol.graph.rotational_bond_keys(gra, lin_keys=lin_keys)


def torsion_groups(gra, axis):
    """ Generate torsion groups make generalizable to multiple axes
    """
    return automol.graph.rotational_groups(gra, *axis)


def torsion_symmetry(gra, axis, lin_keys):
    """ Obtain the symmetry number for the torsion
    """
    return automol.graph.rotational_symmetry_number(
        gra, *axis, lin_keys=lin_keys)


# Helpers
def _graph(zma):
    """ Generate the graph
    """

    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    gra = automol.geom.graph(geo)
    lin_keys = sorted(gdummy_key_dct.keys())

    return gra, lin_keys
