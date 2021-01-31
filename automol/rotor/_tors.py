"""
 Build the rotor and torsion objects

 Rotor: (tors_obj_1, tors_obj_2, tors_obj_3)
"""

# import numpy
import automol
from automol.rotor._util import sort_tors_names
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
        # self.span = span(symmetry)
        # self.indices = Torsion._indices(zma, name)

        # Attributes defaulted to none
        self.pot = None
        self.grid = None

    # @staticmethod
    # def span(symmetry):
    #     """ Obtain the torsional span
    #     """
    #     return (2.0 * numpy.pi) / symmetry

    # @staticmethod
    # def _indices(zma, name):
    #     """ Build indices for the torsion
    #     """
    #     mode_idxs = automol.zmat.coord_idxs(zma, name)
    #     mode_idxs = tuple((idx+1 for idx in mode_idxs))


# Build a converter from object to a dictionary


def torsion_lst(zma, gra, lin_keys):
    """  Build a list of torsion objects
    """

    # Build the torsion objects
    _name_axis_dct = name_axis_dct(zma, gra, lin_keys)

    # Get the sorted tors names for building the list
    sorted_tors_names = sort_tors_names(tuple(_name_axis_dct.keys()))

    tors_obj_lst = ()
    for name in sorted_tors_names:
        axis = _name_axis_dct[name]
        grps = torsion_groups(gra, axis)
        symm = torsion_symmetry(gra, axis, lin_keys)
        tors_obj_lst += (Torsion(zma, name, axis, grps, symm),)

    return tors_obj_lst


def name_axis_dct(zma, gra, lin_keys):
    """ Generate the bond keys for the torsion

        or just get the torsion names and keys?
        build a dictionary
        (build full dct for all tors in rotor
         split dcts into subdcts for groupings, including single torsion)
    """

    tors_axes = all_torsion_axes(gra, lin_keys)

    tors_names = tuple(automol.zmat.torsion_coordinate_name(zma, *keys)
                       for keys in tors_axes)

    return dict(zip(tors_names, tors_axes))


def all_torsion_axes(gra, lin_keys):
    """ Build the torsion axes
    """
    return automol.graph.rotational_bond_keys(gra, lin_keys=lin_keys)


def all_torsion_groups(gra, lin_keys):
    """ Generate torsion groups make generalizable to multiple axes
    """

    axes = all_torsion_axes(gra, lin_keys)
    grps = ()
    for axis in axes:
        grps += (torsion_groups(gra, axis),)

    return grps


def all_torsion_symmetries(gra, lin_keys):
    """ Generate torsion groups make generalizable to multiple axes
    """

    axes = all_torsion_axes(gra, lin_keys)
    syms = ()
    for axis in axes:
        syms += (torsion_symmetry(gra, axis, lin_keys),)

    return syms


def torsion_groups(gra, axis):
    """ Generate torsion groups make generalizable to multiple axes
    """
    return automol.graph.rotational_groups(gra, *axis)


def torsion_symmetry(gra, axis, lin_keys):
    """ Obtain the symmetry number for the torsion
    """
    return automol.graph.rotational_symmetry_number(
        gra, *axis, lin_keys=lin_keys)
