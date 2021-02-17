"""
 Build the rotor and torsion objects.

 Rotor: (tors_obj_1, tors_obj_2, tors_obj_3)
"""

import copy
import numpy
import automol
from automol.rotor._util import sort_tors_names
# from automol.rotor._par import TorsionParam


class Torsion:
    """ Describes a torsion, which one or more make up a rotor
    """

    def __init__(self, zma, name, axis, groups, symmetry,
                 indices=None):
        """ constructor
        """

        self.name = name
        self.zma = zma
        self.symmetry = symmetry
        self.groups = groups
        self.axis = axis
        self.set_span()
        if indices is None:
            self.set_indices()
        else:
            self.indices = indices
        self.pot = None

    def set_span(self):
        """ Obtain the torsional span
        """
        self.span = (2.0 * numpy.pi) / self.symmetry

    def set_indices(self):
        """ Build indices for the torsion
        """
        self.indices = automol.zmat.coord_idxs(self.zma, self.name)

    def copy(self):
        """ return a copy of this Reaction
        """
        return Torsion(
            self.zma, self.name, self.axis, self.groups,
            self.symmetry, self.span, self.indices)


# Build functions
def torsion_lst(zma, gra, lin_keys):
    """  Build a list of torsion objects
    """

    # Build the torsion objects
    _name_axis_dct = name_axis_dct(zma, gra, lin_keys)

    # Get the sorted tors names for building the list
    sorted_tors_names = sort_tors_names(tuple(_name_axis_dct.keys()))

    tors_obj_lst = ()
    for name in sorted_tors_names:

        # Determine constituent rotor pieces in graph system
        axis = _name_axis_dct[name]
        grps = torsion_groups(gra, axis)
        symm = torsion_symmetry(gra, axis, lin_keys)

        # Shift the axis and groups to be in the zma system
        # zaxis = automol.util.dummy.shift_up(zma, axis)
        # zgrps = tuple(automol.util.dummy.shift_up(zma, grp)
        #               for grp in grps)

        # Build the torsion object and add to the list
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
    print('tors_axes', tors_axes)
    tors_names = tuple(automol.zmat.torsion_coordinate_name(zma, *keys)
                       for keys in tors_axes)

    return dict(zip(tors_names, tors_axes))


def all_torsion_axes(gra, lin_keys):
    """ Build the torsion axes
    """
    tors_keys = automol.graph.rotational_bond_keys(gra, lin_keys=lin_keys)
    return tuple(tuple(keys) for keys in tors_keys)


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


# Manipulate the torsion objects
def copy_pot(pot):
    """ dumb function to copy pot
    """
    return copy.deepcopy(pot)


def relabel_for_geometry(torsion):
    """ relabel the torsion objec tto correspond with a geometry converted
        from a z-matrix
    """

    name = torsion.name
    zma = torsion.zma
    symmetry = torsion.symmetry
    ggrps = tuple(automol.util.dummy.shift_down(zma, grp)
                  for grp in torsion.groups)
    gaxis = automol.util.dummy.shift_down(zma, torsion.axis)
    gindices = automol.util.dummy.shift_down(zma, torsion.indices)

    gtorsion = Torsion(zma, name, gaxis, ggrps, symmetry, indices=gindices)
    gtorsion.pot = copy_pot(torsion.pot)

    return gtorsion
