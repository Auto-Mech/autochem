"""
 Build the rotor and torsion objects.

 Rotor: (tors_obj_1, tors_obj_2, tors_obj_3)
"""

import copy
import numpy
from automol.rotor._util import sort_tors_names
import automol.graph
import automol.zmat
import automol.reac


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
        self.scan_geos = None

    def set_span(self):
        """ Obtain the torsional span
        """
        self.span = (2.0 * numpy.pi) / self.symmetry

    def set_indices(self):
        """ Build indices for the torsion
        """
        self.indices = automol.zmat.coord_idxs(self.zma, self.name)


# Build functions
def torsion_lst(zma):
    """  Build a list of torsion objects
    """

    # Get the necessary graph and lin keys
    gra = automol.zmat.graph(zma, stereo=True, dummy=True)
    lin_keys = sorted(
        automol.graph.dummy_atoms_parent_key(gra).values())

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

        # Build the torsion object and add to the list
        tors_obj_lst += (Torsion(zma, name, axis, grps, symm),)

    return tors_obj_lst


def reaction_torsion_lst(zma, zrxn):
    """ torsions from zrxn obj
    """

    zbnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
    tors_axes = tuple(tuple(keys) for keys in zbnd_keys)
    tors_names = tuple(automol.zmat.torsion_coordinate_name(zma, *keys)
                       for keys in tors_axes)

    # grxn = automol.reac.relabel_for_geometry(zrxn)
    # gbnd_keys = automol.reac.rotational_bond_keys(grxn)
    # tors_axes = tuple(tuple(keys) for keys in gbnd_keys)

    # Get the sorted tors names for building the list
    name_dct = dict(zip(tors_names, tors_axes))

    tors_obj_lst = ()
    sorted_tors_names = sort_tors_names(tuple(name_dct.keys()))
    for name in sorted_tors_names:

        # gaxis = name_dct[name]
        # gaxis = tuple(sorted(gaxis))
        # ggrps = automol.reac.rotational_groups(grxn, *gaxis)
        # symm = automol.reac.rotational_symmetry_number(grxn, *gaxis)
        axis = name_dct[name]
        axis = tuple(sorted(axis))
        grps = automol.reac.rotational_groups(zrxn, *axis)
        symm = automol.reac.rotational_symmetry_number(zrxn, *axis)

        # Build the torsion object and add to the list
        # tors_obj_lst += (Torsion(zma, name, gaxis, ggrps, symm),)
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
    ggrps = tuple(automol.zmat.shift_down(zma, grp)
                  for grp in torsion.groups)
    gaxis = automol.zmat.shift_down(zma, torsion.axis)
    gindices = automol.zmat.shift_down(zma, torsion.indices)

    gtorsion = Torsion(zma, name, gaxis, ggrps, symmetry, indices=gindices)
    gtorsion.pot = copy_pot(torsion.pot)

    return gtorsion
