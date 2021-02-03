"""
 Build the rotor and torsion objects

 Rotor: (tors_obj_1, tors_obj_2, tors_obj_3)
"""

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
        # self.pot = None
        # self.grid = None

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


# Manipulate the torsion objects
def relabel_for_geometry(torsion):
    """ relabel the torsion objec tto correspond with a geometry converted
        from a z-matrix
    """

    # Build the geom that will be returned
    geo = automol.zmat.geometry(torsion.zma)
    
    # Build a remdummy list that tells how to shift the groups
    remdummy = _remove_dummy_lst(torsion.zma)
    if remdummy is not None:
        torsion = _shift_remove_dummy_atom(torsion, remdummy)

    return geo, torsion


def _shift_remove_dummy_atom(torsion, dummy_key, product=False):
    """ shift the values of the torsion groups
    """

    name = torsion.name
    zma = torsion.zma
    symmetry = torsion.symmetry
    groups = tuple(_shift_vals(grp, remdummy) for grp in torsion.groups)
    axis = _shift_vals(torsion.axis, remdummy)
    indices = _shift_vals(torsion.indices, remdummy)

    torsion = Torsion(zma, name, axis, groups, symmetry, indices=indices)

    return torsion


def _remove_dummy_lst(zma):
    """
    Build a remdummy list that tells how to shift the groups
    """

    dummy_idxs = sorted(automol.zmat.atom_indices(zma, 'X', match=True))

    if dummy_idxs:
        remdummy = [0 for _ in range(automol.zmat.count(zma))]
        for dummy in dummy_idxs:
            for idx, _ in enumerate(remdummy):
                if dummy < idx:
                    remdummy[idx] += 1
    else:
        remdummy = None

    return remdummy


def _shift_vals(vals, remdummy):
    """ Shift the values using the remdummy list
    """
    return tuple(val-remdummy[val-1] for val in vals)
