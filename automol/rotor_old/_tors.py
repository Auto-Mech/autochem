"""
 Build the rotor and torsion objects.

 Rotor: (tors_obj_1, tors_obj_2, tors_obj_3)
"""
import copy

import numpy

from automol import graph, reac, zmat
from automol.rotor_old._util import sort_tors_names


class OldTorsion:
    """Describes a torsion, which one or more make up a rotor"""

    def __init__(self, zma, name, axis, groups, symmetry, indices=None):
        """constructor"""

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
        """Obtain the torsional span"""
        self.span = (2.0 * numpy.pi) / self.symmetry

    def set_indices(self):
        """Build indices for the torsion"""
        self.indices = zmat.coord_idxs(self.zma, self.name)


# Build functions
def torsion_lst(zma):
    """Build a list of torsion objects"""

    # Get the necessary graph and lin keys
    gra = zmat.graph(zma, stereo=True, dummy=True)
    lin_keys = sorted(graph.dummy_parent_dict(gra).values())
    tors_axes = tuple(map(tuple, graph.rotational_bond_keys(gra, lin_keys=lin_keys)))
    tors_names = tuple(zmat.torsion_coordinate_name(zma, *keys) for keys in tors_axes)

    # Build the torsion objects
    name_dct = dict(zip(tors_names, tors_axes))

    # Get the sorted tors names for building the list
    sorted_tors_names = sort_tors_names(tuple(name_dct.keys()))

    tors_obj_lst = ()
    for name in sorted_tors_names:
        # Determine constituent rotor pieces in graph system
        axis = name_dct[name]
        grps = graph.rotational_groups(gra, *axis)
        symm = graph.rotational_symmetry_number(gra, *axis, lin_keys=lin_keys)

        # Build the torsion object and add to the list
        tors_obj_lst += (OldTorsion(zma, name, axis, grps, symm),)

    return tors_obj_lst


def reaction_torsion_lst(zma, zrxn):
    """torsions from zrxn obj"""

    ts_zgra = reac.ts_graph(zrxn)
    zbnd_keys = graph.rotational_bond_keys(ts_zgra)
    tors_axes = tuple(tuple(keys) for keys in zbnd_keys)
    tors_names = tuple(zmat.torsion_coordinate_name(zma, *keys) for keys in tors_axes)

    # Get the sorted tors names for building the list
    name_dct = dict(zip(tors_names, tors_axes))

    tors_obj_lst = ()
    sorted_tors_names = sort_tors_names(tuple(name_dct.keys()))
    for name in sorted_tors_names:
        axis = name_dct[name]
        axis = tuple(sorted(axis))
        grps = graph.rotational_groups(ts_zgra, *axis)
        symm = graph.rotational_symmetry_number(ts_zgra, *axis)

        # Build the torsion object and add to the list
        tors_obj_lst += (OldTorsion(zma, name, axis, grps, symm),)

    return tors_obj_lst


# Manipulate the torsion objects
def relabel_for_geometry(torsion):
    """relabel the torsion objec tto correspond with a geometry converted
    from a z-matrix
    """

    name = torsion.name
    zma = torsion.zma
    symmetry = torsion.symmetry
    ggrps = tuple(zmat.shift_down(zma, grp) for grp in torsion.groups)
    gaxis = zmat.shift_down(zma, torsion.axis)
    gindices = zmat.shift_down(zma, torsion.indices)

    gtorsion = OldTorsion(zma, name, gaxis, ggrps, symmetry, indices=gindices)
    gtorsion.pot = copy.deepcopy(torsion.pot)

    return gtorsion
