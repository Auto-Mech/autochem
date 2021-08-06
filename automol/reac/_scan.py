""" scan info for specific reaction classes
"""

import math
import numpy
import more_itertools as mit
from phydat import phycon, bnd
from automol.graph import ts
from automol.par import ReactionClass
import automol.zmat
from automol.util import dict_
# from automol.util import numpy_to_float
from automol.reac._util import hydrogen_migration_atom_keys
from automol.reac._util import ring_forming_scission_chain
from automol.reac._util import insertion_forming_bond_keys


# Wrapper function to obtain all of the scan data for a reaction
def build_scan_info(zrxn, zma):
    """ Build all of the scan information
    """

    # Obtain the reactions scan and constraint coordinates
    scan_names = scan_coordinate(zrxn, zma)
    const_names = constraint_coordinates(zrxn, zma)

    constraint_dct = automol.zmat.constraint_dct(zma, const_names)

    # Build the grid
    grids = scan_grid(zrxn, zma)

    # Set the update guess
    update_guess = UPDATE_GUESS_DCT[zrxn.class_]

    return scan_names, constraint_dct, grids, update_guess


# SCAN AND CONSTRAINT COORDINATES #
# Unimolecular reactions
# 1. Hydrogen migrations
def hydrogen_migration_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for a hydrogen migration.

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return (scan_name,)


def hydrogen_migration_constraint_coordinates(rxn, zma):
    """ Obtain the constraint coordinates for a hydrogen migration

    :param rxn: a Reaction object
    :returns: the names of the constraint coordinates in the z-matrix
    :rtype: tuple[str]
    """
    att_key, _, _, ngb_key = hydrogen_migration_atom_keys(rxn)
    dist_name = automol.zmat.distance_coordinate_name(zma, att_key, ngb_key)
    return (dist_name,)


def hydrogen_migration_grid(zrxn, zma, npoints=(18,)):
    """ Build forward 1D grid  for addition reaction
    """

    # Obtain the reactions scan and constraint coordinates
    scan_name, = hydrogen_migration_scan_coordinate(zrxn, zma)
    
    # Build the scan grid
    npoints1 = npoints[0]
    interval = 0.3 * phycon.ANG2BOHR

    frm_bnd_len = _ts_bnd_len(zma, scan_name)
    rmin1 = 2.0 * phycon.ANG2BOHR
    rmin2 = frm_bnd_len + (0.05 * phycon.ANG2BOHR)
    rmax = frm_bnd_len

    if rmax > rmin1:
        npoints = math.ceil((rmax-rmin1)/interval)
        if npoints < 1:
            grid1 = []
        else:
            grid1 = numpy.linspace(rmax, rmin1, npoints)
    else:
        grid1 = []

    grid2 = numpy.linspace(rmin1, rmin2, npoints1)
    grid = numpy.concatenate((grid1, grid2), axis=None)

    return (grid,)


# 2. Beta scissions
def beta_scission_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for a beta scission

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *brk_bnd_key)
    return (scan_name,)


def beta_scission_grid(zrxn, zma, npoints=(14,)):
    """ Build forward 1D grid for a beta scission reaction
    """

    # Obtain the reactions scan and constraint coordinates
    scan_name, = beta_scission_scan_coordinate(zrxn, zma)

    # Build the scan grid
    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(zma, scan_name)
    if frm_bnd_len is not None:
        rmin = frm_bnd_len + (0.1 * phycon.ANG2BOHR)
        rmax = frm_bnd_len + (0.8 * phycon.ANG2BOHR)
    else:
        rmin = 1.4 * phycon.ANG2BOHR
        rmax = 2.0 * phycon.ANG2BOHR
    grid = numpy.linspace(rmin, rmax, npoints1)

    return (grid,)


# 3. Ring-forming scissions
def ring_forming_scission_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for a ring-forming scission

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *brk_bnd_key)
    return (scan_name,)


def ring_forming_scission_constraint_coordinates(rxn, zma):
    """ Obtain the constraint coordinates for a ring-forming scission

    :param rxn: a Reaction object
    :returns: the names of the constraint coordinates in the z-matrix
    :rtype: str
    """
    chain_keys = ring_forming_scission_chain(rxn)
    ang_keys_lst = sorted(mit.windowed(chain_keys[1:], 3))
    dih_keys_lst = sorted(mit.windowed(chain_keys, 4))
    ang_names = [automol.zmat.central_angle_coordinate_name(zma, *ks)
                 for ks in ang_keys_lst]
    dih_names = [automol.zmat.dihedral_angle_coordinate_name(zma, *ks)
                 for ks in dih_keys_lst]
    const_names = tuple(ang_names + dih_names)
    return const_names


def ring_forming_scission_grid(zrxn, zma, npoints=(7,)):
    """ Build forward WD grid for a ring forming scission reaction
        # the following allows for a 2-d grid search in the initial ts_search
        # for now try 1-d grid and see if it is effective
    """

    # Obtain the scan coordinate
    scan_name, = ring_forming_scission_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1 = npoints[0]

    brk_bnd_len = _ts_bnd_len(zma, scan_name)
    if brk_bnd_len is not None:
        r1min = brk_bnd_len + (0.1 * phycon.ANG2BOHR)
        r1max = brk_bnd_len + (0.7 * phycon.ANG2BOHR)
    else:
        r1min = (1.54 + 0.1) * phycon.ANG2BOHR
        r1max = (1.54 + 0.7) * phycon.ANG2BOHR

    grid1 = numpy.linspace(r1min, r1max, npoints1)
    grid = tuple(val.item() for val in grid1)

    return (grid,)


# 4. Eliminations
def elimination_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for an elimination

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    brk_bnd_keys = ts.breaking_bond_keys(rxn.forward_ts_graph)
    brk_bnd_key1, brk_bnd_key2 = brk_bnd_keys

    # Two bonds breaking in eliminations, need to choose brk bnd
    if len(frm_bnd_key | brk_bnd_key1 | brk_bnd_key2) > 3:
        # for ring_size > 3: use brk-bnd that doesn't involve atoms in frm bond
        scn_brk_bnd_key = None
        for brk_bnd_key in brk_bnd_keys:
            if not frm_bnd_key & brk_bnd_key:
                scn_brk_bnd_key = brk_bnd_key
    else:
        # if one brk bnd doesn't have H, use that, else use arbitrary brk bnd
        scn_brk_bnd_key = None
        symbs = automol.zmat.symbols(zma)
        for brk_bnd_key in brk_bnd_keys:
            brk_symbs = tuple(symbs[key] for key in brk_bnd_key1)
            if 'H' not in brk_symbs:
                scn_brk_bnd_key = brk_bnd_key
                break
        if scn_brk_bnd_key is None:
            scn_brk_bnd_key = brk_bnd_key1

    frm_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    brk_name = automol.zmat.distance_coordinate_name(zma, *scn_brk_bnd_key)

    return (frm_name, brk_name)


def elimination_grid(zrxn, zma, npoints=(8, 4)):
    """ Build forward 2D grid for elimination reaction
    """

    # Obtain the scan coordinate
    frm_name, brk_name = elimination_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1, npoints2 = npoints

    frm_bnd_len = _ts_bnd_len(zma, frm_name)
    brk_bnd_len = _ts_bnd_len(zma, brk_name)
    print('frm len', frm_name, frm_bnd_len)
    print('brk len', brk_name, brk_bnd_len)
    if frm_bnd_len is not None:
        r1min = frm_bnd_len + (0.2 * phycon.ANG2BOHR)
        r1max = frm_bnd_len + (1.4 * phycon.ANG2BOHR)
    else:
        r1min = (1.54 + 0.2) * phycon.ANG2BOHR
        r1max = (1.54 + 1.4) * phycon.ANG2BOHR
    if brk_bnd_len is not None:
        r2min = brk_bnd_len + (0.2 * phycon.ANG2BOHR)
        r2max = brk_bnd_len + (0.8 * phycon.ANG2BOHR)
    else:
        r2min = (0.74 + 0.2) * phycon.ANG2BOHR
        r2max = (0.74 + 0.8) * phycon.ANG2BOHR

    print('rmax', r1min, r1max, r2min, r2max)

    grid1 = numpy.linspace(r1min, r1max, npoints1)
    grid2 = numpy.linspace(r2min, r2max, npoints2)

    return (grid1, grid2)


# Bimolecular reactions
# 1. Hydrogen abstractions
def hydrogen_abstraction_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for a hydrogen abstraction

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return (scan_name,)


def hydrogen_abstraction_grid(zrxn, zma, npoints=(8,)):
    """ Build forward 1D grid for hydrogen abstraction reaction
    """

    # Obtain the scan coordinate
    scan_name, = hydrogen_abstraction_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(zma, scan_name)
    if frm_bnd_len is not None:
        rmin = frm_bnd_len + (0.1 * phycon.ANG2BOHR)
        rmax = frm_bnd_len + (1.0 * phycon.ANG2BOHR)
    else:
        rmin = 0.7 * phycon.ANG2BOHR
        rmax = 2.2 * phycon.ANG2BOHR
    grid = numpy.linspace(rmin, rmax, npoints1)

    return (grid,)


def radrad_hydrogen_abstraction_grid(npoints=(8, 4)):
    """ Build forward 1D grid for elimination reaction
    """

    # Build the grid
    npoints1, npoints2 = npoints

    rstart = 2.4 * phycon.ANG2BOHR
    rend1 = 1.4 * phycon.ANG2BOHR
    rend2 = 3.0 * phycon.ANG2BOHR

    grid1 = numpy.linspace(rstart, rend1, npoints1)
    grid2 = numpy.linspace(rstart, rend2, npoints2)
    grid2 = numpy.delete(grid2, 0)

    return (grid1, grid2)


# 2. Additions
def addition_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for an addition

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return (scan_name,)


def addition_grid(zrxn, zma, npoints=(14,)):
    """ Build forward 1D grid for addition reaction
    """

    # Obtain the scan coordinate
    scan_name, = addition_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(zma, scan_name)
    if frm_bnd_len is not None:
        rmin = frm_bnd_len + (0.1 * phycon.ANG2BOHR)
        rmax = frm_bnd_len + (1.2 * phycon.ANG2BOHR)
    else:
        rmin = 1.6 * phycon.ANG2BOHR
        rmax = 2.8 * phycon.ANG2BOHR

    grid = _geometric_progression(
        rmin, rmax, npoints1, gfact=1.1, rstp=0.05)

    return (grid,)


def radrad_addition_grid(npoints=(5, 6)):
    """ Build forward 1D grid for a beta scission reaction
    """

    # Obtain the scan coordinate
    # scan_name = addition_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1, npoints2 = npoints

    rstart = 2.6 * phycon.ANG2BOHR
    rend1 = 1.8 * phycon.ANG2BOHR
    rend2 = 3.85 * phycon.ANG2BOHR

    grid1 = numpy.linspace(rstart, rend1, npoints1)
    grid2 = numpy.linspace(rstart, rend2, npoints2)
    # grid2 = numpy.delete(grid2, 0)

    return (grid1, grid2)


# 3. Insertions
def insertion_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for an insertion

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, _ = insertion_forming_bond_keys(rxn)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return (scan_name,)


def insertion_grid(zrxn, zma, npoints=(16,)):
    """ Build forward 1D grid for insertion reaction
    """

    # Obtain the scan coordinate
    scan_name, = insertion_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(zma, scan_name)
    if frm_bnd_len is not None:
        rmin = frm_bnd_len
        rmax = frm_bnd_len + (1.4 * phycon.ANG2BOHR)
    else:
        rmin = 1.4 * phycon.ANG2BOHR
        rmax = 2.4 * phycon.ANG2BOHR

    grid = numpy.linspace(rmin, rmax, npoints1)

    return (grid,)


# 4. Substitution
def substitution_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for a substitution

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return (scan_name,)


def substitution_grid(zrxn, zma, npoints=(14,)):
    """ Build forward 1D grid for substitution reaction
    """

    # Obtain the scan coordinate
    scan_name, = substitution_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(zma, scan_name)
    if frm_bnd_len is not None:
        rmin = frm_bnd_len
        rmax = frm_bnd_len + (1.4 * phycon.ANG2BOHR)
    else:
        rmin = 0.7 * phycon.ANG2BOHR
        rmax = 2.4 * phycon.ANG2BOHR

    grid = numpy.linspace(rmin, rmax, npoints1)

    return (grid,)


# Wrapper functions to handle rxn obj and zma for any reaction class
def scan_coordinate(rxn, zma):
    """ Obtain the scan coordinates

    :param rxn: a hydrogen migration Reaction object
    """
    function_dct = {
        # unimolecular
        ReactionClass.Typ.HYDROGEN_MIGRATION:
        hydrogen_migration_scan_coordinate,
        ReactionClass.Typ.BETA_SCISSION: beta_scission_scan_coordinate,
        ReactionClass.Typ.RING_FORM_SCISSION:
        ring_forming_scission_scan_coordinate,
        ReactionClass.Typ.ELIMINATION: elimination_scan_coordinate,
        # bimolecular
        ReactionClass.Typ.HYDROGEN_ABSTRACTION:
        hydrogen_abstraction_scan_coordinate,
        ReactionClass.Typ.ADDITION: addition_scan_coordinate,
        ReactionClass.Typ.INSERTION: insertion_scan_coordinate,
        ReactionClass.Typ.SUBSTITUTION: substitution_scan_coordinate,
    }

    fun_ = function_dct[rxn.class_]
    ret = fun_(rxn, zma)
    return ret


def constraint_coordinates(rxn, zma):
    """ Obtain the constraint coordinates

    :param rxn: a hydrogen migration Reaction object
    """

    def _return_empty_tuple(*_):
        return ()

    function_dct = {
        # unimolecular
        ReactionClass.Typ.HYDROGEN_MIGRATION:
        hydrogen_migration_constraint_coordinates,
        ReactionClass.Typ.BETA_SCISSION: _return_empty_tuple,
        ReactionClass.Typ.RING_FORM_SCISSION:
        ring_forming_scission_constraint_coordinates,
        ReactionClass.Typ.ELIMINATION: _return_empty_tuple,
        # bimolecular
        ReactionClass.Typ.HYDROGEN_ABSTRACTION:
        _return_empty_tuple,
        ReactionClass.Typ.ADDITION: _return_empty_tuple,
        ReactionClass.Typ.INSERTION: _return_empty_tuple,
        ReactionClass.Typ.SUBSTITUTION: _return_empty_tuple,
    }

    fun_ = function_dct[rxn.class_]
    ret = fun_(rxn, zma)
    return ret


def scan_grid(zrxn, zma):
    """ Set the grid for a transition state search

        rclass = (typ, spin, radrad)
        # Pass npoints as a 2-element list
    """

    tight_ts_grid_builder_dct = {
        ReactionClass.Typ.BETA_SCISSION: beta_scission_grid,
        ReactionClass.Typ.ADDITION: addition_grid,
        ReactionClass.Typ.HYDROGEN_MIGRATION: hydrogen_migration_grid,
        ReactionClass.Typ.ELIMINATION: elimination_grid,
        ReactionClass.Typ.RING_FORM_SCISSION: ring_forming_scission_grid,
        ReactionClass.Typ.HYDROGEN_ABSTRACTION: hydrogen_abstraction_grid,
        ReactionClass.Typ.SUBSTITUTION: substitution_grid,
        ReactionClass.Typ.INSERTION: insertion_grid
    }

    # var_ts_grid_builder_dct = {
    #  ReactionClass.Typ.ADDITION: radrad_addition_grid,
    #  ReactionClass.Typ.HYDROGEN_ABSTRACTION: radrad_hydrogen_abstraction_grid
    # }

    grid = tight_ts_grid_builder_dct[zrxn.class_](zrxn, zma)

    # Set the main type
    # if radrad and spin == 'low':
    #     grid, update_guess = VAR_TS_GRID_BUILDER_DCT[rtyp](
    #         ts_zma, trans)
    # else:
    #     grid, update_guess = TIGHT_TS_GRID_BUILDER_DCT[rtyp](
    #         ts_zma, trans)

    return grid


# UPDATE GUESS DICTIONARY #
UPDATE_GUESS_DCT = {
    ReactionClass.Typ.BETA_SCISSION: False,
    ReactionClass.Typ.ADDITION: False,
    ReactionClass.Typ.HYDROGEN_MIGRATION: True,
    ReactionClass.Typ.ELIMINATION: False,
    ReactionClass.Typ.RING_FORM_SCISSION: False,
    ReactionClass.Typ.HYDROGEN_ABSTRACTION: False,
    ReactionClass.Typ.SUBSTITUTION: False,
    ReactionClass.Typ.INSERTION: False
}
# ReactionClass.Typ.ADDITION: True,
# ReactionClass.Typ.HYDROGEN_ABSTRACTION: True
# These should be for radical-radical cases of the two - add later


# Helper functions
def _ts_bnd_len(zma, scan_coord):
    """ Obtain the current value of the bond defined by the scam coordinate
    """

    symbs = automol.zmat.symbols(zma)
    dist_coo, = automol.zmat.coordinates(zma)[scan_coord]
    ts_bnd_symbs = tuple(sorted(map(symbs.__getitem__, dist_coo)))
    ts_bnd_len = dict_.values_by_unordered_tuple(bnd.LEN_DCT, ts_bnd_symbs)

    return ts_bnd_len


def _geometric_progression(rmin, rmax, npoints, gfact=1.1, rstp=0.05):
    """ Build a grid using a geometric progresion
    """
    grid = [rmin]
    rgrid = rmin
    for _ in range(npoints):
        rgrid += rstp
        if rgrid == rmax:
            break
        grid.append(rgrid)
        rstp = rstp * gfact
    grid = numpy.array(grid)

    return grid
