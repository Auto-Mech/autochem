""" scan info for specific reaction classes
"""

import math

import more_itertools as mit
import numpy

from phydat import bnd, phycon

from .. import zmat
from ..const import ReactionClass
from ..graph import ts
from ..util import dict_
from ._0core import Reaction, class_, ts_graph
from ._1util import (
    elimination_breaking_bond_keys,
    hydrogen_migration_atom_keys,
    hydrogen_migration_might_dissociate,
    insertion_forming_bond_keys,
    ring_forming_scission_chain,
)


# Wrapper function to obtain all of the scan data for a reaction
def build_scan_info(zrxn: Reaction, zma, var=False):
    """Build all of the scan information"""

    # Obtain the reactions scan and constraint coordinates
    scan_names = scan_coordinate_name(zrxn, zma)
    const_names = constraint_coordinate_names(zrxn, zma)

    constraint_dct = zmat.constraint_dict(zma, const_names)

    # Build the grid
    grids = scan_grid(zrxn, zma, var=var)

    # Set the update guess
    update_guess = scan_update_guess(zrxn, var=var)

    return scan_names, constraint_dct, grids, update_guess


# SCAN AND CONSTRAINT COORDINATES #
# Unimolecular reactions
# 1. Hydrogen migrations
def hydrogen_migration_scan_coordinate(rxn: Reaction, zma):
    """Obtain the scan coordinate for a hydrogen migration.

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    (frm_bnd_key,) = ts.forming_bond_keys(ts_graph(rxn))
    scan_name = zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return (scan_name,)


def hydrogen_migration_constraint_coordinates(rxn: Reaction, zma):
    """Obtain the constraint coordinates for a hydrogen migration

    :param rxn: a Reaction object
    :returns: the names of the constraint coordinates in the z-matrix
    :rtype: tuple[str]
    """
    att_key, _, don_key, ngb_key = hydrogen_migration_atom_keys(rxn)
    dist_name = zmat.distance_coordinate_name(zma, att_key, ngb_key)
    diss_keys = hydrogen_migration_might_dissociate(rxn, att_key, ngb_key, don_key)
    if diss_keys:
        diss_name = zmat.distance_coordinate_name(zma, *diss_keys)
        ret = (
            dist_name,
            diss_name,
        )
    else:
        ret = (dist_name,)
    return ret


def hydrogen_migration_grid(zrxn: Reaction, zma, npoints=(18,)):
    """Build forward 1D grid  for addition reaction"""

    # Obtain the reactions scan and constraint coordinates
    (scan_name,) = hydrogen_migration_scan_coordinate(zrxn, zma)

    # Build the scan grid
    npoints1 = npoints[0]
    interval = 0.3 * phycon.ANG2BOHR

    frm_bnd_len = _ts_bnd_len(zma, scan_name)
    rmin1 = 2.0 * phycon.ANG2BOHR
    rmin2 = frm_bnd_len + (0.05 * phycon.ANG2BOHR)
    rmax = frm_bnd_len

    if rmax > rmin1:
        npoints = math.ceil((rmax - rmin1) / interval)
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
def beta_scission_scan_coordinate(rxn: Reaction, zma):
    """Obtain the scan coordinate for a beta scission

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    (brk_bnd_key,) = ts.breaking_bond_keys(ts_graph(rxn))
    scan_name = zmat.distance_coordinate_name(zma, *brk_bnd_key)
    return (scan_name,)


def beta_scission_grid(zrxn: Reaction, zma, npoints=(14,)):
    """Build forward 1D grid for a beta scission reaction"""

    # Obtain the reactions scan and constraint coordinates
    (scan_name,) = beta_scission_scan_coordinate(zrxn, zma)

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
def ring_forming_scission_scan_coordinate(rxn: Reaction, zma):
    """Obtain the scan coordinate for a ring-forming scission

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    (brk_bnd_key,) = ts.breaking_bond_keys(ts_graph(rxn))
    scan_name = zmat.distance_coordinate_name(zma, *brk_bnd_key)
    return (scan_name,)


def ring_forming_scission_constraint_coordinates(rxn: Reaction, zma):
    """Obtain the constraint coordinates for a ring-forming scission

    :param rxn: a Reaction object
    :returns: the names of the constraint coordinates in the z-matrix
    :rtype: str
    """
    chain_keys = ring_forming_scission_chain(rxn)
    ang_keys_lst = sorted(mit.windowed(chain_keys[1:], 3))
    dih_keys_lst = sorted(mit.windowed(chain_keys, 4))
    ang_names = [zmat.central_angle_coordinate_name(zma, *ks) for ks in ang_keys_lst]
    dih_names = [zmat.dihedral_angle_coordinate_name(zma, *ks) for ks in dih_keys_lst]
    const_names = tuple(ang_names + dih_names)
    return const_names


def ring_forming_scission_grid(zrxn: Reaction, zma, npoints=(7,)):
    """Build forward WD grid for a ring forming scission reaction
    # the following allows for a 2-d grid search in the initial ts_search
    # for now try 1-d grid and see if it is effective
    """

    # Obtain the scan coordinate
    (scan_name,) = ring_forming_scission_scan_coordinate(zrxn, zma)

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
def elimination_scan_coordinate(rxn: Reaction, zma):
    """Obtain the scan coordinate for an elimination

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    (frm_bnd_key,) = ts.forming_bond_keys(ts_graph(rxn))

    brk_bnd_key1, _ = elimination_breaking_bond_keys(rxn)

    frm_name = zmat.distance_coordinate_name(zma, *frm_bnd_key)
    brk_name = zmat.distance_coordinate_name(zma, *brk_bnd_key1)

    return (frm_name, brk_name)


# def elimination_grid(zrxn, zma, npoints=(10, 10)):
def elimination_grid(zrxn: Reaction, zma, npoints=(7, 5)):
    """Build forward 2D grid for elimination reaction"""

    # Obtain the scan coordinate
    frm_name, brk_name = elimination_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1, npoints2 = npoints

    frm_bnd_len = _ts_bnd_len(zma, frm_name)
    brk_bnd_len = _ts_bnd_len(zma, brk_name)
    if frm_bnd_len is not None:
        r1min = frm_bnd_len + (0.1 * phycon.ANG2BOHR)
        r1max = frm_bnd_len + (0.6 * phycon.ANG2BOHR)
    else:
        r1min = (0.85 + 0.1) * phycon.ANG2BOHR
        r1max = (0.85 + 0.8) * phycon.ANG2BOHR
    if brk_bnd_len is not None:
        r2min = brk_bnd_len + (0.3 * phycon.ANG2BOHR)
        r2max = brk_bnd_len + (1.2 * phycon.ANG2BOHR)
    else:
        r2min = (1.50 + 0.3) * phycon.ANG2BOHR
        r2max = (1.50 + 1.2) * phycon.ANG2BOHR

    grid1 = numpy.linspace(r1min, r1max, npoints1)
    grid2 = numpy.linspace(r2min, r2max, npoints2)

    return (grid1, grid2)


# Bimolecular reactions
# 1. Hydrogen abstractions
def hydrogen_abstraction_scan_coordinate(rxn: Reaction, zma):
    """Obtain the scan coordinate for a hydrogen abstraction

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    (frm_bnd_key,) = ts.forming_bond_keys(ts_graph(rxn))
    scan_name = zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return (scan_name,)


def hydrogen_abstraction_grid(zrxn: Reaction, zma, npoints=(8,)):
    """Build forward 1D grid for hydrogen abstraction reaction"""

    # Obtain the scan coordinate
    (scan_name,) = hydrogen_abstraction_scan_coordinate(zrxn, zma)

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


def radrad_hydrogen_abstraction_grid(zrxn: Reaction, zma, npoints=(8, 4)):
    """Build forward 1D grid for elimination reaction"""

    # Obtain the scan coordinate
    (scan_name,) = hydrogen_abstraction_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1, npoints2 = npoints

    # Get the first grid from close to the mid
    frm_bnd_len = _ts_bnd_len(zma, scan_name)
    if frm_bnd_len is not None:
        rmin = frm_bnd_len + (0.1 * phycon.ANG2BOHR)
        rmax = frm_bnd_len + (1.0 * phycon.ANG2BOHR)
    else:
        rmin = 0.7 * phycon.ANG2BOHR
        rmax = 2.2 * phycon.ANG2BOHR

    grid1 = numpy.linspace(rmin, rmax, npoints1)
    grid1 = numpy.flip(grid1)

    # Get the outer grid from mid to long-distance
    rend2 = 4.0 * phycon.ANG2BOHR
    grid2 = numpy.linspace(rmax, rend2, npoints2 + 1)
    # grid2 = numpy.delete(grid2, 0)

    return (grid1, grid2)


# 2. Additions
def addition_scan_coordinate(rxn: Reaction, zma):
    """Obtain the scan coordinate for an addition

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    (frm_bnd_key,) = ts.forming_bond_keys(ts_graph(rxn))
    scan_name = zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return (scan_name,)


def addition_grid(zrxn: Reaction, zma, npoints=(14,)):
    """Build forward 1D grid for addition reaction"""

    # Obtain the scan coordinate
    (scan_name,) = addition_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(zma, scan_name)
    if frm_bnd_len is not None:
        rmin = frm_bnd_len + (0.1 * phycon.ANG2BOHR)
        rmax = frm_bnd_len + (1.2 * phycon.ANG2BOHR)
    else:
        rmin = 1.6 * phycon.ANG2BOHR
        rmax = 2.8 * phycon.ANG2BOHR

    grid = _geometric_progression(rmin, rmax, npoints1, gfact=1.1, rstp=0.05)

    return (grid,)


def radrad_addition_grid(zrxn: Reaction, zma, npoints=(8, 4)):
    """Build forward 1D grid for a beta scission reaction"""

    # Obtain the scan coordinate
    (scan_name,) = addition_scan_coordinate(zrxn, zma)

    # Build the grid
    npoints1, npoints2 = npoints

    # Get the first grid from close to the mid
    frm_bnd_len = _ts_bnd_len(zma, scan_name)
    if frm_bnd_len is not None:
        rmin = frm_bnd_len + (0.2 * phycon.ANG2BOHR)
        rmax = frm_bnd_len + (1.2 * phycon.ANG2BOHR)
    else:
        rmin = 1.5 * phycon.ANG2BOHR
        rmax = 2.8 * phycon.ANG2BOHR

    grid1 = numpy.linspace(rmin, rmax, npoints1)
    grid1 = numpy.flip(grid1)
    # grid1 = _geometric_progression(
    #     rmin, rmax, npoints1, gfact=1.1, rstp=0.05)

    # Get the outer grid from mid to long-distance
    # Add extra point since initial will be dropped
    rend2 = 4.00 * phycon.ANG2BOHR
    grid2 = numpy.linspace(rmax, rend2, npoints2 + 1)
    # grid2 = numpy.delete(grid2, 0)

    return (grid1, grid2)


# 3. Insertions
def insertion_scan_coordinate(rxn: Reaction, zma):
    """Obtain the scan coordinate for an insertion

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, _ = insertion_forming_bond_keys(rxn)
    scan_name = zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return (scan_name,)


def insertion_grid(zrxn: Reaction, zma, npoints=(16,)):
    """Build forward 1D grid for insertion reaction"""

    # Obtain the scan coordinate
    (scan_name,) = insertion_scan_coordinate(zrxn, zma)

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
def substitution_scan_coordinate(rxn: Reaction, zma):
    """Obtain the scan coordinate for a substitution

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    (frm_bnd_key,) = ts.forming_bond_keys(ts_graph(rxn))
    scan_name = zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return (scan_name,)


def substitution_grid(zrxn: Reaction, zma, npoints=(14,)):
    """Build forward 1D grid for substitution reaction"""

    # Obtain the scan coordinate
    (scan_name,) = substitution_scan_coordinate(zrxn, zma)

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


# Aux function to return empty tuple when info not require for class
def _return_empty_tuple(*_):
    """Return empty tuple"""
    return ()


# Wrapper functions to handle rxn obj and zma for any reaction class
SCAN_COORD_DCT = {
    # unimolecular
    ReactionClass.HYDROGEN_MIGRATION: hydrogen_migration_scan_coordinate,
    ReactionClass.BETA_SCISSION: beta_scission_scan_coordinate,
    ReactionClass.RING_FORM_SCISSION: ring_forming_scission_scan_coordinate,
    ReactionClass.ELIMINATION: elimination_scan_coordinate,
    # bimolecular
    ReactionClass.HYDROGEN_ABSTRACTION: hydrogen_abstraction_scan_coordinate,
    ReactionClass.ADDITION: addition_scan_coordinate,
    ReactionClass.INSERTION: insertion_scan_coordinate,
    ReactionClass.SUBSTITUTION: substitution_scan_coordinate,
}


def scan_coordinate_name(rxn: Reaction, zma):
    """Obtain the scan coordinates

    :param rxn: a hydrogen migration Reaction object
    """
    return SCAN_COORD_DCT[class_(rxn)](rxn, zma)


CONSTRAINT_COORD_DCT = {
    # unimolecular
    ReactionClass.HYDROGEN_MIGRATION: hydrogen_migration_constraint_coordinates,
    ReactionClass.BETA_SCISSION: _return_empty_tuple,
    ReactionClass.RING_FORM_SCISSION: ring_forming_scission_constraint_coordinates,
    ReactionClass.ELIMINATION: _return_empty_tuple,
    # bimolecular
    ReactionClass.HYDROGEN_ABSTRACTION: _return_empty_tuple,
    ReactionClass.ADDITION: _return_empty_tuple,
    ReactionClass.INSERTION: _return_empty_tuple,
    ReactionClass.SUBSTITUTION: _return_empty_tuple,
}


def constraint_coordinate_names(rxn: Reaction, zma):
    """Obtain the constraint coordinates

    :param rxn: a hydrogen migration Reaction object
    """
    return CONSTRAINT_COORD_DCT[class_(rxn)](rxn, zma)


TIGHT_TS_GRID_DCT = {
    ReactionClass.BETA_SCISSION: beta_scission_grid,
    ReactionClass.ADDITION: addition_grid,
    ReactionClass.HYDROGEN_MIGRATION: hydrogen_migration_grid,
    ReactionClass.ELIMINATION: elimination_grid,
    ReactionClass.RING_FORM_SCISSION: ring_forming_scission_grid,
    ReactionClass.HYDROGEN_ABSTRACTION: hydrogen_abstraction_grid,
    ReactionClass.SUBSTITUTION: substitution_grid,
    ReactionClass.INSERTION: insertion_grid,
}
VAR_TS_GRID_DCT = {
    ReactionClass.ADDITION: radrad_addition_grid,
    ReactionClass.HYDROGEN_ABSTRACTION: radrad_hydrogen_abstraction_grid,
}


def scan_grid(zrxn: Reaction, zma, var=False):
    """Set the grid for a transition state search"""

    if not var:
        grid = TIGHT_TS_GRID_DCT[class_(zrxn)](zrxn, zma)
    else:
        grid = VAR_TS_GRID_DCT[class_(zrxn)](zrxn, zma)

    return grid


# UPDATE GUESS DICTIONARY #
TIGHT_TS_UPDATE_GUESS_DCT = {
    ReactionClass.BETA_SCISSION: False,
    ReactionClass.ADDITION: False,
    ReactionClass.HYDROGEN_MIGRATION: True,
    ReactionClass.ELIMINATION: False,
    ReactionClass.RING_FORM_SCISSION: False,
    ReactionClass.HYDROGEN_ABSTRACTION: False,
    ReactionClass.SUBSTITUTION: False,
    ReactionClass.INSERTION: False,
}
VAR_TS_UPDATE_GUESS_DCT = {
    ReactionClass.ADDITION: True,
    ReactionClass.HYDROGEN_ABSTRACTION: True,
}


def scan_update_guess(zrxn: Reaction, var=False):
    """Set boolean to control whether the initial guess structure for
    optimization updates along scan, i.e., uses optimized geometry
    from previous grid point
    """

    if not var:
        _update = TIGHT_TS_UPDATE_GUESS_DCT[class_(zrxn)]
    else:
        _update = VAR_TS_UPDATE_GUESS_DCT[class_(zrxn)]

    return _update


# Helper functions
def _ts_bnd_len(zma, scan_coord):
    """Obtain the current value of the bond defined by the scam coordinate"""

    symbs = zmat.symbols(zma)
    (dist_coo,) = zmat.coordinates(zma)[scan_coord]
    ts_bnd_symbs = tuple(sorted(map(symbs.__getitem__, dist_coo)))
    ts_bnd_len = dict_.value_by_unordered_key(bnd.LEN_DCT, ts_bnd_symbs)

    return ts_bnd_len


def _geometric_progression(rmin, rmax, npoints, gfact=1.1, rstp=0.05):
    """Build a grid using a geometric progresion"""
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
