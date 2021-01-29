""" Build the grid for a transition state search
"""

import math
import numpy
import automol
from automol.par import ReactionClass
from phydat import phycon
from phydat import bnd


# Functions to build the grid for searching for TSs
def build_grid(rclass, ts_zma, trans):
    """ Set the grid for a transition state search

        rclass = (typ, spin, radrad)
        # Pass npoints as a 2-element list
    """

    # Unpack the rclass
    rtyp, spin, radrad = rclass

    # Set the main type
    if radrad and spin == 'low':
        grid, update_guess = VAR_TS_GRID_BUILDER_DCT[rtyp](
            ts_zma, trans)
    else:
        grid, update_guess = TIGHT_TS_GRID_BUILDER_DCT[rtyp](
            ts_zma, trans)

    return grid, update_guess


# Tight TS grid
def ring_forming_scission_grid(ts_zma, trans,
                               npoints=(7,), update_guess=False):
    """ Build forward WD grid for a ring forming scission reaction
        # the following allows for a 2-d grid search in the initial ts_search
        # for now try 1-d grid and see if it is effective
    """

    npoints1 = npoints[0]

    brk_bnd_len = _ts_bnd_len(ts_zma, trans, choice='brk')
    if brk_bnd_len is not None:
        r1min = brk_bnd_len + 0.1 * phycon.ANG2BOHR
        r1max = brk_bnd_len + 0.7 * phycon.ANG2BOHR
    else:
        r1min = (1.54 + 0.1) * phycon.ANG2BOHR
        r1max = (1.54 + 0.7) * phycon.ANG2BOHR

    grid1 = numpy.linspace(r1min, r1max, npoints1)
    grid = grid1

    return grid, update_guess


def beta_scission_grid(ts_zma, trans,
                       npoints=(14,), update_guess=False):
    """ Build forward 1D grid for a beta scission reaction
    """

    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(ts_zma, trans, choice='frm')
    if frm_bnd_len is not None:
        rmin = frm_bnd_len + 0.1 * phycon.ANG2BOHR
        rmax = frm_bnd_len + 0.8 * phycon.ANG2BOHR
    else:
        rmin = 1.4 * phycon.ANG2BOHR
        rmax = 2.0 * phycon.ANG2BOHR
    grid = numpy.linspace(rmin, rmax, npoints1)

    return grid, update_guess


def addition_grid(ts_zma, trans,
                  npoints=(14,), update_guess=False):
    """ Build forward 1D grid for addition reaction
    """

    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(ts_zma, trans, choice='frm')
    if frm_bnd_len is not None:
        rmin = frm_bnd_len + 0.1 * phycon.ANG2BOHR
        rmax = frm_bnd_len + 1.2 * phycon.ANG2BOHR
    else:
        rmin = 1.6 * phycon.ANG2BOHR
        rmax = 2.8 * phycon.ANG2BOHR

    grid = _geometric_progression(
        rmin, rmax, npoints1, gfact=1.1, rstp=0.05)

    return grid, update_guess


def hydrogen_migration_grid(ts_zma, trans,
                            npoints=(18,), update_guess=True):
    """ Build forward 1D grid  for addition reaction
    """

    interval = 0.3*phycon.ANG2BOHR

    frm_bnd_len = _ts_bnd_len(ts_zma, trans, choice='frm')
    rmin1 = 2.0*phycon.ANG2BOHR
    rmin2 = frm_bnd_len + 0.05 * phycon.ANG2BOHR
    rmax = frm_bnd_len

    if rmax > rmin1:
        npoints = math.ceil((rmax-rmin1)/interval)
        if npoints < 1:
            grid1 = []
        else:
            grid1 = numpy.linspace(rmax, rmin1, npoints)
    else:
        grid1 = []

    grid2 = numpy.linspace(rmin1, rmin2, 18)
    grid = numpy.concatenate((grid1, grid2), axis=None)

    return grid, update_guess


def unimolecular_elimination_grid(ts_zma, trans,
                                  npoints=(8, 4), update_guess=False):
    """ Build forward 2D grid for elimination reaction
    """

    npoints1, npoints2 = npoints

    frm_bnd_len = _ts_bnd_len(ts_zma, trans, choice='frm')
    brk_bnd_len = _ts_bnd_len(ts_zma, trans, choice='brk')
    if frm_bnd_len is not None and brk_bnd_len is not None:
        r1min = frm_bnd_len + 0.2
        r1max = frm_bnd_len + 1.4
        r2min = brk_bnd_len + 0.2
        r2max = brk_bnd_len + 0.8
    else:
        r1min = (1.54 + 0.2) * phycon.ANG2BOHR
        r1max = (1.54 + 1.4) * phycon.ANG2BOHR
        r2min = (0.74 + 0.2) * phycon.ANG2BOHR
        r2max = (0.74 + 0.8) * phycon.ANG2BOHR

    grid1 = numpy.linspace(r1min, r1max, npoints1) * phycon.ANG2BOHR
    grid2 = numpy.linspace(r2min, r2max, npoints2) * phycon.ANG2BOHR
    grid = (grid1, grid2)

    return grid, update_guess


def hydrogen_abstraction(ts_zma, trans,
                         npoints=(8,), update_guess=False):
    """ Build forward 1D grid for hydrogen abstraction reaction
    """

    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(ts_zma, trans, choice='frm')
    if frm_bnd_len is not None:
        rmin = frm_bnd_len + 0.2
        rmax = frm_bnd_len + 1.0 * phycon.ANG2BOHR
    else:
        rmin = 0.7 * phycon.ANG2BOHR
        rmax = 2.2 * phycon.ANG2BOHR
    grid = numpy.linspace(rmin, rmax, npoints1)

    return grid, update_guess


def substitution(ts_zma, trans,
                 npoints=(14,), update_guess=False):
    """ Build forward 1D grid for substitution reaction
    """

    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(ts_zma, trans, choice='frm')
    if frm_bnd_len is not None:
        rmin = frm_bnd_len
        rmax = frm_bnd_len + 1.4 * phycon.ANG2BOHR
    else:
        rmin = 0.7 * phycon.ANG2BOHR
        rmax = 2.4 * phycon.ANG2BOHR

    grid = numpy.linspace(rmin, rmax, npoints1)

    return grid, update_guess


def insertion(ts_zma, trans,
              npoints=(16,), update_guess=False):
    """ Build forward 1D grid for insertion reaction
    """

    npoints1 = npoints[0]

    frm_bnd_len = _ts_bnd_len(ts_zma, trans, choice='frm')
    if frm_bnd_len is not None:
        rmin = frm_bnd_len
        rmax = frm_bnd_len + 1.4 * phycon.ANG2BOHR
    else:
        rmin = 1.4 * phycon.ANG2BOHR
        rmax = 2.4 * phycon.ANG2BOHR

    grid = numpy.linspace(rmin, rmax, npoints1)

    return grid, update_guess


# Barrierless TS grid
def radrad_addition_grid(npoints=(5, 6), update_guess=True):
    """ Build forward 1D grid for a beta scission reaction
    """

    npoints1, npoints2 = npoints

    rstart = 2.6 * phycon.ANG2BOHR
    rend1 = 1.8 * phycon.ANG2BOHR
    rend2 = 3.85 * phycon.ANG2BOHR

    grid1 = numpy.linspace(rstart, rend1, npoints1)
    grid2 = numpy.linspace(rstart, rend2, npoints2)
    # grid2 = numpy.delete(grid2, 0)
    grid = [grid1, grid2]
    # grid = numpy.concatenate((grid1, grid2), axis=None)

    return grid, update_guess


def radrad_hydrogen_abstraction_grid(npoints=(8, 4), update_guess=True):
    """ Build forward 1D grid for elimination reaction
    """

    npoints1, npoints2 = 8, 4

    rstart = 2.4 * phycon.ANG2BOHR
    rend1 = 1.4 * phycon.ANG2BOHR
    rend2 = 3.0 * phycon.ANG2BOHR

    grid1 = numpy.linspace(rstart, rend1, npoints1)
    grid2 = numpy.linspace(rstart, rend2, npoints2)
    grid2 = numpy.delete(grid2, 0)
    # grid = numpy.concatenate((grid1, grid2), axis=None)
    grid = [grid1, grid2]

    return grid, update_guess


TIGHT_TS_GRID_BUILDER_DCT = {
    ReactionClass.BETA_SCISSION: beta_scission_grid,
    ReactionClass.ADDITION: addition_grid,
    ReactionClass.HYDROGEN_MIGRATION: hydrogen_migration_grid,
    ReactionClass.ELIMINATION: unimolecular_elimination_grid,
    ReactionClass.RING_FORM_SCISSION: ring_forming_scission_grid,
    ReactionClass.HYDROGEN_ABSTRACTION: hydrogen_abstraction,
    ReactionClass.SUBSTITUTION: substitution,
    ReactionClass.INSERTION: insertion,
}

VAR_TS_GRID_BUILDER_DCT = {
    ReactionClass.ADDITION: radrad_addition_grid,
    ReactionClass.HYDROGEN_ABSTRACTION: radrad_hydrogen_abstraction_grid
}


# Get the ts bond length
def _ts_bnd_len(ts_zma, trans, choice):
    """ get the length of a forming or a breaking bond
    """

    assert choice in ('frm', 'brk')
    if choice == 'frm':
        ts_keys = automol.graph.trans.formed_bond_keys(trans)
    else:
        ts_keys = automol.graph.trans.broken_bond_keys(trans)

    # Get the symbols of the atoms in the TS bond
    ts_bnd_name = automol.zmat.bond_key_from_idxs(ts_zma, ts_keys)
    dist_coo, = automol.zmat.coordinates(ts_zma)[ts_bnd_name]
    symbs = automol.zmat.symbols(ts_zma)
    ts_bnd_symbs = tuple(sorted(map(symbs.__getitem__, dist_coo)))

    # Obtain the standard X-Y bond length for atoms X,Y in the TS Bond
    ts_bnd_len = bnd.read_len(ts_bnd_symbs)

    return ts_bnd_len


# Special grid progressions
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
