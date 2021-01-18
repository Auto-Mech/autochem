""" Build the grid for a transition state search
"""

import numpy
from scipy.signal import argrelextrema
import automol
import autofile


# Functions for locating maxima
def find_max_1d(typ, grid, ts_zma, dist_name, scn_save_fs,
                mod_thy_info, constraint_dct):
    """ Find the maxmimum of the grid along one dimension
    """

    # Find the maximum along the scan
    locs_list = []
    locs_lst = []
    enes = []
    for grid_val_i in grid:
        if constraint_dct is None:
            locs_list.append([[dist_name], [grid_val_i]])
        else:
            locs_list.append([constraint_dct, [dist_name], [grid_val_i]])
    for locs in locs_list:
        if scn_save_fs[-1].exists(locs):
            scn_path = scn_save_fs[-1].path(locs)
            sp_save_fs = autofile.fs.single_point(scn_path)
            enes.append(sp_save_fs[-1].file.energy.read(mod_thy_info[1:4]))
            locs_lst.append(locs)
    max_ene = max(enes)
    max_idx = enes.index(max_ene)

    # Build lst of guess zmas
    guess_zmas = []

    # Get zma at maximum
    max_locs = locs_lst[max_idx]
    print(scn_save_fs[-1].path(max_locs))
    max_zma = scn_save_fs[-1].file.zmatrix.read(max_locs)
    guess_zmas.append(max_zma)

    # # Add second guess zma for migrations
    if 'migration' in typ:
        max_grid_val = grid[max_idx]
        mig_zma = automol.zmatrix.set_values(
            ts_zma, {dist_name: max_grid_val})
        guess_zmas.append(mig_zma)

    return guess_zmas


def find_max_2d(grid1, grid2, dist_name, brk_name, scn_save_fs,
                mod_thy_info, constraint_dct):
    """ Find the maxmimum of the grid along two dimensions
    """
    enes_lst = []
    locs_lst_lst = []
    for grid_val_j in grid2:
        locs_list = []
        for grid_val_i in grid1:
            if constraint_dct is None:
                locs_list.append([[dist_name, brk_name],
                                  [grid_val_i, grid_val_j]])
            else:
                locs_list.append([constraint_dct, [dist_name, brk_name],
                                  [grid_val_i, grid_val_j]])
        enes = []
        locs_lst = []
        for locs in locs_list:
            print('locs', locs)
            if scn_save_fs[-1].exists(locs):
                scn_path = scn_save_fs[-1].path(locs)
                sp_save_fs = autofile.fs.single_point(scn_path)
                enes.append(sp_save_fs[-1].file.energy.read(mod_thy_info[1:4]))
                locs_lst.append(locs)
        locs_lst_lst.append(locs_lst)
        if enes:
            enes_lst.append(enes)
        print('enes_lst', enes_lst)
    max_enes = []
    max_locs = []
    for idx_j, enes in enumerate(enes_lst):
        max_ene = -10000.
        max_loc = ''
        for idx_i, ene in enumerate(enes):
            print('ene max_ene', ene, max_ene)
            if ene > max_ene:
                max_ene = ene
                max_loc = locs_lst_lst[idx_j][idx_i]
                print('new max', max_ene, max_loc)
        max_enes.append(max_ene)
        max_locs.append(max_loc)
    print('max enes', max_enes)
    min_ene = 10000.
    locs = []
    for idx_j, ene in enumerate(max_enes):
        print('ene min_ene', ene, min_ene)
        if ene < min_ene:
            min_ene = ene
            locs = max_locs[idx_j]
            print('locs', locs)
    max_locs = locs
    max_ene = min_ene
    print('min max loc', max_ene, max_locs)
    print('min max loc', scn_save_fs[-1].path(max_locs))
    max_zma = scn_save_fs[-1].file.zmatrix.read(max_locs)

    # print('geometry for maximum along scan:', max_zma)
    # print('energy for maximum along scan:', max_ene)

    return max_zma


# Functions to build lists potential sadpts
def vtst_max(grid, dist_name, scn_save_fs,
             mod_thy_info, constraint_dct, ethresh=0.3):
    """ Look along a vtst potential and determine if sadpt there
        (need to make the generic version)
    """

    # Get the locs and energies along the grid
    locs_lst, enes_lst = _grid_vals(
        grid, dist_name, scn_save_fs,
        mod_thy_info, constraint_dct)

    # Locate all potential sadpts
    sadpt_idxs, sadpt_enes = _potential_sadpt(enes_lst, ethresh=ethresh)

    if sadpt_idxs and sadpt_enes:
        # For now, find the greatest max for the saddle point
        max_idx = sadpt_enes.index(max(sadpt_enes))
        sadpt_idx = sadpt_idxs[max_idx][1]

        # Get the locs for the maximum
        sadpt_locs = locs_lst[sadpt_idx]

        # Get the max zma
        sadpt_zma = scn_save_fs[-1].file.zmatrix.read(sadpt_locs)
    else:
        sadpt_zma = None

    return sadpt_zma


def _grid_vals(grid, dist_name, scn_save_fs,
               mod_thy_info, constraint_dct):
    """ efef
    """

    # Initialize the lists
    locs_lst = []
    enes_lst = []

    # Build the lists of all the locs for the grid
    grid_locs = []
    for grid_val_i in grid:
        if constraint_dct is None:
            grid_locs.append([[dist_name], [grid_val_i]])
        else:
            grid_locs.append([constraint_dct, [dist_name], [grid_val_i]])

    # Get the energies along the grid
    for locs in grid_locs:
        if scn_save_fs[-1].exists(locs):
            scn_path = scn_save_fs[-1].path(locs)
            sp_save_fs = autofile.fs.single_point(scn_path)
            enes_lst.append(sp_save_fs[-1].file.energy.read(mod_thy_info[1:4]))
            locs_lst.append(locs)

    return locs_lst, enes_lst


def _potential_sadpt(evals, ethresh=0.3):
    """ Determine points on a 1D-grid that could correspond to
        a saddle point
    """

    # Determine the idxs for all of the local extrema
    loc_max, loc_min = _local_extrema(evals)

    # Find min1-max-min2 triplets for each local maxima
    final_idx = len(evals) - 1
    extrema = _extrema_triplets(loc_max, loc_min, final_idx)

    # Determine which local maxima should be considered for a sadpt search
    sadpt_idxs, sadpt_enes = _potential_sadpt_triplets(
        extrema, evals, ethresh=ethresh)

    return sadpt_idxs, sadpt_enes


def _potential_sadpt_triplets(extrema_trips, evals, ethresh=0.3):
    """ Find triplets to look for sadpts
    """

    sadpt_idxs = tuple()
    sadpt_enes = tuple()
    for trip in extrema_trips:
        minidx1, maxidx, minidx2 = trip
        emin1, emax, emin2 = evals[minidx1], evals[maxidx], evals[minidx2]
        edif1 = abs(emax - emin1)
        edif2 = abs(emax - emin2)
        if edif1 >= ethresh or edif2 >= ethresh:
            sadpt_idxs += (trip,)
            sadpt_enes += (emax,)
            # edifs += ((edif1, edif2),)

    return sadpt_idxs, sadpt_enes


def _extrema_triplets(loc_max, loc_min, final_idx):
    """ Build a tuple of triplets for each local maxima where
        each triplet consists of the idx of the maxima and the
        idxs of the minima (or grid endpoint) the maxima is connected
        to.

        param final_idx: index for the final point on the grid (in 0-index)
    """

    trips = tuple()
    for lmax in loc_max:

        # Find inner min connected to emax
        in_idx = lmax
        while in_idx not in loc_min and in_idx != 0:
            in_idx -= 1

        # Find outer min connected to emax
        out_idx = lmax
        while out_idx not in loc_min and out_idx != final_idx:
            out_idx += 1

        # Add to the triplet lst
        trips += ((in_idx, lmax, out_idx),)

    return trips


def _local_extrema(grid):
    """ Find local min and max on a 1D grid
    """

    loc_max = tuple(argrelextrema(numpy.array(grid), numpy.greater)[0])
    loc_min = tuple(argrelextrema(numpy.array(grid), numpy.less)[0])

    return loc_max, loc_min
