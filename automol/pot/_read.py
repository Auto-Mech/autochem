""" Build the grid for a transition state search
"""

import numpy
from scipy.signal import argrelextrema
from phydat import phycon


# Functions to build lists potential sadpts
def find_max1d(enes_lst, ethresh=0.01*phycon.KCAL2EH, include_endpts=True):
    """ Look along a vtst potential and determine if sadpt there
        (need to make the generic version)

        Right now simply takes max of all possible saddple-points
        If no saddlepoint found, take one of the endpoints if requested
    """

    # Locate all potential sadpts
    sadpt_idxs, sadpt_enes = _potential_sadpt(enes_lst, ethresh=ethresh)

    if sadpt_idxs and sadpt_enes:
        # For now, find the greatest max for the saddle point
        sadpt_idx = sadpt_enes.index(max(sadpt_enes))
        max_idx = sadpt_idxs[sadpt_idx][1]
    else:
        if include_endpts:
            max_idx = 0 if enes_lst[0] > enes_lst[-1] else -1
        else:
            max_idx = None

    return max_idx


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


def _potential_sadpt_triplets(extrema_trips, evals,
                              ethresh=0.3*phycon.KCAL2EH):
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
