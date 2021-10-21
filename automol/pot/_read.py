""" Build the grid for a transition state search
"""

import numpy
from scipy.signal import argrelextrema
from phydat import phycon


# Functions to build lists potential sadpts
def find_max1d(enes_lst,
               max_type='global',
               ethresh=0.01*phycon.KCAL2EH,
               include_endpts=True):
    """ Find the desired local maximum along a 1D-potential. This maximum
        could correspond to the innmermost local maximum or simply the
        global maximum.
    """

    def _global_maximum(idxs, enes):
        """ Determine the index of the global max of a 1D-potential from a
            list of (idx, ene) pairs of local maxima for that potential.

            :param idxs: triplet sets for each local maximum
                (left endpt, max, right endpt), where each is idx for pot lst
            :type idxs: tuple(tuple(int))
            :param enes: energies for each local maximum
            :type enes: tuple(float)
        """
        idx = enes.index(max(enes))
        max_idx = idxs[idx][1]
        return max_idx

    def _innermost_maximum(idxs):
        """ Determine the index of the innermost max of a 1D-potential from a
            list of (idx, ene) pairs of local maxima for that potential.

            :param idxs: triplet sets for each local maximum
                (left endpt, max, right endpt), where each is idx for pot lst
            :type idxs: tuple(tuple(int))
        """
        # for each local maximum idx_set (endpt, max, endpt):
        # find distance b/w middle of potential and max
        mid_idx = int(len(idxs)/2) - 1  # mid idx in 0-indexing
        dist_from_mid_idx = tuple(abs(idx_set[1]-mid_idx)
                                  for idx_set in idxs)
        # determine which local maximum idx_set is closest to middle of pot
        set_idx_of_low_dist = dist_from_mid_idx.index(min(dist_from_mid_idx))
        # Grab the value for max in (endpt, max, endpt) set, this is the idx
        # of the innermost local maximum on potential (idx for whole pot lst)
        max_idx = idxs[set_idx_of_low_dist][1]
        # print('idx_test')
        # print(idxs)
        # print(dist_from_mid_idx)
        # print(set_idx_of_low_dist)
        # print(max_idx)
        return max_idx

    # Locate all potential sadpts
    sadpt_idxs, sadpt_enes = _potential_sadpt(enes_lst, ethresh=ethresh)

    # Determine which of local maxima located meet desired specs; return this
    if sadpt_idxs and sadpt_enes:
        if max_type == 'global':
            max_idx = _global_maximum(sadpt_idxs, sadpt_enes)
        elif max_type == 'innermost':
            # print('sadpt info')
            # print(sadpt_idxs)
            # print(sadpt_enes)
            max_idx = _innermost_maximum(sadpt_idxs)
    else:
        if include_endpts:
            if enes_lst:
                max_idx = 0 if enes_lst[0] > enes_lst[-1] else -1
            else:
                max_idx = None
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
