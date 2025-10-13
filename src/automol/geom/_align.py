""" Align two geometries
"""

import operator

import numpy

from ..graph.base import _01networkx
from ._1conv import graph as graph_conv
from .base import distance_matrix, is_atom, reorder


def align(geo1, geo2):
    """ For a given species, align geo1 to geo2 by reordering the
        the atom ordering of geo1 so that it matches that of geo2.
    """

    if not is_atom(geo1):
        # Convert the two geoms to graphs (need explicit graph)
        gra1, gra2 = graph_conv(geo1), graph_conv(geo2)

        # Generate an isomorphism so that we can map gra1 and gra2
        # Will only be able to map most heavy atoms.
        # Can't figure out what to do when heavy atom has two similar R groups
        # Use a distance matrix evaluation to handle that
        igr1 = _01networkx.from_graph(gra1)
        igr2 = _01networkx.from_graph(gra2)
        iso_dcts = _01networkx.all_isomorphisms(igr1, igr2)

        if any(iso_dcts):
            # Calculate distance matrix of geo2
            geo2_distmat = distance_matrix(geo2)

            # Find isomorphism that minimizes diff in dist mat
            aligned_geo1s, diffs = [], []

            for iso_dct in iso_dcts:
                # Reorder the geo1 using the isomoprhism
                _geo1 = reorder(geo1, iso_dct)
                aligned_geo1s.append(_geo1)

                # Calculate metric for difference between geo1 and geo2
                # Distance Matrices by summing all elements of matrix:
                # |G2DistMat - G2DistMat|
                diffs.append(
                    numpy.sum(
                        numpy.abs(geo2_distmat - distance_matrix(_geo1))))

            # Assigned aligned geo1 to the one with smalled dist matrix diff
            min_idx, _ = min(enumerate(diffs), key=operator.itemgetter(1))
            align_geo1 = aligned_geo1s[min_idx]
        else:
            align_geo1 = None
    else:
        align_geo1 = geo1

    return align_geo1
