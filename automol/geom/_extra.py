""" extra high-level geometry library functions
"""
from typing import List, Optional, Tuple

import numpy

from ..graph import base as graph_base
from ._1conv import graph, inchi
from .base import (
    almost_equal_coulomb_spectrum,
    almost_equal_dist_matrix,
    central_angle,
    count,
    dihedral_angle,
    distance_matrix,
)

CHECK_DEFAULT_DCT = {"dist": 3.5e-1, "coulomb": 1.5e-2, "stereo": None, "tors": None}


def _similar_dist(geo, geoi, arg=3e-1):
    """Compare the distance matrices of two geometries."""
    return almost_equal_dist_matrix(geo, geoi, thresh=arg)


def _similar_tors(geo, geoi, arg=None):
    """Compare the torsions of two geometries"""
    _ = arg  # Added just to make wrapper function work
    return are_torsions_same(geo, geoi)


def _similar_stereo(geo, geoi, arg=None):
    """Compare the stereochemistry of two geometries"""
    _ = arg  # Added just to make wrapper function work
    ich = inchi(geo)
    ichi = inchi(geoi)
    return bool(ich == ichi)


def _similar_coulomb(geo, geoi, arg=1e-2):
    """Compare the Coulomb spectrum of geometries."""
    return almost_equal_coulomb_spectrum(geo, geoi, rtol=arg)


CHECK_FXN_DCT = {
    "dist": _similar_dist,
    "tors": _similar_tors,
    "stereo": _similar_stereo,
    "coulomb": _similar_coulomb,
}


def are_torsions_same(
    geo1,
    geo2,
    tol: float = 0.09,
    with_h_rotors: bool = True,
    idxs_lst: Optional[List[Tuple[int, int, int, int]]] = None,
):
    """Compare torsional angle values

    :param geo1: A geometry
    :type geo1: automol geom data structure
    :param geo2: Another geometry
    :type geo2: automol geom data structure
    :param tol: Tolerance for the angle difference
    :type tol: float, optional
    :param with_h_rotors: Include H rotors in the comparison?
    :type with_h_rotors: bool, optional
    :param idxs_lst: Specify the exact coordinates to compare
    :type idxs_lst: Optional[List[Tuple[int, int, int, int]]]
    """
    if idxs_lst is None:
        gra = graph(geo1, stereo=False)
        idxs_lst = graph_base.rotational_coordinates(
            gra, segment=True, with_h_rotors=with_h_rotors
        )

    same_dihed = True
    for idxs in idxs_lst:
        ang1 = dihedral_angle(geo1, *idxs)
        ang2 = dihedral_angle(geo2, *idxs)
        diff = numpy.pi - abs(abs(ang1 - ang2) - numpy.pi)
        if diff > tol:
            same_dihed = False
    return same_dihed


# Checks
def is_unique(geo, geo_lst, check_dct=None):
    """Compare one of many structure features of a geometry to that of
    a list of geometries to see if it is unique.

    order of atoms also impacts the comparison as well
    """

    # Set default check values if none are provided
    if check_dct is None:
        check_dct = CHECK_DEFAULT_DCT

    unique = True
    like_idx = None
    for idx, geoi in enumerate(geo_lst):
        # Perform all of the desired comparison checks for similarity
        sim_chk_results = []
        for key, val in check_dct.items():
            kwargs = {"arg": val} if val is not None else {}
            sim_chk_results.append(CHECK_FXN_DCT[key](geo, geoi, **kwargs))
        # If all checks come back as True, than geoms are the same
        if all(sim_chk_results):
            unique = False
            like_idx = idx
            break

    return unique, like_idx


def hydrogen_bonded_structure(geo, dist_thresh=4.82, angle_thresh=1.92, tsg=None):
    """Compare bond lengths in structure to determine if there
    is a hydrogen bond

    :param geo: geometry object
    :type geo: geo object (tuple of tuples)
    :param grxn: reaction object (geo indexing)
    :type grxn: automol.reac.Reaction object
    :param dist_thresh: cutoff value for hbond length (Bohr)
    :type dist_thresh: float
    :param angle_thresh: cutoff value for hbond angle (Radian)
    :type angle_thresh: float
    :rtype: boolean
    """
    hydrogen_bond = hydrogen_bonded_idxs(geo, dist_thresh, angle_thresh, tsg)
    return hydrogen_bond is not None


def hydrogen_bonded_idxs(geo, dist_thresh=5.3, angle_thresh=1.92, tsg=None):
    """Compare bond lengths in structure to determine if there
    is a hydrogen bond.

    :param geo: geometry object
    :type geo: geo object (tuple of tuples)
    :param tsg: A TS graph (geo indexing)
    :type tsg: automol graph data structure
    :param dist_thresh: cutoff value for hbond length (Bohr)
    :type dist_thresh: float
    :param angle_thresh: cutoff value for hbond angle (Radian)
    :type angle_thresh: float
    :rtype: tuple
    """
    # Initialize the hydrogen bond list to None
    hydrogen_bond = None
    if count(geo) > 1:
        # Get the forming/breaking bond idxs if possible
        if tsg is not None:
            frm_bnd_keys = graph_base.ts.forming_bond_keys(tsg)
            brk_bnd_keys = graph_base.ts.breaking_bond_keys(tsg)
            rxn_keys = set()
            for key in frm_bnd_keys:
                rxn_keys = rxn_keys | key
            for key in brk_bnd_keys:
                rxn_keys = rxn_keys | key
            rxn_h_idxs = tuple(rxn_keys)
        else:
            rxn_h_idxs = ()

        # Get all potential indices for HB interactions
        gra = graph(geo)
        dist_mat = distance_matrix(geo)
        adj_atm_dct = graph_base.atoms_neighbor_atom_keys(gra)
        h_idxs = graph_base.atom_keys(gra, symb="H")
        acceptor_idxs = list(graph_base.radical_atom_keys(gra))
        acceptor_idxs.extend(list(graph_base.atom_keys(gra, symb="O")))
        # Loop over indices, ignoring H-idxs in reacting bonds
        hb_idxs = tuple(idx for idx in h_idxs if idx not in rxn_h_idxs)
        for h_idx in hb_idxs:
            for acceptor_idx in acceptor_idxs:
                donor_idx = list(adj_atm_dct[h_idx])[0]
                if acceptor_idx in adj_atm_dct[donor_idx]:
                    continue
                if dist_mat[h_idx][acceptor_idx] < dist_thresh:
                    ang = central_angle(geo, donor_idx, h_idx, acceptor_idx)
                    if ang > angle_thresh:
                        hydrogen_bond = (
                            donor_idx,
                            h_idx,
                            acceptor_idx,
                        )
                        dist_thresh = dist_mat[h_idx][acceptor_idx]
    return hydrogen_bond
