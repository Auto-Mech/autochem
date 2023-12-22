"""low-level stereochemistry functions
"""
import numbers
from typing import Any, Callable, Dict, Optional, Tuple, Union

import numpy
from automol import util
from automol.graph.base._00core import (
    AtomKey,
    AtomKeys,
    BondKey,
    CenterKey,
    CenterKeys,
    is_ts_graph,
    local_stereo_priorities,
    stereo_parities,
    substitution_atom_transfers,
    tetrahedral_atom_keys,
    tetrahedral_atoms,
)
from automol.graph.base._03kekule import rigid_planar_bond_keys, rigid_planar_bonds
from automol.graph.base._07geom import geometry_atom_parity, geometry_bond_parity
from automol.util import dict_

AtomNeighborDict = Dict[AtomKey, AtomKeys]
BondNeighborDict = Dict[BondKey, Tuple[AtomKeys, AtomKeys]]
CenterNeighborDict = Dict[CenterKey, Union[AtomKeys, Tuple[AtomKeys, AtomKeys]]]
ParityEvaluator = Callable[
    [Any, CenterKeys, Dict[AtomKey, int], CenterNeighborDict, bool],
    Dict[CenterKey, int],
]


def stereocenter_candidates(
    gra, atom: bool = True, bond: bool = True
) -> CenterNeighborDict:
    """Get keys to stereocenter candidates in the graph

    Stereocenter candidates are atoms and bonds which are potentially stereogenic.
    The only thing left to check are the canonical priorities (atom symmetry classes)
    of their neighbors.

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param atom: Include atom stereocenter candidates? defaults to True
    :type atom: bool, optional
    :param bond: Include bond stereocenter candidates? defaults to True
    :type bond: bool, optional
    :returns: A mapping of stereocenter candidates onto their neighbor keys
    :rtype: Dict[Union[int, frozenset[int]], tuple[tuple, tuple]]
    """
    cand_cent_dct = {}

    # 1. Atom stereocenter candidates: tetrahedral atoms
    # These need to be calculated no matter what
    cand_atm_dct = tetrahedral_atoms(gra, min_ncount=3)
    if atom:
        cand_cent_dct.update(cand_atm_dct)

    if bond:
        # 2. Bond stereocenter candidates: (a.) rigid, planar bonds, (b.) at least one
        # neighbor on each side, (c.) not locked in ring with <8 atoms
        cand_bnd_dct = rigid_planar_bonds(gra, min_ncount=1, min_ring_size=8)

        # 3. For TS graphs, remove redundant bond stereocenter candidates, where both
        # atoms are already atom stereocenter candidates
        cand_atm_keys = set(cand_atm_dct.keys())
        cand_bnd_dct = dict_.filter_by_key(
            cand_bnd_dct, lambda bk: not bk <= cand_atm_keys
        )

        cand_cent_dct.update(cand_bnd_dct)

    return cand_cent_dct


def stereocenter_candidates_grouped_and_sorted(
    cand_cent_dct: CenterNeighborDict, pri_dct: Optional[Dict[AtomKey, int]] = None
) -> Tuple[AtomNeighborDict, BondNeighborDict]:
    """Group stereocenter candidates by type

    :param cand_cent_dct: A mapping of stereocenter candidates onto their neighbor keys
    :type cand_cent_dct: CenterNeighborDict
    :param pri_dct: An optional atom priority mapping for sorting the neighbor keys
    :type pri_dct: Optional[Dict[AtomKey, int]]
    :return: The atom and bond stereocenter candidates as separate dictionaries
    :rtype: Tuple[AtomNeighborDict, BondNeighborDict]
    """
    sort_value = (
        None if pri_dct is None else dict_.sort_value_(pri_dct, missing_val=-numpy.inf)
    )

    def sort_atom_neighbors(nkeys):
        return tuple(sorted(nkeys, key=sort_value))

    def sort_bond_neighbors(bnkeys):
        return tuple(tuple(sorted(nks, key=sort_value)) for nks in bnkeys)

    cand_atm_dct = {}
    cand_bnd_dct = {}

    for key, nkeys in cand_cent_dct.items():
        if isinstance(key, numbers.Number):
            cand_atm_dct[key] = sort_atom_neighbors(nkeys)
        else:
            cand_bnd_dct[key] = sort_bond_neighbors(nkeys)

    return cand_atm_dct, cand_bnd_dct


def stereocenter_candidate_keys(
    gra, atom: bool = True, bond: bool = True
) -> frozenset[Union[int, frozenset[int]]]:
    """Get keys to stereocenter candidates in the graph

    Stereocenter candidates are atoms and bonds which are potentially stereogenic.
    The only thing left to check are the canonical priorities (atom symmetry classes)
    of their neighbors.

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param atom: Include atom stereocenter candidates? defaults to True
    :type atom: bool, optional
    :param bond: Include bond stereocenter candidates? defaults to True
    :type bond: bool, optional
    :returns: The keys of the stereocenter candidates
    :rtype: frozenset[Union[int, frozenset[int]]]
    """
    keys = frozenset()

    # 1. Atom stereocenter candidates: tetrahedral atoms
    atm_keys = tetrahedral_atom_keys(gra)  # These need to be calculated no matter what
    if atom:
        keys |= atm_keys

    if bond:
        # 2. Bond stereocenter candidates: (a.) rigid, planar bonds, (b.) at least one
        # neighbor on each side, (c.) not locked in ring with <8 atoms
        bnd_keys = rigid_planar_bond_keys(gra, min_ncount=1, min_ring_size=8)

        # 3. For TS graphs, remove redundant bond stereocenter candidates, where both
        # atoms are already atom stereocenter candidates
        if is_ts_graph(gra):
            bnd_keys = frozenset({bk for bk in bnd_keys if not bk <= atm_keys})

        keys |= bnd_keys

    return keys


def substitution_reversal_parity_flips(
    tsg,
    keys: AtomKeys,
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
) -> Dict[AtomKey, bool]:
    r"""Identify substitution stereocenters for which the local stereo flips upon TS
    direction reversal


        Case 1: Reversal causes local parity flip (reversal value: `True`):

                  3                     3
                  |                     |
            1-----C  +  2   <=>   1  +  C-----2
                 / \                   / \
                5   4                 5   4

                clockwise             counterclockwise
                ('+')                 ('-')

        Case 2: Reversal leaves local parity unchanged (reversal value: `False`):

                  2                     2
                  |                     |
            1-----C  +  5   <=>   1  +  C-----5
                 / \                   / \
                4   3                 4   3

                clockwise             clockwise
                ('+')                 ('+')

        Given a reversal value, r_flip, the local parity upon reversing TS direction is
        given by the following:

            p_rev = p_forw ^ r_flip

    :param tsg: A TS graph
    :type tsg: automol graph data structure
    :param keys: The stereocenter keys at which the parity is being evaluated
    :type keys: CenterKeys
    :param pri_dct: Canonical priorities for the atoms
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
    :type nkeys_dct: CenterNeighborDict
    :return: A mapping identifying which stereoatoms flip upon reversal
    :rtype: Dict[AtomKey, bool]
    """
    subst_dct = substitution_atom_transfers(tsg)

    subst_flip_dct = {}
    for key in keys:
        if key in subst_dct:
            don_key, acc_key = subst_dct[key]
            nkeys0 = sorted(nkeys_dct[key], key=pri_dct.get)
            nkeys0[nkeys0.index(don_key)] = acc_key
            nkeys = sorted(nkeys0, key=pri_dct.get)

            subst_flip_dct[key] = util.is_even_permutation(nkeys0, nkeys)
    return subst_flip_dct


def parity_evaluator_measure_from_geometry_(
    geo,
    local_stereo: bool = False,
    geo_idx_dct: Optional[Dict[AtomKey, AtomKey]] = None,
) -> ParityEvaluator:
    """Creates a parity evaluator which measures parities from a geometry

    :param geo: A geometry
    :type geo: automol geom data structure
    :param local_stereo: Return local stereo assignments? defaults to False
    :type local_stereo: bool, optional
    :param geo_idx_dct: A mapping of graph keys onto geometry indices
    :type geo_idx_dct: Optional[Dict[AtomKey, AtomKey]]
    :returns: A parity evaluator
    :rtype: ParityEvaluator
    """

    def parity_evaluator(
        gra,
        keys: CenterKeys,
        pri_dct: Dict[AtomKey, int],
        nkeys_dct: CenterNeighborDict,
        is_rev_ts: bool = False,
    ):
        """A parity evaluator which reads parities from a graph

        :param gra: A graph
        :type gra: automol graph data structure
        :param keys: The stereocenter keys at which to evaluate the parity
        :type keys: CenterKeys
        :param pri_dct: A canonical atom priority mapping
        :type pri_dct: Dict[AtomKey, int]
        :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
        :type nkeys_dct: CenterNeighborDict
        :param is_rev_ts: Is this a reversed TS graph?, defaults to False
        :type is_rev_ts: bool, optional
        :return: A parity assignment mapping for the given stereocenter keys
        :rtype: Dict[CenterKey, bool]
        """
        pri_dct = local_stereo_priorities(gra) if local_stereo else pri_dct
        nkeys_dct = dict_.by_key(nkeys_dct, keys)
        atm_nkeys_dct, bnd_nkeys_dct = stereocenter_candidates_grouped_and_sorted(
            nkeys_dct, pri_dct=pri_dct
        )

        atm_par_dct = {
            k: geometry_atom_parity(gra, geo, k, ns, geo_idx_dct=geo_idx_dct)
            for k, ns in atm_nkeys_dct.items()
        }
        bnd_par_dct = {
            k: geometry_bond_parity(gra, geo, sorted(k), ns, geo_idx_dct=geo_idx_dct)
            for k, ns in bnd_nkeys_dct.items()
        }

        par_dct = {**atm_par_dct, **bnd_par_dct}

        if is_rev_ts:
            flip_dct = substitution_reversal_parity_flips(gra, keys, pri_dct, nkeys_dct)
            for key, flip in flip_dct.items():
                par_dct[key] ^= flip

        return par_dct

    return parity_evaluator


def parity_evaluator_read_from_graph(
    gra,
    keys: CenterKeys,
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
    is_rev_ts: bool = False,
) -> Dict[CenterKey, bool]:
    """A parity evaluator which reads parities from a graph

    :param gra: A graph
    :type gra: automol graph data structure
    :param keys: The stereocenter keys at which to evaluate the parity
    :type keys: CenterKeys
    :param pri_dct: Canonical priorities for the atoms
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
    :type nkeys_dct: CenterNeighborDict
    :param is_rev_ts: Is this a reversed TS graph?, defaults to False
    :type is_rev_ts: bool, optional
    :return: A mapping of these stereocenter keys onto their parity assignments
    :rtype: Dict[CenterKey, bool]
    """
    assert nkeys_dct or not nkeys_dct
    assert pri_dct or not pri_dct
    assert is_rev_ts or not is_rev_ts

    return dict_.by_key(stereo_parities(gra), keys)


def parity_evaluator_flip_from_graph(
    gra,
    keys: CenterKeys,
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
    is_rev_ts: bool = False,
) -> Dict[CenterKey, bool]:
    """A parity evaluator which reads parities from a graph and flips them to convert
    canonical to local parities (or vice versa)

    :param gra: A graph
    :type gra: automol graph data structure
    :param keys: The stereocenter keys at which the parity is being evaluated
    :type keys: CenterKeys
    :param pri_dct: Canonical priorities for the atoms
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
    :type nkeys_dct: CenterNeighborDict
    :param is_rev_ts: Is this a reversed TS graph?, defaults to False
    :type is_rev_ts: bool, optional
    :return: A mapping of these stereocenter keys onto their parity assignments
    :rtype: Dict[CenterKey, bool]
    """
    loc_pri_dct = local_stereo_priorities(gra)

    nkeys_dct = dict_.by_key(nkeys_dct, keys)
    can_nkeys_dct, can_bnkeys_dct = stereocenter_candidates_grouped_and_sorted(
        nkeys_dct, pri_dct=pri_dct
    )
    loc_nkeys_dct, loc_bnkeys_dct = stereocenter_candidates_grouped_and_sorted(
        nkeys_dct, pri_dct=loc_pri_dct
    )

    par_dct = dict_.by_key(stereo_parities(gra), keys)

    # Flip atom parities
    for key, can_nkeys in can_nkeys_dct.items():
        loc_nkeys = loc_nkeys_dct[key]

        par_dct[key] ^= util.is_odd_permutation(loc_nkeys, can_nkeys)

    # Flip bond parities
    for bkey, (can_nkeys1, can_nkeys2) in can_bnkeys_dct.items():
        (loc_nkeys1, loc_nkeys2) = loc_bnkeys_dct[bkey]

        # Get the maximum-priority neighbor for each
        can_nmax1, can_nmax2, loc_nmax1, loc_nmax2 = (
            nks[-1] for nks in (can_nkeys1, can_nkeys2, loc_nkeys1, loc_nkeys2)
        )

        par_dct[bkey] ^= (loc_nmax1 != can_nmax1) ^ (loc_nmax2 != can_nmax2)

    # Flip substition stereocenters again, if this is a reverse TS graph
    if is_rev_ts:
        flip_dct = substitution_reversal_parity_flips(gra, keys, pri_dct, nkeys_dct)
        for key, flip in flip_dct.items():
            par_dct[key] ^= flip

    return par_dct
