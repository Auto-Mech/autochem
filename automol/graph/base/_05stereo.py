"""low-level stereochemistry functions
"""
import numbers
from collections import abc
from typing import Any, Callable, Dict, Optional, Tuple, Union

import numpy
from automol import util
from automol.geom import base as geom_base
from automol.graph.base._00core import (
    AtomKey,
    AtomKeys,
    BondKey,
    CenterKey,
    CenterKeys,
    atom_keys,
    explicit,
    is_ts_graph,
    local_stereo_priorities,
    stereo_parities,
    tetrahedral_atom_keys,
    tetrahedral_atoms,
)
from automol.graph.base._03kekule import (
    rigid_planar_bond_keys,
    rigid_planar_bonds,
    vinyl_radical_atom_bond_keys,
)
from automol.graph.base._04class import (
    atom_stereo_sorted_neighbor_keys,
    bond_stereo_sorted_neighbor_keys,
    substitution_atom_transfers,
)
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


def local_parity_flips(keys, pri_dct, loc_pri_dct, nkeys_dct):
    """Parity flips for converting canonical to local parities (or vice versa)

    :param keys: The stereocenter keys at which the parity is being evaluated
    :type keys: CenterKeys
    :param pri_dct: The canonical priority mapping
    :type pri_dct: Dict[AtomKey, int]
    :param loc_pri_dct: The local priority mapping
    :type loc_pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
    :type nkeys_dct: CenterNeighborDict
    """
    nkeys_dct = dict_.by_key(nkeys_dct, keys)
    can_nkeys_dct, can_bnkeys_dct = stereocenter_candidates_grouped_and_sorted(
        nkeys_dct, pri_dct=pri_dct
    )
    loc_nkeys_dct, loc_bnkeys_dct = stereocenter_candidates_grouped_and_sorted(
        nkeys_dct, pri_dct=loc_pri_dct
    )

    flip_dct = {}

    # Determine atom parity flips
    for key, can_nkeys in can_nkeys_dct.items():
        loc_nkeys = loc_nkeys_dct[key]

        flip_dct[key] = util.is_odd_permutation(loc_nkeys, can_nkeys)

    # Determine bond parity flips
    for bkey, (can_nkeys1, can_nkeys2) in can_bnkeys_dct.items():
        (loc_nkeys1, loc_nkeys2) = loc_bnkeys_dct[bkey]

        # Get the maximum-priority neighbor for each atom
        can_nmax1, can_nmax2, loc_nmax1, loc_nmax2 = (
            nks[-1] for nks in (can_nkeys1, can_nkeys2, loc_nkeys1, loc_nkeys2)
        )

        flip_dct[bkey] = (loc_nmax1 != can_nmax1) ^ (loc_nmax2 != can_nmax2)

    return flip_dct


def substitution_reversal_parity_flips(
    tsg,
    keys: CenterKeys,
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
) -> Dict[AtomKey, bool]:
    r"""Parity flips for reversing the direction of a substitution TS graph


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
    :param pri_dct: The priority mapping
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


def vinyl_addition_reactant_parity_flips(
    gra,
    keys: CenterKeys,
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
    ts_nkeys_dct: CenterNeighborDict,
) -> Dict[BondKey, bool]:
    """Parity flips for the reactant parities of a vinyl addition TS graph

    :param gra: A graph of the reactants
    :type gra: automol graph data structure
    :param keys: The stereocenter keys at which the parity is being evaluated
    :type keys: CenterKeys
    :param pri_dct: The priority mapping
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for reactants graph
    :type nkeys_dct: CenterNeighborDict
    :param ts_nkeys_dct: Stereo-determining neighbors for the TS graph
    :type ts_nkeys_dct: CenterNeighborDict
    :return: A mapping identifying which stereoatoms flip upon reversal
    :rtype: Dict[AtomKey, bool]
    """
    vin_dct = vinyl_radical_atom_bond_keys(gra)

    vin_flip_dct = {}
    for vin_key, bkey in vin_dct.items():
        if bkey in keys:
            vin_idx = sorted(bkey).index(vin_key)
            nkeys0 = ts_nkeys_dct[bkey][vin_idx]
            nkeys = nkeys_dct[bkey][vin_idx]

            nmax0 = max(nkeys0, key=pri_dct.get)
            nmax = max(nkeys, key=pri_dct.get)

            vin_flip_dct[bkey] = nmax0 != nmax

    return vin_flip_dct


# def type_1_2_insertion_reactant_parity_flips(
#     gra,
#     keys: CenterKeys,
#     pri_dct: Dict[AtomKey, int],
#     nkeys_dct: CenterNeighborDict,
#     ts_nkeys_dct: CenterNeighborDict,
# ) -> Dict[BondKey, bool]:
#     """Parity flips for the reactant parities of a vinyl addition TS graph

#     :param gra: A graph of the reactants
#     :type gra: automol graph data structure
#     :param keys: The stereocenter keys at which the parity is being evaluated
#     :type keys: CenterKeys
#     :param pri_dct: The priority mapping
#     :type pri_dct: Dict[AtomKey, int]
#     :param nkeys_dct: Stereo-determining neighbors for reactants graph
#     :type nkeys_dct: CenterNeighborDict
#     :param ts_nkeys_dct: Stereo-determining neighbors for the TS graph
#     :type ts_nkeys_dct: CenterNeighborDict
#     :return: A mapping identifying which stereoatoms flip upon reversal
#     :rtype: Dict[AtomKey, bool]
#     """


# stereo parity evaluations
def geometry_atom_parity(gra, geo, atm_key, nkeys=None, geo_idx_dct=None):
    r""" Calculate an atom parity directly from a geometry

    Neighboring atom keys (`nkeys`) must be passed in as a priority-sorted
    list. If `None`, a local parity calculation will occur based on the
    atom keys in the molecular graph.

    Atom parity is defined as follows:

    The four keys passed in are apices of a tetrahedron. Looking at 2, 3,
    and 4 from 1, they will either ascend in clockwise or counterclockwise
    order.

    If ascending in counterclockwise order, the parity is False ('-').
    If ascending in clockwise order, the parity is True ('+').

            2                   2
           /1\                 /1\
          3---4               4---3

        counterclockwise    clockwise
        False               True
        '-'                 '+'

    (Viewed looking down from 1)

    If only three keys are passed in, they will be treated as keys 2, 3,
    and 4 above and it will be assumed that there is a lone pair at 1.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param atm_key: the atom key whose parity is being evaluated
    :type atm_key: int
    :param nkeys: the neighboring atom keys, pre-sorted by priority
    :type nkeys: list[int]
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    """
    assert gra == explicit(
        gra
    ), "Explicit graph should be used when getting parities from geometry."

    keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )

    nkeys = atom_stereo_sorted_neighbor_keys(gra, atm_key) if nkeys is None else nkeys

    # If there are only three groups, use the stereo atom itself as
    # the top apex of the tetrahedron.
    if len(nkeys) == 4:
        keys = nkeys
    else:
        assert len(nkeys) == 3
        keys = [atm_key] + list(nkeys)

    idxs = list(map(geo_idx_dct.__getitem__, keys))
    xyzs = geom_base.coordinates(geo, idxs=idxs)
    det_mat = numpy.ones((4, 4))
    det_mat[:, 1:] = xyzs
    det_val = numpy.linalg.det(det_mat)
    assert det_val != 0.0  # for now, assume no four-atom planes
    par = bool(det_val > 0.0)
    return par


def geometry_bond_parity(gra, geo, bnd_key, bnd_nkeys=None, geo_idx_dct=None):
    r""" Calculate a bond parity directly from a geometry

    Neighboring bond keys (`bnd_nkeys`) must be passed in as a pair of
    priority-sorted lists corresponding to the first and second atoms in
    `bnd_key`. Note that the latter must be an *ordered list* in this case!
    If `None`, a local parity calculation will occur based on the atom keys
    in the molecular graph.

    Bond parity is defined as follows:

    For each atom in the double bond, find the heavy-atom neighbor with the
    higher canonical number. Although hydrogen atoms have higher canonical
    numbers, they are always given lowest priority.

    If the neighbors are cis to each other, the parity is False ('-').
    If the neighbors are trans to each other, the parity is True ('+').

        max     max      max     min
           \   /            \   /
            A=B              A=B
           /   \            /   \
        min     min      min     max

        cis              trans
        False            True
        '-'              '+'

    If one side only has a single neighbor, then it is compared with the
    maximum neighbor on the other side.

        max     nei      max
           \   /            \
            A=B              A=B
           /                /   \
        min              min     nei

        cis              trans
        False            True
        '-'              '+'

    If both sides have only single neighbors, then they are compared to
    each other.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param bnd_key: the bond key. If using `bnd_nkeys`, this must be an
        ordered list!
    :type bnd_key: list[int]
    :param bnd_nkeys: a pair of lists of neighboring keys for the first and
        second atoms in `bnd_key`, respectively.
    :type bnd_nkeys: list[list[int]]
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    """
    assert gra == explicit(
        gra
    ), "Explicit graph should be used when getting parities from geometry."

    assert (
        isinstance(bnd_key, abc.Collection) and len(bnd_key) == 2
    ), f"{bnd_key} is not a valid bond key."
    key1, key2 = bnd_key

    keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )

    if bnd_nkeys is None:
        nkey1s, nkey2s = bond_stereo_sorted_neighbor_keys(gra, key1, key2)
    else:
        assert (
            isinstance(bnd_nkeys, abc.Collection)
            and len(bnd_nkeys) == 2
            and all(isinstance(nk, abc.Collection) for nk in bnd_nkeys)
        ), f"Bond neighbor keys should be a pair of lists: {bnd_nkeys}"
        nkey1s, nkey2s = bnd_nkeys

    idx1 = geo_idx_dct[key1]
    idx2 = geo_idx_dct[key2]
    nidx1 = geo_idx_dct[nkey1s[-1]]
    nidx2 = geo_idx_dct[nkey2s[-1]]

    (xyz1,) = geom_base.coordinates(geo, idxs=(idx1,))
    (xyz2,) = geom_base.coordinates(geo, idxs=(idx2,))
    (nxyz1,) = geom_base.coordinates(geo, idxs=(nidx1,))
    (nxyz2,) = geom_base.coordinates(geo, idxs=(nidx2,))

    bnd1_vec = numpy.subtract(nxyz1, xyz1)
    bnd2_vec = numpy.subtract(nxyz2, xyz2)

    dot_val = numpy.vdot(bnd1_vec, bnd2_vec)
    assert dot_val != 0.0  # for now, assume not collinear
    par = bool(dot_val < 0.0)
    return par


def parity_evaluator_measure_from_geometry_(
    geo,
    local_stereo: bool = False,
    geo_idx_dct: Optional[Dict[AtomKey, AtomKey]] = None,
) -> ParityEvaluator:
    """Creates a parity evaluator that measures parities from a geometry

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
        """A parity evaluator that measures parities from a geometry

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
        # 0. Determine local priories, if requesting local stereo
        pri_dct = local_stereo_priorities(gra) if local_stereo else pri_dct

        # 1. Get sorted neighbor keys for stereocenters
        nkeys_dct = dict_.by_key(nkeys_dct, keys)
        atm_nkeys_dct, bnd_nkeys_dct = stereocenter_candidates_grouped_and_sorted(
            nkeys_dct, pri_dct=pri_dct
        )

        # 2. Measure parities for atoms and bonds
        atm_par_dct = {
            k: geometry_atom_parity(gra, geo, k, ns, geo_idx_dct=geo_idx_dct)
            for k, ns in atm_nkeys_dct.items()
        }
        bnd_par_dct = {
            k: geometry_bond_parity(gra, geo, sorted(k), ns, geo_idx_dct=geo_idx_dct)
            for k, ns in bnd_nkeys_dct.items()
        }

        # 3. Combine parity mappings
        par_dct = {**atm_par_dct, **bnd_par_dct}

        # 4. Apply substitution reversal flips, if this is a reverse TS graph and we are
        # measuring local parities
        if local_stereo and is_rev_ts:
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
    """A parity evaluator that reads parities from a graph

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
    """A parity evaluator that reads parities from a graph and flips them to convert
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
    # 0. Read in the parities
    par_dct = dict_.by_key(stereo_parities(gra), keys)

    # 1. Determine local priorities
    loc_pri_dct = local_stereo_priorities(gra)

    # 2. Apply local parity flips to convert to/from canonical stereo
    loc_flip_dct = local_parity_flips(keys, pri_dct, loc_pri_dct, nkeys_dct)
    for key, flip in loc_flip_dct.items():
        par_dct[key] ^= flip

    # 3. Apply substitution reversal flips, if this is a reverse TS graph
    if is_rev_ts:
        subst_flip_dct = substitution_reversal_parity_flips(
            gra, keys, pri_dct, nkeys_dct
        )
        for key, flip in subst_flip_dct.items():
            par_dct[key] ^= flip

    return par_dct


def parity_evaluator_reactants_from_local_ts_graph_(
    loc_tsg,
    local_stereo: bool = False,
) -> ParityEvaluator:
    """Creates a parity evaluator that evaluates reactant parities from a local TS graph

    :param loc_tsg: A TS graph with local parities
    :type loc_tsg: automol graph data structure
    :param local_stereo: Return local stereo assignments? defaults to False
    :type local_stereo: bool, optional
    :returns: A parity evaluator
    :rtype: ParityEvaluator
    """
    ts_loc_par_dct = stereo_parities(loc_tsg)
    ts_nkeys_dct = stereocenter_candidates(loc_tsg)

    def parity_evaluator(
        gra,
        keys: CenterKeys,
        pri_dct: Dict[AtomKey, int],
        nkeys_dct: CenterNeighborDict,
        is_rev_ts: bool = False,
    ):
        """A parity evaluator that evaluates reactant parities from a local TS graph

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
        assert is_rev_ts or not is_rev_ts

        # 0. Read in local TS parities
        par_dct = dict_.by_key(ts_loc_par_dct, keys)

        # 1. Determine local priorities
        loc_pri_dct = local_stereo_priorities(gra)

        # 2. Correct vinyl addition parity flips
        vin_flip_dct = vinyl_addition_reactant_parity_flips(
            gra, keys, loc_pri_dct, nkeys_dct, ts_nkeys_dct=ts_nkeys_dct
        )
        for key, flip in vin_flip_dct.items():
            par_dct[key] ^= flip

        # 3. Apply local parity flips to convert to canonical stereo, if requested
        if not local_stereo:
            loc_flip_dct = local_parity_flips(keys, pri_dct, loc_pri_dct, nkeys_dct)
            for key, flip in loc_flip_dct.items():
                par_dct[key] ^= flip

        return par_dct

    return parity_evaluator
