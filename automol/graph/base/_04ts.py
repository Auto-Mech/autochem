"""TS classification and other functions
"""
import itertools
from typing import Dict, List, Tuple

from automol import util
from automol.graph.base._00core import (
    atom_keys,
    atom_stereo_sorted_neighbor_keys,
    atoms_neighbor_atom_keys,
    bond_stereo_keys,
    bond_stereo_sorted_neighbor_keys,
    is_ts_graph,
    local_stereo_priorities,
    sort_by_size,
    stereo_parities,
    tetrahedral_atom_keys,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reacting_bond_keys,
    ts_reagents_graph_without_stereo,
    ts_reverse,
)
from automol.graph.base._02algo import (
    connected_components,
    rings_bond_keys,
    sorted_ring_atom_keys_from_bond_keys,
)
from automol.graph.base._03kekule import (
    rigid_planar_bond_keys,
    vinyl_radical_atom_bond_keys,
)


def is_bimolecular(tsg) -> bool:
    """Is this a TS for a bimolecular reaction?

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: `True` if it is, `False` if it isn't
    :rtype: bool
    """
    return len(connected_components(ts_reagents_graph_without_stereo(tsg))) == 2


def has_reacting_ring(tsg) -> bool:
    """Does this TS graph have a ring which is involved in the reaction?

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return bool(forming_rings_atom_keys(tsg) or breaking_rings_atom_keys(tsg))


def atom_transfers(tsg) -> Dict[int, Tuple[int, int]]:
    """Get a dictionary describing atom transfers; keys are transferring atoms, values
    are donors and acceptors, respectively

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A list of triples containing the donor atom, the transferring atom, and
        the acceptor atom, respectively
    :rtype: Dict[int, Tuple[int, int]]
    """
    brk_bkeys = ts_breaking_bond_keys(tsg)
    frm_bkeys = ts_forming_bond_keys(tsg)

    tra_dct = {}
    for brk_bkey, frm_bkey in itertools.product(brk_bkeys, frm_bkeys):
        if brk_bkey & frm_bkey:
            (tra_key,) = brk_bkey & frm_bkey
            (don_key,) = brk_bkey - frm_bkey
            (acc_key,) = frm_bkey - brk_bkey
            tra_dct[tra_key] = (don_key, acc_key)

    return tra_dct


def sn2_atom_transfers(tsg) -> Dict[int, Tuple[int, int]]:
    """Get a dictionary describing atom transfers for Sn2 reactions; keys are the
    transferring atoms, values are the donors and acceptors, respectively

    Identifies transferring atoms that are tetrahedral

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A list of triples containing the donor atom, the transferring atom, and
        the acceptor atom, respectively
    :rtype: Dict[int, Tuple[int, int]]
    """
    tra_dct = atom_transfers(tsg)
    tra_keys = set(tra_dct.keys())
    tet_keys = tetrahedral_atom_keys(tsg)
    sn2_keys = tra_keys & tet_keys
    return util.dict_.by_key(tra_dct, sn2_keys)


def zmatrix_sorted_reactants_keys(tsg) -> List[List[int]]:
    """For bimolecular reactions without a TS ring, return keys for the reactants in the
    order they should appear in the z-matrix

    For unimolecular reactions or bimolecular reactions with a TS ring, returns None

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: For bimolecular reactions without a TS ring, the keys for each reactant,
        in z-matrix order, i.e. atom donor first or, if there isn't one, larger reactant
        first, as measured by (heavy atoms, total atoms, electrons)
    :rtype: List[int]
    """
    if not is_ts_graph(tsg) or has_reacting_ring(tsg) or not is_bimolecular(tsg):
        return None

    rcts_gra = ts_reagents_graph_without_stereo(tsg, dummy=True)
    rct_gras = connected_components(rcts_gra)

    # 1. If there is an atom transfer, put the donor reagent first
    tra_keys = set(atom_transfers(tsg))
    if tra_keys:
        rct_gras = sorted(rct_gras, key=lambda g: atom_keys(g) & tra_keys, reverse=True)
    # 2. Otherwise, put the larger reagent first
    else:
        rct_gras = sort_by_size(rct_gras)

    rcts_keys = tuple(map(tuple, map(sorted, map(atom_keys, rct_gras))))
    return rcts_keys


def zmatrix_starting_ring_keys(tsg) -> List[int]:
    """Return keys for a TS ring to start from, sorted in z-matrix order

    If there isn't a TS ring, this returns `None`

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: The ring keys, sorted to exclude breaking bonds and include forming bonds
        as late as possible
    :rtype: List[int]
    """
    if not has_reacting_ring(tsg):
        return None

    rngs_keys = reacting_rings_atom_keys(tsg)

    if len(rngs_keys) > 1:
        raise NotImplementedError(f"Not implemented for multiple reacting rings: {tsg}")

    (rng_keys,) = rngs_keys
    brk_bkeys = {bk for bk in ts_breaking_bond_keys(tsg) if bk < set(rng_keys)}
    frm_bkeys = {bk for bk in ts_forming_bond_keys(tsg) if bk < set(rng_keys)}
    frm_keys = list(itertools.chain(*frm_bkeys))

    return util.ring.cycle_to_optimal_split(rng_keys, brk_bkeys, frm_keys)


def reacting_rings_atom_keys(tsg) -> List[List[int]]:
    """Get the atom keys to rings containing breaking or forming bonds

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    rngs_bnd_keys = reacting_rings_bond_keys(tsg)
    rngs_atm_keys = tuple(map(sorted_ring_atom_keys_from_bond_keys, rngs_bnd_keys))
    return rngs_atm_keys


def reacting_rings_bond_keys(tsg):
    """get the bond keys to rings forming in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    bnd_keys = ts_reacting_bond_keys(tsg)
    rngs_bnd_keys = tuple(
        bks for bks in rings_bond_keys(tsg, ts_=True) if bnd_keys & bks
    )
    return rngs_bnd_keys


def forming_rings_atom_keys(tsg) -> List[List[int]]:
    """Get the atom keys to rings forming in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    frm_rngs_bnd_keys = forming_rings_bond_keys(tsg)
    frm_rngs_atm_keys = tuple(
        map(sorted_ring_atom_keys_from_bond_keys, frm_rngs_bnd_keys)
    )
    return frm_rngs_atm_keys


def forming_rings_bond_keys(tsg):
    """get the bond keys to rings forming in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    frm_bnd_keys = ts_forming_bond_keys(tsg)
    frm_rngs_bnd_keys = tuple(
        bks for bks in rings_bond_keys(tsg, ts_=True) if frm_bnd_keys & bks
    )
    return frm_rngs_bnd_keys


def breaking_rings_atom_keys(tsg):
    """get the atom keys to rings breaking in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    brk_rngs_bnd_keys = breaking_rings_bond_keys(tsg)
    brk_rngs_atm_keys = tuple(
        map(sorted_ring_atom_keys_from_bond_keys, brk_rngs_bnd_keys)
    )
    return brk_rngs_atm_keys


def breaking_rings_bond_keys(tsg):
    """get the bond keys to rings breaking in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    brk_bnd_keys = ts_breaking_bond_keys(tsg)
    brk_rngs_bnd_keys = tuple(
        bks for bks in rings_bond_keys(tsg, ts_=True) if brk_bnd_keys & bks
    )
    return brk_rngs_bnd_keys


def sn2_local_stereo_reversal_flips(tsg) -> Dict[int, bool]:
    r"""For Sn2 reaction sites, identifies for which of them the local stereo flips upon
    reversing the TS direction

        Case 1: Reversal causes local parity flip (reversal value: `True`):

                  3                     3
                  |                     |
            5-----1  +  6   <=>   5  +  1-----6
                 / \                   / \
                4   2                 4   2

                clockwise             counterclockwise
                ('+')                 ('-')

        Case 2: Reversal leaves local parity unchanged (reversal value: `False`):

                  3                     3
                  |                     |
            5-----1  +  6   <=>   5  +  1-----6
                 / \                   / \
                4   2                 4   2

                clockwise             counterclockwise
                ('+')                 ('-')

        Given a reversal value, r_flip, the local parity upon reversing TS direction is
        given by the following:

            p_rev = p_forw ^ r_flip

    :param tsg: A TS graph, with or without stereo (stereo is ignored)
    :type tsg: automol graph data structure
    :return: Which Sn2 sites flip stereo upon reversing TS direction
    :rtype: Dict[int, bool]
    """
    rflip_dct = {}
    for tra_key, (don_key, acc_key) in sn2_atom_transfers(tsg).items():
        # Get the forward direction neighboring keys, sorted by local priority
        nks0 = list(atom_stereo_sorted_neighbor_keys(tsg, tra_key))
        # Replace the donor with the acceptor
        nks0[nks0.index(don_key)] = acc_key
        # Get the reverse direction neighboring keys, sorted by local priority
        nks1 = atom_stereo_sorted_neighbor_keys(ts_reverse(tsg), tra_key)
        # If the orderings are related by an even permuation, the parity will flip (due
        # to umbrella inversion of the stereo center, see ASCII diagrams above)
        rflip_dct[tra_key] = util.is_even_permutation(nks0, nks1)
    return rflip_dct


def constrained_1_2_insertion_local_parities(loc_tsg) -> Dict[frozenset, bool]:
    r"""Calculates the local parities of the reactant of a constrained 1,2-insertion, if
    present

    In general, there should only be one, but we allow the possibility of multiples for
    consistency.

        Constraint:

             4         8              3           6
              \   +   /     3          \ +     + /
               1=====2      |  <=>      1-------2
              /       \     6          / \     / \
             5         7              4   5   8   7

            trans                 clockwise  clockwise
            ('+')                 ('+')      ('+')

        Equation: p_b12 = (p_a1 AND p_1) XNOR (p_a2 AND p_2)

        Where p_a1 and p_a2 are the parities of the atoms, p_b is the resulting
        parity of the 1=2 double bond, and p_1 and p_2 are the parities of the
        permutations required to move the other atom in the bond first in the
        list and the leaving atom second for each atom.

    :param loc_tsg: TS graph with local stereo parities
    :type loc_tsg: automol graph data structure
    :return: The local parities of bonds subject to the constraint
    :rtype: Dict[frozenset, bool]
    """
    loc_par_dct = util.dict_.filter_by_value(
        stereo_parities(loc_tsg), lambda x: x is not None
    )
    loc_pri_dct = local_stereo_priorities(loc_tsg)
    nkeys_dct = atoms_neighbor_atom_keys(loc_tsg)

    # Check first reactants, then products
    gra = ts_reagents_graph_without_stereo(loc_tsg, prod=False)
    frm_bkeys = ts_forming_bond_keys(loc_tsg)
    rkeys_lst = list(map(set, forming_rings_atom_keys(loc_tsg)))

    par_dct = {}
    # Conditions:
    # A. Bond is rigid and planar
    for bkey in rigid_planar_bond_keys(gra):
        # Conditions:
        # B. Both atoms in the bond are stereogenic
        if all(k in loc_par_dct for k in bkey):
            key1, key2 = bkey
            (key0,) = next(k for k in frm_bkeys if key1 in k) - {key1}
            (key3,) = next(k for k in frm_bkeys if key2 in k) - {key2}
            keys = {key0, key1, key2, key3}
            # C. The atoms are part of a forming ring
            if any(keys & rkeys for rkeys in rkeys_lst):
                # This is a constrained bond!
                # i. Read out the relevant atom and bond parities
                p_a1 = loc_par_dct[key1]
                p_a2 = loc_par_dct[key2]

                # ii. Calculate the local bond parity
                nk1s = sorted(nkeys_dct[key1], key=loc_pri_dct.__getitem__)
                nk2s = sorted(nkeys_dct[key2], key=loc_pri_dct.__getitem__)
                srt_nk1s = util.move_items_to_front(nk1s, [key2, key0])
                srt_nk2s = util.move_items_to_front(nk2s, [key1, key3])
                p_1 = util.is_even_permutation(nk1s, srt_nk1s)
                p_2 = util.is_even_permutation(nk2s, srt_nk2s)

                # p_b12 = (p_a1 XNOR p_1) XNOR (p_a2 XNOR p_2)
                p_b12 = not ((not (p_a1 ^ p_1)) ^ (not (p_a2 ^ p_2)))

                par_dct[bkey] = p_b12

    return par_dct


def vinyl_addition_local_parities(loc_tsg) -> Dict[frozenset, bool]:
    r""" Calculates the local parity of the reactant or product of a
    vinyl addition, if present

        Constraint:

                 5                        5         6
                  \   +                    \   -   /
                   1=====2   +  6   <=>     1=====2
                  /       \                /       \
                 4         3              4         3

                trans                    cis
                ('+')                    ('-')

    :param loc_tsg: TS graph with local stereo parities
    :type loc_tsg: automol graph data structure
    :return: The local parities of bonds subject to the constraint
    :rtype: Dict[frozenset, bool]
    """
    ste_bkeys = bond_stereo_keys(loc_tsg)
    loc_par_dct = util.dict_.filter_by_value(
        stereo_parities(loc_tsg), lambda x: x is not None
    )
    loc_pri_dct = local_stereo_priorities(loc_tsg)

    gra = ts_reagents_graph_without_stereo(loc_tsg, prod=False)
    vin_dct = vinyl_radical_atom_bond_keys(gra)

    par_dct = {}

    for vin_akey, vin_bkey in vin_dct.items():
        if vin_bkey in ste_bkeys:
            par = loc_par_dct[vin_bkey]

            (vin_nkey,) = vin_bkey - {vin_akey}

            nkeys_ts, _ = bond_stereo_sorted_neighbor_keys(
                loc_tsg, vin_akey, vin_nkey, pri_dct=loc_pri_dct
            )
            nkeys_r, _ = bond_stereo_sorted_neighbor_keys(
                gra, vin_akey, vin_nkey, pri_dct=loc_pri_dct
            )

            # If the maximum priority neighbors don't match, flip
            # the parity
            if nkeys_r and nkeys_r[-1] != nkeys_ts[-1]:
                # Flip the parity
                par = not par

                par_dct[vin_bkey] = par

    return par_dct
