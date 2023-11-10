"""TS classification and other functions
"""
import itertools
from typing import Dict, List, Tuple

from automol import util
from automol.graph.base._00core import (
    atom_keys,
    atoms_neighbor_atom_keys,
    bond_stereo_keys,
    bond_stereo_sorted_neighbor_keys,
    is_ts_graph,
    local_stereo_priorities,
    sort_by_size,
    stereo_parities,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reacting_bond_keys,
    ts_reagents_graph_without_stereo,
)
from automol.graph.base._02algo import (
    connected_components,
    rings_bond_keys,
    sorted_ring_atom_keys_from_bond_keys,
)
from automol.graph.base._03kekule import rigid_planar_bond_keys, vinyl_radical_atom_keys


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


def constrained_1_2_insertion_local_parities(loc_tsg):
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
    :return: The key of the stereo site, its local parity, and a boolean
        indicating whether or not this is for the products
    :rtype: frozenset({int, int}, bool, bool
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


def vinyl_addition_local_parities(loc_tsg):
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
    :return: The key of the stereo site, its local parity, and a boolean
        indicating whether or not this is for the products
    :rtype: frozenset({int, int}, bool, bool
    """
    bkeys = bond_stereo_keys(loc_tsg)
    loc_par_dct = util.dict_.filter_by_value(
        stereo_parities(loc_tsg), lambda x: x is not None
    )
    loc_pri_dct = local_stereo_priorities(loc_tsg)

    gra = ts_reagents_graph_without_stereo(loc_tsg, prod=False)
    vin_keys = vinyl_radical_atom_keys(gra)

    par_dct = {}
    # Identify stereogenic bonds with vinyl radical atoms
    bkeys = [bk for bk in bkeys if bk & vin_keys]
    for bkey in bkeys:
        par = loc_par_dct[bkey]

        # Get sorted neighbors for both atoms
        akeys = sorted(bkey)
        tnkeys_pair = bond_stereo_sorted_neighbor_keys(
            loc_tsg, *akeys, pri_dct=loc_pri_dct
        )
        gnkeys_pair = bond_stereo_sorted_neighbor_keys(gra, *akeys, pri_dct=loc_pri_dct)
        for akey, tnkeys, gnkeys in zip(akeys, tnkeys_pair, gnkeys_pair):
            # Compare sorted neighbors for mismatch
            if gnkeys and tnkeys != gnkeys:
                assert akey in vin_keys, (
                    f"Neighbor mismatch at {akey} is not due to "
                    f"vinyl addition:\n{loc_tsg}"
                )
                # If the maximum priority neighbors don't match, flip
                # the parity
                if tnkeys[-1] != gnkeys[-1]:
                    # Flip the parity
                    par = not par

                    par_dct[bkey] = par

    return par_dct
