"""TS classification and other functions
"""
import functools
import itertools
import warnings
from typing import Dict, Optional, Tuple

import numpy
from automol import util
from automol.graph.base._00core import (
    atom_implicit_hydrogens,
    atom_neighbor_atom_keys,
    atoms_neighbor_atom_keys,
    bond_neighbor_atom_keys,
    bond_stereo_keys,
    is_ts_graph,
    local_stereo_priorities,
    stereo_parities,
    tetrahedral_atom_keys,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reactants_graph_without_stereo,
    ts_reagents_graphs_without_stereo,
    ts_reverse,
    without_bonds_by_orders,
    without_dummy_atoms,
)
from automol.graph.base._02algo import forming_rings_atom_keys, reacting_rings_bond_keys
from automol.graph.base._03kekule import (
    rigid_planar_bond_keys,
    vinyl_radical_atom_bond_keys,
)
from automol.util import dict_


# reaction site classification
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


def substitutions(tsg) -> Dict[int, Tuple[int, int]]:
    """Get a dictionary describing substitution reaction sites

    Maps transferring atoms onto their leaving and entering atoms, respectively

    (Limited to substitutions at tetrahedral atoms)

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A mapping of transferring atoms onto leaving and entering atoms
    :rtype: Dict[int, Tuple[int, int]]
    """
    tra_dct = atom_transfers(tsg)
    tra_keys = set(tra_dct.keys())
    tet_keys = tetrahedral_atom_keys(tsg)
    subst_keys = tra_keys & tet_keys
    return util.dict_.by_key(tra_dct, subst_keys)


def eliminations(tsg) -> Dict[int, Tuple[int, int]]:
    """Get a dictionary describing elimination reaction sites

    Maps bonds across which eliminations occur onto their leaving atoms, along with the
    forming bond key for each, if present

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A mapping of elimination bonds onto leaving atoms and forming bond keys
        (Leaving atoms are sorted in order of the elimination bond atoms)
    :rtype: Dict[int, Tuple[int, int]]
    """
    brk_bkeys_pool = ts_breaking_bond_keys(tsg)
    frm_bkeys_pool = ts_forming_bond_keys(tsg)
    rng_bkeys_lst = reacting_rings_bond_keys(tsg)

    def is_elimination_bond_(brk_bkeys, frm_bkeys):
        """An elimination bond key is a ring key that intersects two breaking bonds and
        is not a forming bond
        """

        def is_elimination_bkey(bkey):
            return all(bkey & bk for bk in brk_bkeys) and bkey not in frm_bkeys

        return is_elimination_bkey

    def common_atom(bkey1, bkey2):
        return next(iter(bkey1 & bkey2))

    # 1. Check reacting rings
    elim_dct = {}
    for rng_bkeys in rng_bkeys_lst:
        brk_bkeys = brk_bkeys_pool & rng_bkeys
        frm_bkeys = frm_bkeys_pool & rng_bkeys

        # 2. Require two breaking bonds within the ring
        if len(brk_bkeys) == 2:
            # 3. Find the elimination bond
            is_elim = is_elimination_bond_(brk_bkeys, frm_bkeys)
            bkey = next(filter(is_elim, rng_bkeys), None)
            if bkey is not None:
                # a. Sort the breaking bonds in order of the elimination bond atom
                brk_bkey1, brk_bkey2 = sorted(
                    brk_bkeys, key=functools.partial(common_atom, bkey)
                )

                # b. Get the corresponding leaving atoms
                (lea_key1,) = brk_bkey1 - bkey
                (lea_key2,) = brk_bkey2 - bkey

                # c. Get the forming bond key, if there is one
                assert len(frm_bkeys) <= 1, "Unexpected multiple forming bonds:{tsg}"
                frm_bkey = next(iter(frm_bkeys), None)

                elim_dct[bkey] = (lea_key1, lea_key2, frm_bkey)

    return elim_dct


def insertions(tsg) -> Dict[int, Tuple[int, int]]:
    """Get a dictionary describing insertion reaction sites

    Maps bonds across which insertions occur onto their entering atoms, along with the
    breaking bond key for each, if present

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A mapping of insertion bonds onto entering atoms and breaking bond keys
        (Entering atoms are sorted in order of the insertion bond atoms)
    :rtype: Dict[int, Tuple[int, int]]
    """
    return eliminations(ts_reverse(tsg))


# vvv DEPRECATED vvv
def atom_stereo_sorted_neighbor_keys(
    gra, key, self_apex: bool = False, pri_dct: Optional[Dict[int, int]] = None
):
    """Get keys for the neighbors of an atom that are relevant for atom
    stereochemistry, sorted by priority (if requested)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param key: the atom key
    :type key: int
    :param self_apex: If there are only 3 neighbors, put this atom as the apex?
    :type self_apex: bool, optional
    :param pri_dct: Priorities to sort by (optional)
    :type pri_dct: Optional[Dict[int, int]]
    :returns: The keys of neighboring atoms
    :rtype: tuple[int]
    """
    gra = without_dummy_atoms(gra)
    nhyd_dct = atom_implicit_hydrogens(gra)
    pri_dct = local_stereo_priorities(gra) if pri_dct is None else pri_dct

    # If this is an Sn2 stereocenter, use the reactants graph
    if key in substitutions(gra):
        gra = ts_reactants_graph_without_stereo(gra)

    # Get the neighboring atom keys
    nkeys = list(atom_neighbor_atom_keys(gra, key))

    # Add Nones for the implicit hydrogens
    nkeys.extend([None] * nhyd_dct[key])

    # Sort them by priority
    nkeys = sorted(nkeys, key=dict_.sort_value_(pri_dct, missing_val=-numpy.inf))

    # Optionally, if there are only three groups, use the stereo atom itself as
    # the top apex of the tetrahedron
    if self_apex and len(nkeys) < 4:
        assert len(nkeys) == 3
        nkeys = [key] + list(nkeys)

    return tuple(nkeys)


def bond_stereo_sorted_neighbor_keys(
    gra, key1, key2, pri_dct: Optional[Dict[int, int]] = None
):
    """Get keys for the neighbors of a bond that are relevant for bond
    stereochemistry, sorted by priority (if requested)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param key1: the first atom in the bond
    :type key1: int
    :param key2: the second atom in the bond
    :type key2: int
    :param pri_dct: Priorities to sort by (optional)
    :type pri_dct: Optional[Dict[int, int]]
    :returns: The keys of neighboring atoms for the first and second atoms
    :rtype: tuple[int], tuple[int]
    """
    gra = without_dummy_atoms(gra)
    nhyd_dct = atom_implicit_hydrogens(gra)
    pri_dct = local_stereo_priorities(gra) if pri_dct is None else pri_dct

    gras = ts_reagents_graphs_without_stereo(gra) if is_ts_graph(gra) else [gra]

    nkeys1 = set()
    nkeys2 = set()
    # For TS graphs, loop over reactants and products
    for gra_ in gras:
        # Check that the bond is rigid and planar on this side of the reaction, by
        # checking for a tetrahedral atom
        tet_keys = tetrahedral_atom_keys(gra_)
        if key1 not in tet_keys and key2 not in tet_keys:
            # Add these neighboring keys to the list
            nkeys1_, nkeys2_ = bond_neighbor_atom_keys(gra_, key1, key2)
            nkeys1.update(nkeys1_)
            nkeys2.update(nkeys2_)

    nkeys1 = list(nkeys1)
    nkeys2 = list(nkeys2)
    nkeys1.extend([None] * nhyd_dct[key1])
    nkeys2.extend([None] * nhyd_dct[key2])

    # Check that we don't have more than two neighbors on either side
    if len(nkeys1) > 2 or len(nkeys2) > 2:
        warnings.warn(
            f"Unusual neighbor configuration at bond {key1}-{key2} may result in "
            f"incorrect bond stereochemistry handling for this graph:\n{gra}"
        )

        # Temporary patch for misidentified substitutions at double bonds, which are
        # really two-step addition-eliminations
        gra_ = without_bonds_by_orders(gra, [0.1, 0.9])
        nkeys1, nkeys2 = map(list, bond_neighbor_atom_keys(gra_, key1, key2))
        nkeys1.extend([None] * nhyd_dct[key1])
        nkeys2.extend([None] * nhyd_dct[key2])

    nkeys1 = sorted(nkeys1, key=dict_.sort_value_(pri_dct, missing_val=-numpy.inf))
    nkeys2 = sorted(nkeys2, key=dict_.sort_value_(pri_dct, missing_val=-numpy.inf))
    return (tuple(nkeys1), tuple(nkeys2))


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
    for tra_key, (don_key, acc_key) in substitutions(tsg).items():
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
    gra = ts_reactants_graph_without_stereo(loc_tsg)
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

    gra = ts_reactants_graph_without_stereo(loc_tsg)
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
