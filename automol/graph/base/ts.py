""" transition state graph data structure

DEPRECATED (under construction to phase out)

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import qcelemental as qcel

import automol.amchi.base  # !!!!
from automol import util
from automol.graph.base._algo import (
    rings_bond_keys,
    sorted_ring_atom_keys_from_bond_keys,
)
from automol.graph.base._canon import (
    calculate_priorities_and_assign_stereo,
    local_priority_dict,
    parity_evaluator_flip_local_,
    to_local_stereo,
)
from automol.graph.base._core import (
    atoms_neighbor_atom_keys,
    bond_stereo_keys,
    bond_stereo_sorted_neighbor_keys,
    has_stereo,
    set_stereo_parities,
    stereo_keys,
    stereo_parities,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_graph,
    ts_reacting_atom_keys,
    ts_reacting_bond_keys,
    ts_reagents_graph_without_stereo,
    ts_reverse,
)
from automol.graph.base._heur import heuristic_bond_distance as _heuristic_bond_distance
from automol.graph.base._kekule import rigid_planar_bond_keys, vinyl_radical_atom_keys
from automol.graph.base._stereo import (
    expand_stereo,
)

# Rename TS-specific functions defined elsewhere
graph = ts_graph
forming_bond_keys = ts_forming_bond_keys
breaking_bond_keys = ts_breaking_bond_keys
reacting_bond_keys = ts_reacting_bond_keys
reacting_atom_keys = ts_reacting_atom_keys
reverse = ts_reverse
reagents_graph_without_stereo = ts_reagents_graph_without_stereo


def forming_rings_atom_keys(tsg):
    """get the atom keys to rings forming in the TS graph

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


def reactants_graph(tsg, stereo=True):
    """Get the reactants from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    return reagents_graph(tsg, prod=False, stereo=stereo)


def products_graph(tsg, stereo=True):
    """Get the products from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    return reagents_graph(tsg, prod=True, stereo=stereo)


def reagents_graph(tsg, prod=False, stereo=True):
    """Get the reactants or products from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param prod: Do this for the products, instead of the reactants?
    :type prod: bool
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    gra = ts_reagents_graph_without_stereo(tsg, prod=prod)
    if stereo and has_stereo(tsg):
        _, gra = calculate_priorities_and_assign_stereo(
            gra,
            backbone_only=False,
            break_ties=False,
            par_eval_=parity_evaluator_reagents_from_ts_(tsg, prod=prod),
        )
    return gra


def parity_evaluator_reagents_from_ts_(tsg, prod=False):
    r"""Determines reactant or product stereochemistry from a TS graph

    (For internal use by the calculate_priorities_and_assign_stereo() function)

    Sn2 constraint is taken care of by parity_evaluator_flip_local_

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param prod: Do this for the products, instead of the reactants?
    :type prod: bool
    :returns: A parity evaluator, `p_`, for which `p_(gra, pri_dct)(key)`
        returns the parity for a given atom, given a set of priorities.
    """
    # Handle Sn2 reactions by reversing before localizing, if getting products
    tsg0 = ts_reverse(tsg) if prod else tsg
    loc_tsg0 = to_local_stereo(tsg0)

    # Handle constrained eliminations
    ckey, cpar, cprod = constrained_elimination_local_parity(loc_tsg0)
    if cprod is False:
        loc_tsg0 = set_stereo_parities(loc_tsg0, {ckey: cpar})

    # Handle vinyl radical additions
    vkey, vpar, vprod = vinyl_addition_local_parity(loc_tsg0)
    if vprod is False:
        loc_tsg0 = set_stereo_parities(loc_tsg0, {vkey: vpar})

    # Now that we have handled the exceptions, the local parities correspond to
    # what they will be for the reactants/products graph, so we can simply flip
    # the local stereo to find canonical assignments
    par_eval_flip_ = parity_evaluator_flip_local_()

    def _evaluator(gra, pri_dct, ts_rev=False):
        """Parity evaluator based on current priorities

        Note: `ts_rev` gets ignored here, because we are *not* assigning stereo
        to a TS graph (We are using TS stereo to assign to a non-TS graph)

        :param gra: molecular graph with canonical stereo parities
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities
        :type pri_dct: dict
        :param ts_rev: Is this a reversed TS graph?
        :type ts_rev: bool
        """

        # Sanity check
        assert gra == ts_reagents_graph_without_stereo(tsg0)

        # Do-nothing line to prevent linting complaint
        assert ts_rev or not ts_rev

        p0_ = par_eval_flip_(loc_tsg0, pri_dct)

        def _parity(key):
            # Stereocenters directly shared with TS can be obtained from the
            # flip local parity evaluator
            par = p0_(key)

            if par is None:
                raise NotImplementedError(
                    f"Reagent determination for this TS graph is not "
                    f"yet implemented:\n{tsg}"
                )

            return par

        return _parity

    return _evaluator


def expand_stereo_for_reaction(tsg, rcts_gra, prds_gra):
    """Expand TS graph stereo that is consistent with reactants and products

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param rcts_gra: reactants graph
    :type rcts_gra: automol graph data structure
    :param prds_gra: products graph
    :type prds_gra: automol graph data structure
    :return: A list of TS graphs with stereo assignments
    """
    # Identify stereo atoms and bonds
    fkeys = list(stereo_keys(rcts_gra))
    rkeys = list(stereo_keys(prds_gra))

    # Convert to local assignments
    loc_fgra = to_local_stereo(rcts_gra)
    loc_rgra = to_local_stereo(prds_gra)

    # Identify local parities
    fpars0 = util.dict_.values_by_key(stereo_parities(loc_fgra), fkeys)
    rpars0 = util.dict_.values_by_key(stereo_parities(loc_rgra), rkeys)
    assert not any(
        p is None for p in fpars0 + rpars0
    ), "Corrupted stereo assignments: {rcts_gra}\n{prds_gra}"

    # Include *all* possibilities for the TS, including symmetry equivalent
    # ones
    ste_tsgs = []
    for ste_tsg in expand_stereo(tsg, enant=True, symeq=True):
        # Convert forward and reverse graphs to local assignments
        loc_ftsg = to_local_stereo(ste_tsg)
        loc_rtsg = to_local_stereo(reverse(ste_tsg))

        # Handle constrained eliminations
        ckey, cpar, cprod = constrained_elimination_local_parity(loc_ftsg)
        if cprod is False:
            loc_ftsg = set_stereo_parities(loc_ftsg, {ckey: cpar})
        if cprod is True:
            loc_rtsg = set_stereo_parities(loc_rtsg, {ckey: cpar})

        # Handle vinyl additions
        vkey, vpar, vprod = vinyl_addition_local_parity(loc_ftsg)
        if vprod is False:
            loc_ftsg = set_stereo_parities(loc_ftsg, {vkey: vpar})
        if vprod is True:
            loc_rtsg = set_stereo_parities(loc_rtsg, {vkey: vpar})

        # Identify local parities
        fpars1 = util.dict_.values_by_key(stereo_parities(loc_ftsg), fkeys)
        rpars1 = util.dict_.values_by_key(stereo_parities(loc_rtsg), rkeys)
        # If they match the reactants and products, include these TS
        # assignments
        if fpars0 == fpars1 and rpars0 == rpars1:
            ste_tsgs.append(ste_tsg)

    return tuple(ste_tsgs)


def constrained_elimination_local_parity(loc_tsg):
    r"""Calculates the local parity of the reactant or product of a
    constrained elimination, if present

        Constraint:

                3           6         4         8
                 \ +     + /           \   +   /     3
                  1-------2      <=>    1=====2      |
                 / \     / \           /       \     6
                4   5   8   7         5         7

            clockwise  clockwise     trans
            ('+')      ('+')         ('+')

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
    loc_pri_dct = local_priority_dict(loc_tsg)
    nkeys_dct = atoms_neighbor_atom_keys(loc_tsg)

    # Check first reactants, then products
    for prod in [False, True]:
        gra = ts_reagents_graph_without_stereo(loc_tsg, prod=prod)
        if prod:
            frm_bkeys = ts_breaking_bond_keys(loc_tsg)
            rkeys_lst = list(map(set, breaking_rings_atom_keys(loc_tsg)))
        else:
            frm_bkeys = ts_forming_bond_keys(loc_tsg)
            rkeys_lst = list(map(set, forming_rings_atom_keys(loc_tsg)))

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

                    # There should only be one, so quit if we find it
                    return bkey, p_b12, prod

    return (None, None, None)


def vinyl_addition_local_parity(loc_tsg):
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
    loc_pri_dct = local_priority_dict(loc_tsg)

    for prod in [False, True]:
        gra = ts_reagents_graph_without_stereo(loc_tsg, prod=prod)
        vin_keys = vinyl_radical_atom_keys(gra)

        # Identify stereogenic bonds with vinyl radical atoms
        bkey = next((bk for bk in bkeys if bk & vin_keys), None)
        if bkey is not None:
            par = loc_par_dct[bkey]

            # Get sorted neighbors for both atoms
            akeys = sorted(bkey)
            tnkeys_pair = bond_stereo_sorted_neighbor_keys(
                loc_tsg, *akeys, pri_dct=loc_pri_dct
            )
            gnkeys_pair = bond_stereo_sorted_neighbor_keys(
                gra, *akeys, pri_dct=loc_pri_dct
            )
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

                        # There should only be one, so quit if we find it
                        return bkey, par, prod

    return (None, None, None)


def heuristic_bond_distance(
    tsg,
    key1: int,
    key2: int,
    fdist_factor: float = 1.1,
    bdist_factor: float = 0.9,
    angstrom: bool = True,
    check: bool = False,
) -> float:
    """The heuristic bond distance between two bonded or reacting atoms

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param key1: The first atom key
    :type key1: int
    :param key2: The second atom key
    :type key2: int
    :param fdist_factor: Set the forming bond distance to this times the average
        van der Waals radius, defaults to 1.1
    :type fdist_factor: float, optional
    :param bdist_factor: Set the breaking bond distance to this times the average
        van der Waals radius, defaults to 0.9
    :param angstrom: Return in angstroms intead of bohr?, defaults to True
    :type angstrom: bool, optional
    :param check: Check that these atoms are in fact bonded/reacting, defaults to False
    :type check: bool, optional
    :return: The heuristic distance
    :rtype: float
    """
    units = "angstrom" if angstrom else "bohr"
    frm_keys = list(automol.graph.base.ts.forming_bond_keys(tsg))
    brk_keys = list(automol.graph.base.ts.breaking_bond_keys(tsg))
    rxn_keys = frm_keys + brk_keys

    key = frozenset({key1, key2})
    if key in rxn_keys:
        symb_dct = automol.graph.base.atom_symbols(tsg)
        dist = sum(qcel.vdwradii.get(symb_dct[k], units=units) for k in key) / 2
        dist *= fdist_factor if key in frm_keys else bdist_factor
    else:
        dist = _heuristic_bond_distance(tsg, key1, key2, angstrom=angstrom, check=check)

    return dist
