""" transition state graph data structure

DEPRECATED (under construction to phase out)

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import automol.amchi.base  # !!!!
from automol import util
from automol.graph.base._0core import (
    atom_neighbor_atom_keys,
    atoms_neighbor_atom_keys,
    bond_stereo_keys,
    bond_stereo_sorted_neighbor_keys,
    has_stereo,
    local_stereo_priorities,
    set_stereo_parities,
    stereo_parities,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_graph,
    ts_reacting_atom_keys,
    ts_reacting_bond_keys,
    ts_reagents_graph_without_stereo,
    ts_reverse,
    ts_transferring_atoms,
    without_dummy_atoms,
    without_reacting_bonds,
)
from automol.graph.base._2algo import (
    rings_bond_keys,
    sorted_ring_atom_keys_from_bond_keys,
)
from automol.graph.base._3kekule import (
    rigid_planar_bond_keys,
    ts_linear_reacting_atom_keys,
    ts_reacting_electron_direction,
    vinyl_radical_atom_keys,
)
from automol.graph.base._4heur import (
    heuristic_bond_distance as _heuristic_bond_distance,
)
from automol.graph.base._6canon import (
    calculate_priorities_and_assign_stereo,
    parity_evaluator_flip_local_,
    to_local_stereo,
)
from automol.graph.base._9stereo import (
    expand_stereo,
)

# Rename TS-specific functions defined elsewhere
graph = ts_graph
forming_bond_keys = ts_forming_bond_keys
breaking_bond_keys = ts_breaking_bond_keys
reacting_bond_keys = ts_reacting_bond_keys
reacting_atom_keys = ts_reacting_atom_keys
reverse = ts_reverse
transferring_atoms = ts_transferring_atoms
reagents_graph_without_stereo = ts_reagents_graph_without_stereo
linear_reacting_atom_keys = ts_linear_reacting_atom_keys
reacting_electron_direction = ts_reacting_electron_direction


def has_reacting_ring(tsg) -> bool:
    """Does this TS graph have a ring which is involved in the reaction?

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return bool(forming_rings_atom_keys(tsg) or breaking_rings_atom_keys(tsg))


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


def reactants_graph(tsg, stereo=True, dummy=False):
    """Get the reactants from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool, optional
    :param dummy: Keep dummy atoms? default False
    :type dummy: bool, optional
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    return reagents_graph(tsg, prod=False, stereo=stereo, dummy=dummy)


def products_graph(tsg, stereo=True, dummy=False):
    """Get the products from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool, optional
    :param dummy: Keep dummy atoms? default False
    :type dummy: bool, optional
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    return reagents_graph(tsg, prod=True, stereo=stereo, dummy=dummy)


def reagents_graph(tsg, prod=False, stereo=True, dummy=False):
    """Get the reactants or products from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param prod: Do this for the products, instead of the reactants?
    :type prod: bool
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool, optional
    :param dummy: Keep dummy atoms? default False
    :type dummy: bool, optional
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    gra = ts_reagents_graph_without_stereo(tsg, prod=prod, dummy=dummy)
    if stereo and has_stereo(tsg):
        _, gra, _ = calculate_priorities_and_assign_stereo(
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
    cpar_dct = constrained_1_2_insertion_local_parities(loc_tsg0)
    loc_tsg0 = set_stereo_parities(loc_tsg0, cpar_dct)

    # Handle vinyl radical additions
    vpar_dct = vinyl_addition_local_parities(loc_tsg0)
    loc_tsg0 = set_stereo_parities(loc_tsg0, vpar_dct)

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
        assert gra == ts_reagents_graph_without_stereo(loc_tsg0)

        # Do-nothing line to prevent linting complaint
        assert ts_rev or not ts_rev

        loc_gra = ts_reagents_graph_without_stereo(loc_tsg0, keep_stereo=True)
        p0_ = par_eval_flip_(loc_gra, pri_dct)

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
    rgra = without_dummy_atoms(rcts_gra)
    pgra = without_dummy_atoms(prds_gra)
    # Allow *all* possibilities for the TS, including symmetry equivalent ones
    ste_tsgs = expand_stereo(tsg, enant=True, symeq=True)

    if has_stereo(rgra) or has_stereo(pgra):
        all_ste_tsgs = ste_tsgs
        ste_tsgs = []
        for ste_tsg in all_ste_tsgs:
            rgra_ = reactants_graph(ste_tsg, dummy=False)
            pgra_ = products_graph(ste_tsg, dummy=False)
            if rgra == rgra_ and pgra == pgra_:
                ste_tsgs.append(ste_tsg)

    return tuple(ste_tsgs)


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
    frm_keys = list(automol.graph.base.ts.forming_bond_keys(tsg))
    brk_keys = list(automol.graph.base.ts.breaking_bond_keys(tsg))
    rxn_keys = frm_keys + brk_keys

    key = frozenset({key1, key2})
    if key in rxn_keys:
        symb_dct = automol.graph.base.atom_symbols(tsg)
        symb1, symb2 = map(symb_dct.__getitem__, key)
        dist = util.heuristic.bond_distance_limit(symb1, symb2, angstrom=angstrom)
        dist *= fdist_factor if key in frm_keys else bdist_factor
    else:
        dist = _heuristic_bond_distance(tsg, key1, key2, angstrom=angstrom, check=check)

    return dist


def plane_keys(tsg, key: int, include_self: bool = True):
    """Keys used to define a plane for forming the TS geometry

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param key: The key of a bond-forming atom
    :type key: int
    :param include_self: Whether to include the key itself; defaults to `True`
    :type include_self: bool, optional
    """
    nrbs_gra = without_reacting_bonds(tsg)
    rcts_gra = reactants_graph(tsg, stereo=False)

    nkeys_rct = atom_neighbor_atom_keys(rcts_gra, key)
    nkeys_nrb = atom_neighbor_atom_keys(nrbs_gra, key)

    pkeys = {key} if include_self else set()
    pkeys |= nkeys_nrb if len(nkeys_rct) > 3 else nkeys_rct

    rp_bkeys = rigid_planar_bond_keys(rcts_gra)
    rp_bkey = next((bk for bk in rp_bkeys if key in bk), None)
    if rp_bkey is not None:
        (key_,) = rp_bkey - {key}
        pkeys |= {key_}
        pkeys |= atom_neighbor_atom_keys(rcts_gra, key_)

    return frozenset(pkeys)
