""" transition state graph data structure

DEPRECATED (under construction to phase out)

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
import numpy
from automol import util
import automol.amchi.base    # !!!!
from automol.graph.base._core import ts_graph
from automol.graph.base._core import ts_breaking_bond_keys
from automol.graph.base._core import ts_forming_bond_keys
from automol.graph.base._core import ts_reacting_bond_keys
from automol.graph.base._core import ts_reacting_atom_keys
from automol.graph.base._core import ts_reverse
from automol.graph.base._core import ts_reagents_graph_without_stereo
from automol.graph.base._core import has_stereo
from automol.graph.base._core import stereo_parities
from automol.graph.base._core import set_stereo_parities
from automol.graph.base._core import stereo_keys
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import string
from automol.graph.base._core import relabel
from automol.graph.base._core import nonbackbone_hydrogen_keys
from automol.graph.base._core import bond_stereo_sorted_neighbor_keys
from automol.graph.base._algo import rings_bond_keys
from automol.graph.base._algo import isomorphic
from automol.graph.base._algo import shortest_path_between_groups
from automol.graph.base._algo import sorted_ring_atom_keys_from_bond_keys
from automol.graph.base._canon import to_local_stereo
from automol.graph.base._canon import local_priority_dict
from automol.graph.base._canon import stereogenic_atom_keys
from automol.graph.base._canon import stereogenic_bond_keys
from automol.graph.base._canon import stereogenic_atom_keys_from_priorities
from automol.graph.base._canon import stereogenic_bond_keys_from_priorities
from automol.graph.base._canon import canonical_priorities
from automol.graph.base._canon import calculate_priorities_and_assign_stereo
from automol.graph.base._canon import parity_evaluator_flip_local_
from automol.graph.base._kekule import vinyl_radical_atom_keys
from automol.graph.base._kekule import rigid_planar_bond_keys
from automol.graph.base._stereo import expand_stereo
from automol.graph.base._stereo import expand_stereo_with_priorities_and_amchis

# Rename TS-specific functions defined elsewhere
graph = ts_graph
forming_bond_keys = ts_forming_bond_keys
breaking_bond_keys = ts_breaking_bond_keys
reacting_bond_keys = ts_reacting_bond_keys
reacting_atom_keys = ts_reacting_atom_keys
reverse = ts_reverse
reagents_graph_without_stereo = ts_reagents_graph_without_stereo


def forming_rings_atom_keys(tsg):
    """ get the atom keys to rings forming in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    frm_rngs_bnd_keys = forming_rings_bond_keys(tsg)
    frm_rngs_atm_keys = tuple(map(sorted_ring_atom_keys_from_bond_keys,
                                  frm_rngs_bnd_keys))
    return frm_rngs_atm_keys


def forming_rings_bond_keys(tsg):
    """ get the bond keys to rings forming in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    frm_bnd_keys = ts_forming_bond_keys(tsg)
    frm_rngs_bnd_keys = tuple(
        bks for bks in rings_bond_keys(tsg, ts_=True)
        if frm_bnd_keys & bks)
    return frm_rngs_bnd_keys


def breaking_rings_atom_keys(tsg):
    """ get the atom keys to rings breaking in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    brk_rngs_bnd_keys = breaking_rings_bond_keys(tsg)
    brk_rngs_atm_keys = tuple(map(sorted_ring_atom_keys_from_bond_keys,
                                  brk_rngs_bnd_keys))
    return brk_rngs_atm_keys


def breaking_rings_bond_keys(tsg):
    """ get the bond keys to rings breaking in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    brk_bnd_keys = ts_breaking_bond_keys(tsg)
    brk_rngs_bnd_keys = tuple(
        bks for bks in rings_bond_keys(tsg, ts_=True)
        if brk_bnd_keys & bks)
    return brk_rngs_bnd_keys


def reactants_graph(tsg, stereo=True):
    """ Get the reactants from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    return reagents_graph(tsg, prod=False, stereo=stereo)


def products_graph(tsg, stereo=True):
    """ Get the products from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    return reagents_graph(tsg, prod=True, stereo=stereo)


def reagents_graph(tsg, prod=False, stereo=True):
    """ Get the reactants or products from a TS graph

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
            gra, backbone_only=False, break_ties=False,
            par_eval_=parity_evaluator_reagents_from_ts_(tsg, prod=prod))
    return gra


def parity_evaluator_reagents_from_ts_(tsg, prod=False):
    r""" Determines reactant or product stereochemistry from a TS graph

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
        """ Parity evaluator based on current priorities

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
                    f"yet implemented:\n{tsg}")

            return par

        return _parity

    return _evaluator


def expand_stereo_for_reaction(tsg, rcts_gra, prds_gra):
    """ Expand TS graph stereo that is consistent with reactants and products

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
    assert not any(p is None for p in fpars0 + rpars0), (
        "Corrupted stereo assignments: {rcts_gra}\n{prds_gra}")

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
    r""" Calculates the local parity of the reactant or product of a
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
        stereo_parities(loc_tsg), lambda x: x is not None)
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
                key0, = next(k for k in frm_bkeys if key1 in k) - {key1}
                key3, = next(k for k in frm_bkeys if key2 in k) - {key2}
                keys = {key0, key1, key2, key3}
                # C. The atoms are part of a forming ring
                if any(keys & rkeys for rkeys in rkeys_lst):
                    # This is a constrained bond!
                    # i. Read out the relevant atom and bond parities
                    p_a1 = loc_par_dct[key1]
                    p_a2 = loc_par_dct[key2]

                    # ii. Calculate the local bond parity
                    nk1s = sorted(
                        nkeys_dct[key1], key=loc_pri_dct.__getitem__)
                    nk2s = sorted(
                        nkeys_dct[key2], key=loc_pri_dct.__getitem__)
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
        stereo_parities(loc_tsg), lambda x: x is not None)
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
                loc_tsg, *akeys, pri_dct=loc_pri_dct)
            gnkeys_pair = bond_stereo_sorted_neighbor_keys(
                gra, *akeys, pri_dct=loc_pri_dct)
            for (akey, tnkeys, gnkeys) in zip(akeys, tnkeys_pair, gnkeys_pair):
                # Compare sorted neighbors for mismatch
                if tnkeys != gnkeys:
                    assert akey in vin_keys, (
                        f"Neighbor mismatch at {akey} is not due to "
                        f"vinyl addition:\n{loc_tsg}")
                    # If the maximum priority neighbors don't match, flip
                    # the parity
                    if tnkeys[-1] != gnkeys[-1]:
                        # Flip the parity
                        par = not par

                        # There should only be one, so quit if we find it
                        return bkey, par, prod

    return (None, None, None)


# vvvvvvvvvvvvvvvvvvvvvvvvv DEPRECTATED vvvvvvvvvvvvvvvvvv
def reaction_stereo_satisfies_elimination_constraint(ftsg_loc, rtsg_loc):
    r""" Check whether a reaction satisfies the elimination stereo constraint

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
    """
    floc_par_dct = stereo_parities(ftsg_loc)
    rloc_par_dct = stereo_parities(rtsg_loc)

    loc_pri_dct = local_priority_dict(reactants_graph(ftsg_loc))

    ste_akeys = atom_stereo_keys(ftsg_loc)
    brk_bkeys = ts_breaking_bond_keys(ftsg_loc)
    nkeys_dct = atoms_neighbor_atom_keys(ftsg_loc)
    rng_akeys_lst = list(map(set, forming_rings_atom_keys(ftsg_loc)))

    satisfies = True
    bkeys = bond_stereo_keys(rtsg_loc)
    for bkey in bkeys:
        # Conditions:
        # A. Both atoms in the bond were stereogenic in the reactants
        if all(k in ste_akeys for k in bkey):
            key1, key2 = bkey
            key0, = next(k for k in brk_bkeys if key1 in k) - {key1}
            key3, = next(k for k in brk_bkeys if key2 in k) - {key2}
            keys = {key0, key1, key2, key3}
            # B. The leaving groups are joined together, forming a TS ring
            if any(keys & ks for ks in rng_akeys_lst):
                # This is a constrained bond!
                # i. Read out the relevant atom and bond parities
                p_a1 = floc_par_dct[key1]
                p_a2 = floc_par_dct[key2]
                p_b12_actual = rloc_par_dct[frozenset({key1, key2})]

                # ii. Calculate what the bond parity should be
                nk1s = sorted(nkeys_dct[key1], key=loc_pri_dct.__getitem__)
                nk2s = sorted(nkeys_dct[key2], key=loc_pri_dct.__getitem__)
                srt_nk1s = util.move_items_to_front(nk1s, [key2, key0])
                srt_nk2s = util.move_items_to_front(nk2s, [key1, key3])
                p_1 = util.is_even_permutation(nk1s, srt_nk1s)
                p_2 = util.is_even_permutation(nk2s, srt_nk2s)

                # p_b12 = (p_a1 XNOR p_1) XNOR (p_a2 XNOR p_2)
                p_b12 = not ((not (p_a1 ^ p_1)) ^ (not (p_a2 ^ p_2)))
                satisfies &= (p_b12 == p_b12_actual)

    return satisfies


def negate_nonbackbone_hydrogen_keys(gra):
    """ Flip the signs of hydrogen keys

    :param gra: molecular graph
    :type gra: automol graph data structure
    :return: molecular graph with hydrogen keys negated
    :rtype: automol graph data structure
    """
    gra = reactants_graph(gra)
    hyd_keys = nonbackbone_hydrogen_keys(gra)
    atm_key_dct = {k: -abs(k) for k in hyd_keys}
    return relabel(gra, atm_key_dct)


def fleeting_stereogenic_atom_keys(tsg, ts_enant=True):
    """ Identify fleeting stereogenic atoms in a TS structure

        A stereocenter is 'fleeting' if it occurs only in the TS, not in the
        reactants or products.

        Setting `ts_enant=False` generates only *dia*stereogenic keys.  The
        atom/bond is deemed diastereogenic if its inversion results in a
        different diastereomer. For bonds, this is always the case. For atoms,
        this is assumed to be the case unless there are no other chiral atoms,
        in which case we assume inversion generates an enantiomer.

        :param tsg: TS graph
        :type tsg: automol graph data structure
        :param ts_enant: Include fleeting enantiomer stereo sites?
        :type ts_enant: bool
    """
    tsg = without_dummy_atoms(tsg)
    pri_dct = canonical_priorities(tsg, backbone_only=False)

    tsg_ste_atm_keys = stereogenic_atom_keys_from_priorities(
        tsg, pri_dct=pri_dct, assigned=True)
    rct_ste_atm_keys = stereogenic_atom_keys(reactants_graph(tsg),
                                             assigned=True)
    prd_ste_atm_keys = stereogenic_atom_keys(products_graph(tsg),
                                             assigned=True)
    ste_atm_keys = (tsg_ste_atm_keys - rct_ste_atm_keys) - prd_ste_atm_keys

    # Handle the case where we only want diastereogenic atoms.
    if not ts_enant:
        # If there was only one fleeting atom stereocenter, inverting it will
        # create an enantiomer rather than a diastereomer. Exclude it.
        if tsg_ste_atm_keys & ste_atm_keys and len(tsg_ste_atm_keys) == 1:
            ste_atm_keys -= tsg_ste_atm_keys

    return ste_atm_keys


def fleeting_stereogenic_bond_keys(tsg):
    """ Identify fleeting stereogenic bonds in a TS structure

        A stereocenter is 'fleeting' if it occurs only in the TS, not in the
        reactants or products.

        :param tsg: TS graph
        :type tsg: automol graph data structure
    """
    tsg = without_dummy_atoms(tsg)
    pri_dct = canonical_priorities(tsg, backbone_only=False)

    tsg_ste_bnd_keys = stereogenic_bond_keys_from_priorities(
        tsg, pri_dct=pri_dct, assigned=True)
    rct_ste_bnd_keys = stereogenic_bond_keys(reactants_graph(tsg),
                                             assigned=True)
    prd_ste_bnd_keys = stereogenic_bond_keys(products_graph(tsg),
                                             assigned=True)
    ste_bnd_keys = (tsg_ste_bnd_keys - rct_ste_bnd_keys) - prd_ste_bnd_keys
    return ste_bnd_keys


def fleeting_stereogenic_keys(tsg, ts_enant=True):
    """ Identify fleeting stereogenic atoms and bonds in a TS structure

        A stereocenter is 'fleeting' if it occurs only in the TS, not in the
        reactants or products.

        Setting `ts_enant=False` generates only *dia*stereogenic keys.  The
        atom/bond is deemed diastereogenic if its inversion results in a
        different diastereomer. For bonds, this is always the case. For atoms,
        this is assumed to be the case unless there are no other chiral atoms,
        in which case we assume inversion generates an enantiomer.

        :param tsg: TS graph
        :type tsg: automol graph data structure
        :param ts_enant: Include fleeting enantiomer stereo sites?
        :type ts_enant: bool
    """
    ste_atm_keys = fleeting_stereogenic_atom_keys(tsg, ts_enant=ts_enant)
    ste_bnd_keys = fleeting_stereogenic_bond_keys(tsg)
    ste_keys = ste_atm_keys | ste_bnd_keys
    return ste_keys


def fleeting_stereosite_sorted_neighbors(tsg, ts_enant=True):
    """ Neighbor atoms at fleeting stereosites, sorted by proximity to the
        reaction site

        For TS graphs which are aligned apart from their breaking and forming
        bond keys, this graph will be different for different TS diastereomers.

        :param tsg: TS graph
        :type tsg: automol graph data structure
        :param ts_enant: Include fleeting enantiomer stereo sites?
        :type ts_enant: bool
        :returns: A dictionary keyed by fleeting stereogenic keys, with values
        of the specific neighbors that are reacting ()
    """
    rxn_atm_keys = ts_reacting_atom_keys(tsg)
    nkeys_dct = atoms_neighbor_atom_keys(tsg)

    def _distance_from_reaction_site(key):
        path = shortest_path_between_groups(tsg, {key}, rxn_atm_keys)
        dist = len(path) if path is not None else numpy.inf
        metric = (dist, key)
        return metric

    def _sorted_nkeys(key, excl_key=None):
        """ Sort neighboring keys based on proximity to the reaction site
        """
        nkeys = nkeys_dct[key] - {excl_key}
        return tuple(sorted(nkeys, key=_distance_from_reaction_site))

    # Add atom stereo sites
    ste_atm_keys = fleeting_stereogenic_atom_keys(tsg, ts_enant=ts_enant)
    srt_ste_nkey_dct = {k: _sorted_nkeys(k) for k in ste_atm_keys}

    # Add bond stereo sites
    ste_bnd_keys = list(map(sorted, fleeting_stereogenic_bond_keys(tsg)))
    srt_ste_nkey_dct.update(
        {(k2, k1): _sorted_nkeys(k1, k2)
         for bk in ste_bnd_keys for k1, k2 in itertools.permutations(bk)}
    )

    return srt_ste_nkey_dct


def are_equivalent(tsg1, tsg2, ts_stereo=True, ts_enant=False):
    """ Are these TS graphs energetically equivalent?

    Requires two TS graphs that are exactly aligned, having identical reactant
    graphs (including their keys) and differing only in the breaking/forming
    bonds. The underlying assumption is that both TS graphs refer to the same
    set of initial initial reactant geometries.

    By default, they are deemed equivalent if they have the same energy, which
    occurs when:
    (a.) their reactant/product graphs are identical (same indexing)
    (b.) their TS graphs are isomorphic, and
    (c.) the neighboring atoms at non-enantiomeric fleeting reaction sites are
    equidistant from the reaction site.

    A stereocenter is 'fleeting' if it occurs only in the TS, not in the
    reactants or products.

    Note that differences at "enantiomeric" fleeting reaction sites are still
    considered energetically equivalent, since the TS structures are mirror
    images. The only time a reaction site is treated as enantiomeric is when
    there are no other fleeting or non-fleeting atom stereosites in the TS.

    This assumption could in principle break down in cases where multiple
    fleeting atom reaction sites are simultaneously formed by the reaction and
    are all simultaneously inverted in the second TS relative to the first. I
    can't currently think of any cases where this would happen. Initially, I
    thought Diels-Alder reactions would be the exception, but in that case the
    stereo sites are not fleeting, so the different possibilities are captured
    by the stereochemistry of the reactants and products. Diels-Alder reactions
    with some topological symmetry may be exceptions.

    :param tsg1: TS graph
    :type tsg1: automol graph data structure
    :param tsg2: TS graph for comparison
    :type tsg2: automol graph data structure
    :param ts_stereo: Treat fleeting TS stereoisomers as distinct TSs?
    :type ts_stereo: bool
    :param ts_enant: Treat fleeting TS enantiomers as distinct TSs?
    :type ts_enant: bool
    :returns: `True` if they are, `False` if they aren't
    """
    assert reactants_graph(tsg1) == reactants_graph(tsg2), (
        f"This function assumes these TS graphs are exactly aligned, apart\n"
        f"from their breaking/forming bonds, but they aren't:"
        f"{string(tsg1)}\n----\n{string(tsg2)}"
    )

    ret = isomorphic(tsg1, tsg2, stereo=True)
    # If they are isomorphic, determine if they are fleeting diastereomers.
    # Since they are fully aligned (and, in context, are formed from the same
    # reactant geometries), this is assumed to be the case if the neighboring
    # atoms at non-enantiomeric reaction sites have different distances to the
    # reaction site (see docstring).
    if ts_stereo and ret:
        srt_ste_nkey_dct1 = fleeting_stereosite_sorted_neighbors(
            tsg1, ts_enant=ts_enant)
        srt_ste_nkey_dct2 = fleeting_stereosite_sorted_neighbors(
            tsg2, ts_enant=ts_enant)
        # For any fleeting stereosites that the two graphs have in common,
        # check whether their neighoring atoms have the same relationship to
        # the reaction site.
        for key in srt_ste_nkey_dct1:
            if key in srt_ste_nkey_dct2:
                ret &= srt_ste_nkey_dct1[key] == srt_ste_nkey_dct2[key]

    return ret


def expand_reaction_stereo(tsg, enant=True, symeq=False, const=True,
                           log=False):
    """ Obtain all possible stereoisomer combinations for a reaction, encoding
        reactant and product stereo assignments in forward and reverse TS
        graphs, respectively

        :param tsg: A TS graph; If it has stereo assignments, they will be
            assumed to be for the reactants, and all compatible products
            assignments will be exapnded
        :type tsg: automol graph data structure
        :param enant: Include all enantiomers, or only canonical ones?
        :type enant: bool
        :param symeq: Include symmetrically equivalent stereoisomers?
        :type symeq: bool
        :param const: Constrain bond stereo based on reactant atom stereo?
        :type const: bool
        :param log: Print information to the screen?
        :type log: bool
        :returns: a series of pairs of forward and reverse graphs containing
            mutually compatible stereo assignments
    """
    rxns = []

    # 1. Expand all possible reactant and product assignments
    ftsgs, fpri_dcts, fchis = zip(
        *expand_stereo_with_priorities_and_amchis(tsg))
    rtsgs, rpri_dcts, rchis = zip(
        *expand_stereo_with_priorities_and_amchis(reverse(tsg)))

    seen_rxn_chis = []
    # 2. Loop over forward TS graphs (reactant assignments)
    for ftsg, fpri_dct, fchi in zip(ftsgs, fpri_dcts, fchis):
        # 3. Convert the forward TS graph to local stereo.
        ftsg_loc = to_local_stereo(ftsg, pri_dct=fpri_dct)

        # 4. Loop over reverse TS graphs (product assignments)
        for rtsg, rpri_dct, rchi in zip(rtsgs, rpri_dcts, rchis):
            # 5. Convert the reverse TS graph to local stereo.
            rtsg_loc = to_local_stereo(rtsg, pri_dct=rpri_dct)

            # 6. Check if the local stereo assignments for reactants and
            # products are mutually compatible.
            if reaction_stereo_is_physical(ftsg_loc, rtsg_loc, const=const):

                # 7. Check if the reaction is canonical.
                if (enant or
                        automol.amchi.base.is_canonical_enantiomer_reaction(
                            fchi, rchi)):

                    # 8. Check if the reaction has been seen before (symmetry)
                    if symeq or (fchi, rchi) not in seen_rxn_chis:
                        rxns.append((ftsg, rtsg))
                        seen_rxn_chis.append((fchi, rchi))

                        if log:
                            print(f'{fchi} =>\n{rchi}\n')

    return tuple(rxns)


def reaction_stereo_is_physical(ftsg_loc, rtsg_loc, const=True):
    """ Does this pair of forward and reverse reactions have compatible stereo?

        :param ftsg_loc: a forward TS graph, with local stereo assignments for
            the reactants
        :type ftsg_loc: automol graph data structure
        :param rtsg_loc: a reverse TS graph, with local stereo assignments for
            the products
        :type rtsg_loc: automol graph data structure
        :param const: Constrain bond stereo based on reactant atom stereo?
        :type const: bool
        :rtype: bool
    """
    # 1. Check conserved stereo sites
    reac_keys = ts_reacting_atom_keys(ftsg_loc)
    fste_keys = stereo_keys(ftsg_loc)
    rste_keys = stereo_keys(rtsg_loc)
    cons_keys = list(fste_keys & rste_keys)

    rcts_gra = negate_nonbackbone_hydrogen_keys(reactants_graph(ftsg_loc))
    prds_gra = negate_nonbackbone_hydrogen_keys(reactants_graph(rtsg_loc))

    fnkeys_dct = atoms_neighbor_atom_keys(rcts_gra)
    rnkeys_dct = atoms_neighbor_atom_keys(prds_gra)
    fvin_keys = vinyl_radical_atom_keys(rcts_gra)
    rvin_keys = vinyl_radical_atom_keys(prds_gra)

    floc_par_dct = stereo_parities(ftsg_loc)
    rloc_par_dct = stereo_parities(rtsg_loc)

    is_consistent = True

    for cons_key in cons_keys:
        # Conserved atom keys should not be involved in the reaction
        if not isinstance(cons_key, frozenset):
            assert cons_key not in reac_keys, (
                f"Assumption fails! Conserved atom stereo site {cons_key} "
                f"is a reacting atom: {reac_keys}.")

            is_consistent &= (floc_par_dct[cons_key] ==
                              rloc_par_dct[cons_key])
        # Conserved bond keys may be involved in the reaction, but we must make
        # sure their local parities don't change
        else:
            key1, key2 = sorted(cons_key, reverse=True,
                                key=lambda k: fnkeys_dct[k] == rnkeys_dct[k])
            fnkey1s = sorted(fnkeys_dct[key1] - {key2})
            fnkey2s = sorted(fnkeys_dct[key2] - {key1})
            rnkey1s = sorted(rnkeys_dct[key1] - {key2})
            rnkey2s = sorted(rnkeys_dct[key2] - {key1})

            if (max(fnkey1s) == max(rnkey1s) and
                    max(fnkey2s) == max(rnkey2s)):
                is_consistent &= (floc_par_dct[cons_key] ==
                                  rloc_par_dct[cons_key])
            else:
                assert ((cons_key & fvin_keys or cons_key & rvin_keys) and
                        (fnkey1s == rnkey1s and fnkey2s[0] == rnkey2s[0])), (
                    f"Assumption fails! Conserved non-vinyl bond stereo site "
                    f"{cons_key} will have its local parity altered."
                    f"\nForward neighbors: {fnkey1s} / {fnkey2s}"
                    f"\nReverse neighbors: {rnkey1s} / {rnkey2s}")

                is_consistent &= (floc_par_dct[cons_key] !=
                                  rloc_par_dct[cons_key])

    # 2. Check for insertion/elimination stereo constraint
    if const:
        # A. Elimination
        is_consistent &= reaction_stereo_satisfies_elimination_constraint(
            ftsg_loc, rtsg_loc)
        # B. Insertion
        is_consistent &= reaction_stereo_satisfies_elimination_constraint(
            rtsg_loc, ftsg_loc)

    return is_consistent
