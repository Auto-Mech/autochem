""" transition state graph data structure

Forming bonds are encoded as 0.1-order bonds
Breaking bonds are encoded as 0.9-order bonds

Otherwise, this is equivalent to any other graph

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
from automol import util
import automol.amchi.base    # !!!!
from automol.graph.base._core import stereo_parities
from automol.graph.base._core import stereo_keys
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import add_bonds
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import from_ts_graph as _from_ts_graph
from automol.graph.base._core import forming_bond_keys
from automol.graph.base._core import breaking_bond_keys
from automol.graph.base._core import reacting_atoms
from automol.graph.base._core import negate_hydrogen_keys
from automol.graph.base._core import string
from automol.graph.base._algo import rings_bond_keys
from automol.graph.base._algo import sorted_ring_atom_keys_from_bond_keys
from automol.graph.base._algo import isomorphic
from automol.graph.base._canon import to_local_stereo
from automol.graph.base._canon import local_priority_dict
from automol.graph.base._canon import stereogenic_atom_keys
from automol.graph.base._canon import stereogenic_bond_keys
from automol.graph.base._canon import stereogenic_atom_keys_from_priorities
from automol.graph.base._canon import stereogenic_bond_keys_from_priorities
from automol.graph.base._canon import canonical_priorities
from automol.graph.base._kekule import vinyl_radical_atom_keys
from automol.graph.base._stereo import expand_stereo_with_priorities_and_amchis


def graph(gra, frm_bnd_keys, brk_bnd_keys):
    """ generate a transition-state graph
    """
    frm_bnd_keys = frozenset(map(frozenset, frm_bnd_keys))
    brk_bnd_keys = frozenset(map(frozenset, brk_bnd_keys))

    frm_ord_dct = {k: 0.1 for k in frm_bnd_keys}
    brk_ord_dct = {k: 0.9 for k in brk_bnd_keys}

    tsg = add_bonds(gra, frm_bnd_keys, ord_dct=frm_ord_dct, check=False)
    tsg = add_bonds(tsg, brk_bnd_keys, ord_dct=brk_ord_dct, check=False)
    return tsg


def reverse(tsg, dummies=True):
    """ reverse a transition state graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    if not dummies:
        tsg = without_dummy_atoms(tsg)

    return graph(gra=tsg,
                 frm_bnd_keys=breaking_bond_keys(tsg),
                 brk_bnd_keys=forming_bond_keys(tsg))


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
    frm_bnd_keys = forming_bond_keys(tsg)
    frm_rngs_bnd_keys = tuple(
        bks for bks in rings_bond_keys(tsg, ts_graph=True)
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
    brk_bnd_keys = breaking_bond_keys(tsg)
    brk_rngs_bnd_keys = tuple(
        bks for bks in rings_bond_keys(tsg, ts_graph=True)
        if brk_bnd_keys & bks)
    return brk_rngs_bnd_keys


def reactants_graph(tsg):
    """ get a graph of the reactants from a transition state graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    gra = _from_ts_graph(tsg)
    return gra


def products_graph(tsg):
    """ get a graph of the products from a transition state graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    return reactants_graph(reverse(tsg))


def fleeting_stereogenic_atom_keys(tsg, enant=True):
    """ Identify fleeting stereogenic atoms in a TS structure

    Setting `enant=False` generates only *dia*stereogenic atom keys.  The atom
    is deemed diastereogenic if its inversion results in a different
    diastereomer. For now, this is assumed to be the case unless it is the only
    chiral atom in the TS, in which case we assume inversion generates an
    enantiomer.

    A stereocenter is 'fleeting' if it occurs only in the TS, not in the
    reactants or products.

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param enant: Include fleeting stereogenic atoms that generate enantiomers
        rather than diastereomers upon inversion?
    :type enant: bool
    """
    tsg = without_dummy_atoms(tsg)
    pri_dct = canonical_priorities(tsg, backbone_only=False, ts_graph=True)

    tsg_ste_atm_keys = stereogenic_atom_keys_from_priorities(
        tsg, pri_dct=pri_dct, assigned=True, ts_graph=True)
    rct_ste_atm_keys = stereogenic_atom_keys(reactants_graph(tsg),
                                             assigned=True)
    prd_ste_atm_keys = stereogenic_atom_keys(products_graph(tsg),
                                             assigned=True)
    ste_atm_keys = (tsg_ste_atm_keys - rct_ste_atm_keys) - prd_ste_atm_keys

    # Handle the case where we only want diastereogenic atoms.
    if not enant:
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
    pri_dct = canonical_priorities(tsg, backbone_only=False, ts_graph=True)

    tsg_ste_bnd_keys = stereogenic_bond_keys_from_priorities(
        tsg, pri_dct=pri_dct, assigned=True, ts_graph=True)
    rct_ste_bnd_keys = stereogenic_bond_keys(reactants_graph(tsg),
                                             assigned=True)
    prd_ste_bnd_keys = stereogenic_bond_keys(products_graph(tsg),
                                             assigned=True)
    ste_bnd_keys = (tsg_ste_bnd_keys - rct_ste_bnd_keys) - prd_ste_bnd_keys
    return ste_bnd_keys


def fleeting_stereogenic_keys(tsg, enant=True):
    """ Identify fleeting stereogenic atoms and bonds in a TS structure

    A stereocenter is 'fleeting' if it occurs only in the TS, not in the
    reactants or products.

    Setting `enant=False` generates only *dia*stereogenic keys.  The atom/bond
    is deemed diastereogenic if its inversion results in a different
    diastereomer. For bonds, this is always the case. For atoms, this is
    assumed to be the case unless there are no other chiral atoms, in which
    case we assume inversion generates an enantiomer.

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param enant: Include fleeting stereogenic atoms that generate enantiomers
        rather than diastereomers upon inversion?
    :type enant: bool
    """
    ste_atm_keys = fleeting_stereogenic_atom_keys(tsg, enant=enant)
    ste_bnd_keys = fleeting_stereogenic_bond_keys(tsg)
    ste_keys = ste_atm_keys | ste_bnd_keys
    return ste_keys


# def fleeting_stereosite_reacting_neighbors(tsg, enant=True):
#     """ Reacting neighbor atoms at fleeting stereosites

#     :param tsg: TS graph
#     :type tsg: automol graph data structure
#     :param enant: Include fleeting stereogenic atoms that generate
#          enantiomers rather than diastereomers upon inversion?
#     :type enant: bool
#     :returns: A dictionary keyed by fleeting stereogenic keys, with values of
#         the specific neighbors that are reacting ()
#     """
#     rxn_atm_keys = reacting_atoms(tsg)

#     def _sorted_nkeys(nkeys):
#         """ Sort neighboring keys based on proximity to the reaction site
#         """

#     nkeys_dct = atoms_neighbor_atom_keys(tsg)
#     ste_atm_keys = fleeting_stereogenic_atom_keys(tsg, enant=enant)
#     ste_bnd_keys = fleeting_stereogenic_bond_keys(tsg)


def are_energetically_equivalent(tsg1, tsg2):
    """ Are these Reaction objects energetically equivalent?

    Requires two TS graphs that are exactly aligned, having identical reactant
    graphs (including their keys) and differing only in the breaking/forming
    bonds.

    They are energetically equivalent if
    (a.) their reactant/product graphs are identical (same indexing)
    (b.) their TS graphs are isomorphic, and
    (c.) they DO NOT differ in their breaking/forming bonds at a fleeting
    stereocenter in a way that makes them diastereomers. (If the difference
    makes them enantiomers, they are still energetically equivalent.)

    A stereocenter is 'fleeting' if it occurs only in the TS, not in the
    reactants or products.

    :param tsg1: TS graph
    :type tsg1: automol graph data structure
    :param tsg2: TS graph for comparison
    :type tsg2: automol graph data structure
    :returns: `True` if they are, `False` if they aren't
    """
    assert reactants_graph(tsg1) == reactants_graph(tsg2), (
        f"This function assumes these TS graphs are exactly aligned, apart\n"
        f"from their breaking/forming bonds, but they aren't:"
        f"{string(tsg1)}\n----\n{string(tsg2)}"
    )

    ret = isomorphic(tsg1, tsg2, stereo=True)
    # If they are isomorphic, determine if they are fleeting diastereomers
    if ret:
        pass

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
    reac_keys = reacting_atoms(ftsg_loc)
    fste_keys = stereo_keys(ftsg_loc)
    rste_keys = stereo_keys(rtsg_loc)
    cons_keys = list(fste_keys & rste_keys)

    rcts_gra = negate_hydrogen_keys(reactants_graph(ftsg_loc))
    prds_gra = negate_hydrogen_keys(reactants_graph(rtsg_loc))

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
    brk_bkeys = breaking_bond_keys(ftsg_loc)
    ngb_keys_dct = atoms_neighbor_atom_keys(ftsg_loc)
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
                nk1s = sorted(ngb_keys_dct[key1], key=loc_pri_dct.__getitem__)
                nk2s = sorted(ngb_keys_dct[key2], key=loc_pri_dct.__getitem__)
                srt_nk1s = util.move_items_to_front(nk1s, [key2, key0])
                srt_nk2s = util.move_items_to_front(nk2s, [key1, key3])
                p_1 = util.is_even_permutation(nk1s, srt_nk1s)
                p_2 = util.is_even_permutation(nk2s, srt_nk2s)

                # p_b12 = (p_a1 XNOR p_1) XNOR (p_a2 XNOR p_2)
                p_b12 = not ((not (p_a1 ^ p_1)) ^ (not (p_a2 ^ p_2)))
                satisfies &= (p_b12 == p_b12_actual)

    return satisfies
