""" transition state graph data structure

Forming bonds are encoded as 0.1-order bonds
Breaking bonds are encoded as 0.9-order bonds

Otherwise, this is equivalent to any other graph

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
from automol import util
from automol.graph.base._core import bond_orders
from automol.graph.base._core import stereo_parities
from automol.graph.base._core import stereo_keys
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import add_bonds
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import has_stereo
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import from_ts_graph as _from_ts_graph
from automol.graph.base._algo import rings_bond_keys
from automol.graph.base._algo import sorted_ring_atom_keys_from_bond_keys
from automol.graph.base._canon import to_local_stereo as _to_local_stereo
from automol.graph.base._canon import from_local_stereo as _from_local_stereo
from automol.graph.base._canon import stereogenic_keys
from automol.graph.base._canon import local_priority_dict
from automol.graph.base._stereo import expand_stereo as _expand_stereo
from automol.graph.base._stereo import expand_stereo_with_priorities_and_amchis
from automol.graph.base._stereo import canonical_assignment_representation


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


def forming_bond_keys(tsg):
    """ get the forming bonds from a transition state graph
    """
    ord_dct = bond_orders(tsg)
    frm_bnd_keys = [k for k, o in ord_dct.items() if round(o, 1) == 0.1]
    return frozenset(map(frozenset, frm_bnd_keys))


def breaking_bond_keys(tsg):
    """ get the forming bonds from a transition state graph
    """
    ord_dct = bond_orders(tsg)
    brk_bnd_keys = [k for k, o in ord_dct.items() if round(o, 1) == 0.9]
    return frozenset(map(frozenset, brk_bnd_keys))


def reacting_atoms(tsg):
    """ get all of the atoms involved in the reaction
    """
    bnd_keys = forming_bond_keys(tsg) | breaking_bond_keys(tsg)
    atm_keys = frozenset(itertools.chain(*bnd_keys))
    return atm_keys


def reverse(tsg, dummies=True):
    """ reverse a transition state graph
    """
    if not dummies:
        tsg = without_dummy_atoms(tsg)

    return graph(gra=tsg,
                 frm_bnd_keys=breaking_bond_keys(tsg),
                 brk_bnd_keys=forming_bond_keys(tsg))


def forming_rings_atom_keys(tsg):
    """ get the atom keys to rings forming in the TS graph
    """
    frm_rngs_bnd_keys = forming_rings_bond_keys(tsg)
    frm_rngs_atm_keys = tuple(map(sorted_ring_atom_keys_from_bond_keys,
                                  frm_rngs_bnd_keys))
    return frm_rngs_atm_keys


def forming_rings_bond_keys(tsg):
    """ get the bond keys to rings forming in the TS graph
    """
    frm_bnd_keys = forming_bond_keys(tsg)
    frm_rngs_bnd_keys = tuple(bks for bks in rings_bond_keys(tsg)
                              if frm_bnd_keys & bks)
    return frm_rngs_bnd_keys


def breaking_rings_atom_keys(tsg):
    """ get the atom keys to rings breaking in the TS graph
    """
    brk_rngs_bnd_keys = breaking_rings_bond_keys(tsg)
    brk_rngs_atm_keys = tuple(map(sorted_ring_atom_keys_from_bond_keys,
                                  brk_rngs_bnd_keys))
    return brk_rngs_atm_keys


def breaking_rings_bond_keys(tsg):
    """ get the bond keys to rings breaking in the TS graph
    """
    brk_bnd_keys = breaking_bond_keys(tsg)
    brk_rngs_bnd_keys = tuple(bks for bks in rings_bond_keys(tsg)
                              if brk_bnd_keys & bks)
    return brk_rngs_bnd_keys


def reactants_graph(tsg):
    """ get a graph of the reactants from a transition state graph
    """
    gra = _from_ts_graph(tsg)
    return gra


def products_graph(tsg):
    """ get a graph of the products from a transition state graph
    """
    return reactants_graph(reverse(tsg))


def to_local_stereo(tsg):
    """ Convert a TS graph to local stereo assignments, where parities are
        defined relative to the ordering of indices rather than the canonical
        stereo priority.

    :param tsg: a TS graph with canonical stereo assignments
    :returns: a TS graph with local stereo assignments
    """
    rgra = reactants_graph(tsg)
    frm_bnd_keys = forming_bond_keys(tsg)
    brk_bnd_keys = breaking_bond_keys(tsg)

    rgra = _to_local_stereo(rgra)
    loc_tsg = graph(rgra, frm_bnd_keys, brk_bnd_keys)
    return loc_tsg


def from_local_stereo(loc_tsg):
    """ Convert a TS graph from local stereo assignments back to canonical
    stereo assignments, where parities are independent of atom ordering.

    :param loc_tsg: a TS graph with local stereo assignments
    :returns: a TS graph with canonical stereo assignments
    """
    rgra = reactants_graph(loc_tsg)
    frm_bnd_keys = forming_bond_keys(loc_tsg)
    brk_bnd_keys = breaking_bond_keys(loc_tsg)

    rgra = _from_local_stereo(rgra)
    tsg = graph(rgra, frm_bnd_keys, brk_bnd_keys)
    return tsg


def expand_stereo(tsg, symeq=False):
    """ Obtain all possible stereoisomers of a TS graph, ignoring its assignments

        :param tsg: A TS graph, without stereo
        :type tsg: automol graph data structure
        :param symeq: Include symmetrically equivalent stereoisomers?
        :type symeq: bool
        :returns: a series of molecular graphs for the stereoisomers
    """
    return _expand_stereo(tsg, symeq=symeq)


def expand_reaction_stereo(tsg, symeq=False, const=True):
    """ Obtain all possible stereoisomer combinations for a reaction, encoding
        reactant and product stereo assignments in forward and reverse TS
        graphs, respectively

        :param tsg: A TS graph; If it has stereo assignments, they will be
            assumed to be for the reactants, and all compatible products
            assignments will be exapnded
        :type tsg: automol graph data structure
        :param symeq: Include symmetrically equivalent stereoisomers?
        :type symeq: bool
        :param const: Constrain bond stereo based on reactant atom stereo?
        :type const: bool
        :returns: a series of pairs of forward and reverse graphs containing
            mutually compatible stereo assignments
    """
    # If the input graph has stereo, this fixes the stereo of the reactants.
    # A limited expansion over compatible products will be performed.
    if has_stereo(tsg):
        assert not stereogenic_keys(tsg)
        ftsgs = [tsg]
    # Otherwise, the stereo of the reactants is unspecified. Expand all unique
    # stereo assignments for the reactants.
    # A complete expansion over all unique reactions will be performed.
    else:
        ftsgs = expand_stereo(tsg, symeq=False)

    # 1. Expand all possible reverse (product) stereo assignments, along with
    # their priorities (for symmetry filtering)
    rtsgs, rpri_dcts = zip(*expand_stereo_with_priorities(reverse(tsg)))

    reacs = []
    # 2. Loop over forward TS graphs (reactant stereo assignments)
    for ftsg in ftsgs:
        # 3. Convert the forward TS graph to local stereo.
        ftsg_loc = to_local_stereo(ftsg)

        seen_rreps = []
        # 4. Loop over reverse TS graphs (product stereo assignments)
        for rtsg, rpri_dct in zip(rtsgs, rpri_dcts):
            # 5. Convert the reverse TS graph to local stereo.
            rtsg_loc = to_local_stereo(rtsg)

            # 6. Check if the local stereo assignments match for conserved
            # stereo centers.
            if reaction_stereo_is_physical(ftsg_loc, rtsg_loc, const=const):
                # 7. Check if this is a unique assignment for the reverse
                # reaction. If so, add the forward and reverse graphs to the
                # list of reacs.
                rrep = canonical_assignment_representation(rtsg, rpri_dct)
                if symeq or rrep not in seen_rreps:
                    reacs.append((ftsg, rtsg))
                    seen_rreps.append(rrep)

    return reacs


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
    assert not reac_keys & set(util.flatten(cons_keys)), (
        f"Assumption fails! Conserved stereo sites include reacting atoms:"
        f"\nConserved stereo sites: {cons_keys}"
        f"\nReacting atoms: {reac_keys}")

    floc_par_dct = stereo_parities(ftsg_loc)
    rloc_par_dct = stereo_parities(rtsg_loc)

    floc_cons_pars = list(map(floc_par_dct.__getitem__, cons_keys))
    rloc_cons_pars = list(map(rloc_par_dct.__getitem__, cons_keys))

    is_consistent = floc_cons_pars == rloc_cons_pars

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
