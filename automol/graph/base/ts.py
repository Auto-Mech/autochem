""" transition state graph data structure

Forming bonds are encoded as 0.1-order bonds
Breaking bonds are encoded as 0.9-order bonds

Otherwise, this is equivalent to any other graph

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
from automol import util
from automol.util import dict_
from automol.graph.base._core import bond_orders
from automol.graph.base._core import stereo_parities
from automol.graph.base._core import set_bond_orders
from automol.graph.base._core import stereo_keys
from automol.graph.base._core import add_bonds
from automol.graph.base._core import remove_bonds
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import has_stereo
from automol.graph.base._algo import rings_bond_keys
from automol.graph.base._algo import sorted_ring_atom_keys_from_bond_keys
from automol.graph.base._canon import to_local_stereo as _to_local_stereo
from automol.graph.base._canon import from_local_stereo as _from_local_stereo
from automol.graph.base._canon import stereogenic_keys
from automol.graph.base._stereo import expand_stereo as _expand_stereo
from automol.graph.base._stereo import expand_stereo_with_priorities
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
    frm_bnd_keys = forming_bond_keys(tsg)
    ord_dct = dict_.transform_values(bond_orders(tsg), func=round)
    gra = set_bond_orders(tsg, ord_dct)
    gra = remove_bonds(gra, frm_bnd_keys)
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


def expand_stereo(tsg, sym_filter=True):
    """ Obtain all possible stereoisomers of a TS graph, ignoring its assignments

        :param tsg: A TS graph, without stereo
        :type tsg: automol graph data structure
        :param sym_filter: filter out symmetrically equivalent stereoisomers?
        :type sym_filter: bool
        :returns: a series of molecular graphs for the stereoisomers
    """
    return _expand_stereo(tsg, sym_filter=sym_filter)


def expand_reaction_stereo(tsg, sym_filter=True):
    """ Obtain all possible stereoisomer combinations for a reaction, encoding
        reactant and product stereo assignments in forward and reverse TS
        graphs, respectively

        :param tsg: A TS graph; If it has stereo assignments, they will be
            assumed to be for the reactants, and all compatible products
            assignments will be exapnded
        :type tsg: automol graph data structure
        :param sym_filter: filter out symmetrically equivalent stereoisomers?
        :type sym_filter: bool
        :returns: a series of pairs of forward and reverse graphs containing
            mutually compatible stereo assignments
    """
    if has_stereo(tsg):
        assert not stereogenic_keys(tsg)
        ftsgs = [tsg]
    else:
        ftsgs = expand_stereo(tsg, sym_filter=True)

    rtsgs, rpri_dcts = zip(*expand_stereo_with_priorities(reverse(tsg)))

    pairs = []

    for ftsg in ftsgs:
        ftsg_loc = to_local_stereo(ftsg)

        seen_rreps = []
        for rtsg, rpri_dct in zip(rtsgs, rpri_dcts):
            rrep = canonical_assignment_representation(rtsg, rpri_dct)
            rtsg_loc = to_local_stereo(rtsg)

            if reaction_has_consistent_stereo(ftsg_loc, rtsg_loc):
                if not sym_filter or rrep not in seen_rreps:
                    pairs.append((ftsg, rtsg))
                    seen_rreps.append(rrep)

    return pairs


def reaction_has_consistent_stereo(ftsg_loc, rtsg_loc):
    """ Does this reaction have consistent stereo?

        :param ftsg_loc: a forward TS graph, with local stereo assignments for
            the reactants
        :type ftsg_loc: automol graph data structure
        :param ftsg_loc: a forward TS graph, with local stereo assignments for
            the products
        :type rtsg_loc: automol graph data structure
        :rtype: bool
    """
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

    return floc_cons_pars == rloc_cons_pars
