""" transition state graph data structure

Forming bonds are encoded as 0.1-order bonds
Breaking bonds are encoded as 0.9-order bonds

Otherwise, this is equivalent to any other graph

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
from automol.util import dict_
from automol.graph.base._core import bond_orders
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import set_bond_orders
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import add_bonds
from automol.graph.base._core import remove_bonds
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._algo import rings_bond_keys
from automol.graph.base._algo import sorted_ring_atom_keys_from_bond_keys
from automol.graph.base._stereo import stereomers as _stereomers
from automol.graph.base._stereo import stereogenic_atom_keys
from automol.graph.base._stereo import stereogenic_bond_keys
from automol.graph.base._stereo import (to_index_based_stereo as
                                        _to_index_based_stereo)
from automol.graph.base._stereo import (from_index_based_stereo as
                                        _from_index_based_stereo)


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


def nonconserved_atom_stereo_keys(ste_tsg, check=True):
    """ Determine atom stereo centers which are not conserved by the reaction.

    This includes atoms which are stereogenic for the products but not for the
    reactants ("created" stereo centers) and atoms which are stereogenic for
    the reactants but not for the products ("destroyed" stereo centers).

    :param ste_tsg: The TS graph, with stereo assignments.
    :returns: Created and destroyed atom stereo centers, respectively.
    :rtype: (frozenset, frozenset)
    """
    rcts_gra = reactants_graph(ste_tsg)
    prds_gra = products_graph(ste_tsg)

    if check:
        ste_atm_keys = stereogenic_atom_keys(rcts_gra)
        ste_bnd_keys = stereogenic_bond_keys(rcts_gra)
        assert not ste_atm_keys, (
            "Unassigned atom stereo centers: {}".format(str(ste_atm_keys)))
        assert not ste_bnd_keys, (
            "Unassigned bond stereo centers: {}".format(str(ste_bnd_keys)))

    keys1 = atom_stereo_keys(rcts_gra)
    keys2 = stereogenic_atom_keys(prds_gra, assigned=True)

    cre_ste_atm_keys = frozenset(keys2 - keys1)
    des_ste_atm_keys = frozenset(keys1 - keys2)

    return cre_ste_atm_keys, des_ste_atm_keys


def nonconserved_bond_stereo_keys(ste_tsg, check=True):
    """ Determine bond stereo centers which are not conserved by the reaction.

    This includes bonds which are stereogenic for the products but not for the
    reactants ("created" stereo centers) and bonds which are stereogenic for
    the reactants but not for the products ("destroyed" stereo centers).

    :param ste_tsg: The TS graph, with stereo assignments.
    :returns: Created and destroyed bond stereo centers, respectively.
    :rtype: (frozenset, frozenset)
    """
    rcts_gra = reactants_graph(ste_tsg)
    prds_gra = products_graph(ste_tsg)

    if check:
        ste_atm_keys = stereogenic_atom_keys(rcts_gra)
        ste_bnd_keys = stereogenic_bond_keys(rcts_gra)
        assert not ste_atm_keys, (
            "Unassigned atom stereo centers: {}".format(str(ste_atm_keys)))
        assert not ste_bnd_keys, (
            "Unassigned bond stereo centers: {}".format(str(ste_bnd_keys)))

    keys1 = bond_stereo_keys(rcts_gra)
    keys2 = stereogenic_bond_keys(prds_gra, assigned=True)

    cre_ste_bnd_keys = frozenset(keys2 - keys1)
    des_ste_bnd_keys = frozenset(keys1 - keys2)

    return cre_ste_bnd_keys, des_ste_bnd_keys


def to_index_based_stereo(ste_tsg):
    """ Convert a TS graph to index-based stereo assignments, where parities
    are defined relative to the ordering of indices rather than the absolute
    stereo priority.

    :param ste_tsg: a TS graph with absolute stereo assignments
    :returns: a TS graph with index-based stereo assignments
    """
    rcts_gra = reactants_graph(ste_tsg)
    frm_bnd_keys = forming_bond_keys(ste_tsg)
    brk_bnd_keys = breaking_bond_keys(ste_tsg)

    rcts_gra = _to_index_based_stereo(rcts_gra)
    idx_tsg = graph(rcts_gra, frm_bnd_keys, brk_bnd_keys)
    return idx_tsg


def from_index_based_stereo(idx_tsg):
    """ Convert a TS graph from index-based stereo assignments back to absolute
    stereo assignments, where parities are independent of atom ordering.

    :param idx_tsg: a TS graph with index-based stereo assignments
    :returns: a TS graph with absolute stereo assignments
    """
    rcts_gra = reactants_graph(idx_tsg)
    frm_bnd_keys = forming_bond_keys(idx_tsg)
    brk_bnd_keys = breaking_bond_keys(idx_tsg)

    rcts_gra = _from_index_based_stereo(rcts_gra)
    ste_tsg = graph(rcts_gra, frm_bnd_keys, brk_bnd_keys)
    return ste_tsg


def stereomers(tsg):
    """ Expand all possible stereo assignments for the reactants in this TS
    graph. (Ignores stereo assignments already present, if any.)

    :param tsg: The TS graph, without stereo assignments.
    :returns: All possible TS graphs with stereo assignments for the reactants.
    """
    rcts_gra = reactants_graph(tsg)
    frm_bnd_keys = forming_bond_keys(tsg)
    brk_bnd_keys = breaking_bond_keys(tsg)

    rcts_gra = without_stereo_parities(rcts_gra)
    rcts_sgrs = _stereomers(rcts_gra)
    ste_tsgs = tuple(
        graph(rcts_sgr, frm_bnd_keys, brk_bnd_keys) for rcts_sgr in rcts_sgrs)
    return ste_tsgs


def compatible_reverse_stereomers(ste_tsg):
    """ Given a TS graph with stereo assignments, expand all possible reverse
    graphs compatble with the forward graph.

    :param ste_tsg: The TS graph, with stereo assignments.
    :returns: All possible reverse TS graphs.
    """
    frm_bnd_keys = forming_bond_keys(ste_tsg)
    brk_bnd_keys = breaking_bond_keys(ste_tsg)
    _, des_ste_atm_keys = nonconserved_atom_stereo_keys(ste_tsg)
    _, des_ste_bnd_keys = nonconserved_bond_stereo_keys(ste_tsg)
    cons_atm_keys = sorted(atom_stereo_keys(ste_tsg) - des_ste_atm_keys)
    cons_bnd_keys = sorted(bond_stereo_keys(ste_tsg) - des_ste_bnd_keys)

    # 1. Determine index-based stereo assignments for conserved stereo centers
    idx_tsg = to_index_based_stereo(ste_tsg)
    cons_idx_atm_pars = dict_.values_by_key(
        atom_stereo_parities(idx_tsg), cons_atm_keys)
    cons_idx_bnd_pars = dict_.values_by_key(
        bond_stereo_parities(idx_tsg), cons_bnd_keys)

    # 2. Determine all possible index-based stereo assignments for the reverse
    #    reaction.
    prds_gra = without_stereo_parities(products_graph(ste_tsg))
    prds_sgrs = _stereomers(prds_gra)
    prds_idx_sgrs = list(map(_to_index_based_stereo, prds_sgrs))
    rev_idx_tsgs_pool = [
        graph(p, brk_bnd_keys, frm_bnd_keys) for p in prds_idx_sgrs]

    # 3. Find possibilities which match the assignments for the conserved
    #    stereo centers.
    rev_idx_tsgs = []
    for rev_idx_tsg in rev_idx_tsgs_pool:
        rev_cons_idx_atm_pars = dict_.values_by_key(
            atom_stereo_parities(rev_idx_tsg), cons_atm_keys)
        rev_cons_idx_bnd_pars = dict_.values_by_key(
            bond_stereo_parities(rev_idx_tsg), cons_bnd_keys)
        if (rev_cons_idx_atm_pars == cons_idx_atm_pars and
                rev_cons_idx_bnd_pars == cons_idx_bnd_pars):
            rev_idx_tsgs.append(rev_idx_tsg)

    # 4. Convert the matching reverse graphs back from index-based stereo
    #    assignments to absolute stereo assignments.
    rev_ste_tsgs = list(map(from_index_based_stereo, rev_idx_tsgs))
    return rev_ste_tsgs
