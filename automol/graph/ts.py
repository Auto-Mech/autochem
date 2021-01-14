""" transition state graph data structure

Forming bonds are encoded as 0.1-order bonds
Breaking bonds are encoded as 0.9-order bonds

Otherwise, this is equivalent to any other graph
"""
from automol.util import dict_
from automol.graph._graph_base import bond_orders
from automol.graph._graph_base import set_bond_orders
from automol.graph._graph import add_bonds
from automol.graph._graph import remove_bonds
from automol.graph._ring import rings_bond_keys
from automol.graph._ring import sorted_ring_atom_keys_from_bond_keys


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


def reverse(tsg):
    """ reverse a transition state graph
    """
    return graph(gra=tsg,
                 frm_bnd_keys=breaking_bond_keys(tsg),
                 brk_bnd_keys=forming_bond_keys(tsg))


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


def vmatrix(tsg):
    """ v-matrix for a TS graph
    """
    import automol
    print(automol.graph.string(tsg, one_indexed=False))
