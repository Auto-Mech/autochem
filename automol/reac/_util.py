""" Common utilities for reaction classes
"""
import automol.graph
from automol.graph import ts
from automol import par


def hydrogen_migration_atom_keys(rxn):
    """ Obtain the atoms involved in a hydrogen migration reaction, sorted in
    canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the attacking atom, the transferring atom, the donating atom, and
    a neighbor to the attacking atom along the chain to the donating atom
    :rtype: (int, int, int, int)
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    tra_key, = frm_bnd_key & brk_bnd_key
    att_key, = frm_bnd_key - brk_bnd_key
    don_key, = brk_bnd_key - frm_bnd_key

    gra = ts.reactants_graph(rxn.forward_ts_graph)
    path = automol.graph.shortest_path_between_atoms(gra, att_key, don_key)
    ngb_key = automol.graph.atom_neighbor_atom_key(
        gra, att_key, incl_atm_keys=path)
    return att_key, tra_key, don_key, ngb_key


def ring_forming_scission_atom_keys(rxn):
    """ Obtain the atoms involved in a ring-forming scission reaction, sorted in
    canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the attacking atom, the transferring atom, the donating atom
    :rtype: (int, int, int, int)
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    tra_key, = frm_bnd_key & brk_bnd_key
    att_key, = frm_bnd_key - brk_bnd_key
    don_key, = brk_bnd_key - frm_bnd_key
    return att_key, tra_key, don_key


def ring_forming_scission_chain(rxn):
    """ Obtain the chain in a ring-forming scission reaction from the donating
    attom to the attacking atom.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: atoms along the chain
    :rtype: tuple[int]
    """
    att_key, _, don_key = ring_forming_scission_atom_keys(rxn)
    gra = ts.reactants_graph(rxn.forward_ts_graph)
    path = automol.graph.shortest_path_between_atoms(gra, don_key, att_key)
    return tuple(path)


def insertion_forming_bond_keys(rxn):
    """ Obtain the forming bonds for an insertion reaction, sorted in canonical
    order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the forming bond keys
    :rtype: (frozenset[int], frozenset[int])
    """
    assert rxn.class_ == par.ReactionClass.INSERTION
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    # Choose the forming bond that doesn't intersect with the breaking bond, if
    # one of them does
    frm_bnd_keys = sorted(ts.forming_bond_keys(rxn.forward_ts_graph),
                          key=sorted)
    frm_bnd_keys = sorted(frm_bnd_keys,
                          key=lambda x: len(x & brk_bnd_key))
    return tuple(frm_bnd_keys)
