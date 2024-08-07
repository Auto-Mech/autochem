""" Common utilities for reaction classes
"""

import itertools

from .. import graph
from ..const import ReactionClass
from ..graph import ts
from ._0core import Reaction, class_, ts_graph


def hydrogen_migration_atom_keys(rxn: Reaction):
    """Obtain the atoms involved in a hydrogen migration reaction, sorted in
    canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the attacking atom, the transferring atom, the donating atom, and
    a neighbor to the attacking atom along the chain to the donating atom
    :rtype: (int, int, int, int)
    """
    (frm_bnd_key,) = ts.forming_bond_keys(ts_graph(rxn))
    (brk_bnd_key,) = ts.breaking_bond_keys(ts_graph(rxn))
    (tra_key,) = frm_bnd_key & brk_bnd_key
    (att_key,) = frm_bnd_key - brk_bnd_key
    (don_key,) = brk_bnd_key - frm_bnd_key

    gra = ts.reactants_graph(ts_graph(rxn))
    path = graph.shortest_path_between_atoms(gra, att_key, don_key)
    ngb_key = graph.atom_neighbor_atom_key(gra, att_key, incl_keys=path)
    return att_key, tra_key, don_key, ngb_key


def hydrogen_migration_might_dissociate(rxn: Reaction, att_key, ngb_key, don_key):
    """Obtain the atoms involved when a migration might be similar to an HO2
        elimination

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the dissociation oxygen atom, the dissociation carbon atom
    :rtype: (int, int, int, int)
    """
    gra = ts.reactants_graph(ts_graph(rxn))
    atm_symbs = graph.atom_symbols(gra)
    o_ngb = []
    diss_key = None
    if atm_symbs[att_key] == "O" and atm_symbs[ngb_key] == "O":
        o_ngbs = graph.atom_neighbor_atom_keys(gra, ngb_key)
        o_ngb = [ngb for ngb in o_ngbs if ngb != att_key]
    if len(o_ngb) == 1:
        if don_key in graph.atom_neighbor_atom_keys(gra, o_ngb[0]):
            diss_key = (
                ngb_key,
                o_ngb[0],
            )
    return diss_key


def ring_forming_scission_atom_keys(rxn: Reaction):
    """Obtain the atoms involved in a ring-forming scission reaction, sorted
        in canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the attacking atom, the transferring atom, the donating atom
    :rtype: (int, int, int, int)
    """
    (frm_bnd_key,) = ts.forming_bond_keys(ts_graph(rxn))
    (brk_bnd_key,) = ts.breaking_bond_keys(ts_graph(rxn))
    (tra_key,) = frm_bnd_key & brk_bnd_key
    (att_key,) = frm_bnd_key - brk_bnd_key
    (don_key,) = brk_bnd_key - frm_bnd_key
    return att_key, tra_key, don_key


def ring_forming_scission_chain(rxn: Reaction):
    """Obtain the chain in a ring-forming scission reaction from the donating
    attom to the attacking atom.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: atoms along the chain
    :rtype: tuple[int]
    """
    att_key, _, don_key = ring_forming_scission_atom_keys(rxn)
    gra = ts.reactants_graph(ts_graph(rxn))
    path = graph.shortest_path_between_atoms(gra, don_key, att_key)
    return tuple(path)


def hydrogen_abstraction_atom_keys(rxn: Reaction):
    """Obtain the atoms involved in a hydrogen abstraction, sorted in
    canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the attacking atom, the transferring atom, the donating atom
    :rtype: (int, int, int)
    """
    (frm_bnd_key,) = ts.forming_bond_keys(ts_graph(rxn))
    (brk_bnd_key,) = ts.breaking_bond_keys(ts_graph(rxn))
    (hyd_key,) = frm_bnd_key & brk_bnd_key
    (att_key,) = frm_bnd_key - brk_bnd_key
    (don_key,) = brk_bnd_key - frm_bnd_key
    return att_key, hyd_key, don_key


def elimination_breaking_bond_keys(rxn: Reaction):
    """Obtain the breaking bonds for an elimination reaction, sorted in
        canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the breaking bond keys
    :rtype: (frozenset[int], frozenset[int])
    """
    assert class_(rxn) == ReactionClass.ELIMINATION
    tsg = ts_graph(rxn)
    (frm_bnd_key,) = ts.forming_bond_keys(tsg)
    brk_bnd_keys = ts.breaking_bond_keys(tsg)
    brk_bnd_key1, brk_bnd_key2 = brk_bnd_keys

    symbs = graph.atom_symbols(tsg)

    if len(frm_bnd_key | brk_bnd_key1 | brk_bnd_key2) > 3:
        # for ring_size > 3: use brk-bnd that doesn't involve atoms in frm bond
        scn_brk_bnd_key = None
        for brk_bnd_key in brk_bnd_keys:
            if not frm_bnd_key & brk_bnd_key:
                scn_brk_bnd_key = brk_bnd_key
    else:
        # if one brk bnd doesn't have H, use that, else use arbitrary brk bnd
        scn_brk_bnd_key = None
        for brk_bnd_key in brk_bnd_keys:
            brk_symbs = tuple(symbs[key] for key in brk_bnd_key1)
            if "H" not in brk_symbs:
                scn_brk_bnd_key = brk_bnd_key
                break
        if scn_brk_bnd_key is None:
            scn_brk_bnd_key = brk_bnd_key1

    brk_bnd_key1 = scn_brk_bnd_key
    (brk_bnd_key2,) = brk_bnd_keys - {brk_bnd_key1}

    return brk_bnd_key1, brk_bnd_key2


def insertion_forming_bond_keys(rxn: Reaction):
    """Obtain the forming bonds for an insertion reaction, sorted in canonical
    order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the forming bond keys
    :rtype: (frozenset[int], frozenset[int])
    """

    assert class_(rxn) == ReactionClass.INSERTION
    # Choose the forming bond with the fewest neighbors, to get the terminal
    # atom if there is one.
    tsg = ts_graph(rxn)
    # For tricky reasons, we need to sort in descending order here.
    # The reason is that the ordering of bonds here determines the ordering of
    # atoms in the z-matrix, which means this function could otherwise yield
    # inconsistent results (i.e. the two bonds will reverse order).
    # If needed, we could add sorting that is based on the symbols of the
    # atoms instead.
    frm_bnd_keys = reversed(sorted(ts.forming_bond_keys(tsg), key=sorted))
    frm_bnd_keys = sorted(
        frm_bnd_keys, key=lambda x: graph.atom_count(graph.bond_neighborhood(tsg, x))
    )

    return tuple(frm_bnd_keys)


def substitution_atom_keys(rxn: Reaction):
    """Obtain the atoms involved in a substitution reaction, sorted in
        canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the attacking atom, the transferring atom, the leaving atom
    :rtype: (int, int, int)
    """
    (frm_bnd_key,) = ts.forming_bond_keys(ts_graph(rxn))
    (brk_bnd_key,) = ts.breaking_bond_keys(ts_graph(rxn))
    (tra_key,) = frm_bnd_key & brk_bnd_key
    (att_key,) = frm_bnd_key - brk_bnd_key
    (lea_key,) = brk_bnd_key - frm_bnd_key
    return att_key, tra_key, lea_key


def assert_is_valid_reagent_graph_list(gras):
    """Assert that a sequence of graphs has the appropriate form for reactants
    or products in a reaction object

    The sequence is appropriate if every graph is explicit and none of them
    have overlapping atom keys.

    :param gras: the graphs
    :type gras: list

    """
    gras_str = "\n---\n".join(map(graph.string, gras))
    assert _are_all_explicit(
        gras
    ), f"Implicit hydrogens are not allowed here!\nGraphs:\n{gras_str}"
    assert _have_no_common_atom_keys(
        gras
    ), f"Overlapping atom keys are not allowed here!\nGraphs:\n{gras_str}"
    assert _have_no_partial_stereo(
        gras
    ), f"If present, stereo must be complete!\nGraphs:\n{gras_str}"


def _are_all_explicit(gras):
    return all(gra == graph.explicit(gra) for gra in gras)


def _have_no_common_atom_keys(gras):
    atm_keys = list(itertools.chain(*map(graph.atom_keys, gras)))
    return len(atm_keys) == len(set(atm_keys))


def _have_no_partial_stereo(gras):
    if any(map(graph.has_stereo, gras)):
        ret = not any(map(graph.unassigned_stereocenter_keys, gras))
    else:
        ret = True
    return ret
