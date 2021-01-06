""" igraph interface
"""

import itertools
from qcelemental import periodictable as pt
import igraph
import automol
import automol.create.graph
from automol.util import dict_


def from_graph(gra):
    """ igraph object from a molecular graph
    """
    atm_labels = sorted(automol.graph.atom_keys(gra))
    bnd_labels = list(sorted(map(sorted, automol.graph.bond_keys(gra))))

    atm_keys = atm_labels
    bnd_keys = list(map(frozenset, bnd_labels))

    atm_vals = dict_.values_by_key(automol.graph.atoms(gra), atm_keys)
    bnd_vals = dict_.values_by_key(automol.graph.bonds(gra), bnd_keys)

    atm_colors = list(itertools.starmap(_encode_vertex_attributes, atm_vals))
    bnd_colors = list(itertools.starmap(_encode_edge_attributes, bnd_vals))

    igr = igraph.Graph(bnd_labels)

    igr.vs['label'] = atm_labels
    igr.es['label'] = bnd_labels

    igr.vs['color'] = atm_colors
    igr.es['color'] = bnd_colors

    return igr


def to_graph(igr):
    """ igraph object from a molecular graph
    """
    atm_labels = igr.vs['label']
    bnd_labels = igr.es['label']

    atm_colors = igr.vs['color']
    bnd_colors = igr.es['color']

    atm_keys = atm_labels
    bnd_keys = list(map(frozenset, bnd_labels))

    atm_vals = list(map(_decode_vertex_attributes, atm_colors))
    bnd_vals = list(map(_decode_edge_attributes, bnd_colors))

    atms = dict(zip(atm_keys, atm_vals))
    bnds = dict(zip(bnd_keys, bnd_vals))
    gra = automol.create.graph.from_atoms_and_bonds(atms, bnds)
    return gra


def _encode_vertex_attributes(sym, imp_hyd_vlc, par):
    """ encode vertex attributes as an integer (or "color")

    scheme:
        atomic symbol               <=> hundreds place, as atomic number
        implicit hydrogen valence   <=> tens place
        parity                      <=> ones place (None->0, False->1, True->2)
    """
    id3 = pt.to_Z(sym)
    id2 = imp_hyd_vlc
    id1 = 0 if par is None else 1 + int(par)

    color = id3 * 100 + id2 * 10 + id1 * 1
    return color


def _encode_edge_attributes(order, par):
    """ encode bond attributes as an integer (or "color")

    scheme:
        bond order   <=> tens place
        parity       <=> ones place (None->0, False->1, True->2)
    """
    id2 = order
    id1 = 0 if par is None else 1 + int(par)

    color = id2 * 10 + id1 * 1
    return color


def _decode_vertex_attributes(color):
    """ decode vertex attributes from integer values

    scheme:
        atomic symbol               <=> hundreds place, as atomic number
        implicit hydrogen valence   <=> tens place
        parity                      <=> ones place (None->0, False->1, True->2)
    """
    id3 = color // 100
    color -= id3 * 100
    id2 = color // 10
    color -= id2 * 10
    id1 = color // 1

    sym = pt.to_E(id3)
    imp_hyd_vlc = id2
    assert id1 in (0, 1, 2)
    par = None if id1 == 0 else bool(id1 - 1)

    return sym, imp_hyd_vlc, par


def _decode_edge_attributes(color):
    """ decode vertex attributes from integer values

    scheme:
        bond order   <=> tens place
        parity       <=> ones place (None->0, False->1, True->2)
    """
    id2 = color // 10
    color -= id2 * 10
    id1 = color // 1

    order = id2
    assert id1 in (0, 1, 2)
    par = None if id1 == 0 else bool(id1 - 1)

    return order, par


def automorphisms(igr):
    """ get the list of automorphisms for an igraph object
    """
    atm_keys = igr.vs['label']
    atm_colors = igr.vs['color']
    bnd_colors = igr.es['color']
    auts = igr.get_automorphisms_vf2(color=atm_colors, edge_color=bnd_colors)
    aut_dcts = [dict(zip(atm_keys, aut)) for aut in auts]
    return aut_dcts


def canonical_permutation(igr):
    """ get the list of automorphisms for an igraph object

    (currently does not consider bond stereo, because the igraph/BLISS function
    doesn't consider edge colors)
    """
    atm_keys = igr.vs['label']
    atm_colors = igr.vs['color']
    perm = igr.canonical_permutation(color=atm_colors)
    perm_dct = dict(zip(atm_keys, perm))
    return perm_dct
