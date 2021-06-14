""" igraph interface

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools
import igraph
from phydat import ptab
from automol.util import dict_
from automol.graph.base._core import from_atoms_and_bonds
from automol.graph.base._core import atoms
from automol.graph.base._core import bonds
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bond_keys


def from_graph(gra):
    """ igraph object from a molecular graph
    """
    atm_keys = sorted(atom_keys(gra))
    bnd_keys = sorted(bond_keys(gra), key=sorted)

    atm_vals = dict_.values_by_key(atoms(gra), atm_keys)
    bnd_vals = dict_.values_by_key(bonds(gra), bnd_keys)

    atm_colors = list(itertools.starmap(_encode_vertex_attributes, atm_vals))
    bnd_colors = list(itertools.starmap(_encode_edge_attributes, bnd_vals))

    atm_idx_dct = dict(map(reversed, enumerate(atm_keys)))
    bnd_idxs = [sorted(map(atm_idx_dct.__getitem__, k)) for k in bnd_keys]
    igr = igraph.Graph(bnd_idxs)

    igr.vs['keys'] = atm_keys
    igr.vs['color'] = atm_colors
    igr.es['color'] = bnd_colors

    return igr


def to_graph(igr):
    """ igraph object from a molecular graph
    """
    atm_keys = igr.vs['keys']
    atm_key_dct = dict(enumerate(atm_keys))

    bnd_idxs = [e.tuple for e in igr.es]
    bnd_keys = [frozenset(map(atm_key_dct.__getitem__, k)) for k in bnd_idxs]

    atm_colors = igr.vs['color']
    bnd_colors = igr.es['color']

    atm_vals = list(map(_decode_vertex_attributes, atm_colors))
    bnd_vals = list(map(_decode_edge_attributes, bnd_colors))

    atms = dict(zip(atm_keys, atm_vals))
    bnds = dict(zip(bnd_keys, bnd_vals))
    gra = from_atoms_and_bonds(atms, bnds)
    return gra


def _encode_vertex_attributes(sym, imp_hyd_vlc, par):
    """ encode vertex attributes as an integer (or "color")

    scheme:
        atomic symbol               <=> hundreds place, as atomic number
        implicit hydrogen valence   <=> tens place
        parity                      <=> ones place (None->0, False->1, True->2)
    """
    id3 = ptab.to_number(sym)
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

    sym = ptab.to_symbol(id3)
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


def isomorphisms(igr1, igr2):
    """ get the list of automorphisms for an igraph object
    """
    atm_colors1 = igr1.vs['color']
    bnd_colors1 = igr1.es['color']
    atm_colors2 = igr2.vs['color']
    bnd_colors2 = igr2.es['color']
    isos = igr2.get_isomorphisms_vf2(
        other=igr1,
        color1=atm_colors2, color2=atm_colors1,
        edge_color1=bnd_colors2, edge_color2=bnd_colors1,
    )

    atm_keys1 = igr1.vs['keys']
    atm_keys2 = igr2.vs['keys']
    atm_keys2_dct = dict(enumerate(atm_keys2))
    iso_dcts = [dict(zip(atm_keys1, map(atm_keys2_dct.__getitem__, i)))
                for i in isos]
    return iso_dcts


def automorphisms(igr):
    """ get the list of automorphisms for an igraph object
    """
    return isomorphisms(igr, igr)


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
#
#
# if __name__ == '__main__':
#     GRA1 = ({0: ('C', 0, None), 1: ('H', 0, None), 2: ('H', 0, None),
#              3: ('H', 0, None), 4: ('H', 0, None)},
#             {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
#              frozenset({0, 4}): (1, None), frozenset({0, 3}): (1, None)})
#     GRA2 = automol.graph.relabel(GRA1, {0: 1, 1: 0})
#     GRA2 = automol.graph.relabel(GRA2, dict(enumerate(range(5, 10))))
#     print(automol.graph.string(GRA2, one_indexed=False))
#     IGR1 = from_graph(GRA1)
#     IGR2 = from_graph(GRA2)
#     isomorphisms(IGR1, IGR2)
#     # print(IGR1)
#     # print(IGR2)
