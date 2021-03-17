""" molecular graph
"""

import operator
import itertools
import functools
import numpy
import future.moves.itertools as fmit
from phydat import ptab
import automol.formula
from automol.util import dict_
from automol.graph._graph_base import atoms
from automol.graph._graph_base import bonds
from automol.graph._graph_base import atom_keys
from automol.graph._graph_base import bond_keys
from automol.graph._graph_base import atom_symbols
from automol.graph._graph_base import atom_implicit_hydrogen_valences
from automol.graph._graph_base import atom_stereo_parities
from automol.graph._graph_base import bond_orders
from automol.graph._graph_base import bond_stereo_parities
from automol.graph._graph_base import set_atom_implicit_hydrogen_valences
from automol.graph._graph_base import set_atom_stereo_parities
from automol.graph._graph_base import set_bond_orders
from automol.graph._graph_base import set_bond_stereo_parities
from automol.graph._graph_base import relabel
from automol.graph import _networkx
from automol.graph import _igraph
import automol.create.graph as _create


# setters
def standard_keys(gra):
    """ replace the current atom keys with standard indices, counting from zero
    """
    atm_key_dct = dict(map(reversed, enumerate(sorted(atom_keys(gra)))))
    return relabel(gra, atm_key_dct)


def standard_keys_for_sequence(gras):
    """ assigns non-overlapping keys to a sequence of graphs

    (returns a series of key maps for each)
    """
    atm_key_dcts = []

    shift = 0
    for gra in gras:
        natms = atom_count(gra, dummy=True, with_implicit=False)

        atm_key_dct = {atm_key: idx+shift
                       for idx, atm_key in enumerate(sorted(atom_keys(gra)))}
        atm_key_dcts.append(atm_key_dct)

        shift += natms

    gras = [relabel(gra, atm_key_dct)
            for gra, atm_key_dct in zip(gras, atm_key_dcts)]

    return gras, atm_key_dcts


def transform_keys(gra, atm_key_func):
    """ transform atom keys with a function
    """
    atm_keys = atom_keys(gra)
    atm_key_dct = dict(zip(atm_keys, map(atm_key_func, atm_keys)))
    return relabel(gra, atm_key_dct)


def relabel_for_zmatrix(gra, zma_keys, dummy_key_dct):
    """ relabel a geometry graph to line up with a z-matrix

    Inserts dummy atoms and sorts/relabels in z-matrix order.

    Graph keys should correspond to the geometry used for conversion.

    :param gra: the graph
    :param zma_keys: graph keys in the order they appear in the z-matrix
    :param dummy_key_dct: dummy keys introduced on z-matrix conversion, by atom
        they are attached to
    """
    gra = add_dummy_atoms(gra, dummy_key_dct)
    key_dct = dict(map(reversed, enumerate(zma_keys)))
    gra = relabel(gra, key_dct)
    return gra


def relabel_for_geometry(gra):
    """ relabel a z-matrix graph to line up with a geometry

    The result will line up with a geometry converted from the z-matrix, with
    dummy atoms removed.

    Removes dummy atoms and relabels in geometry order.

    Graph keys should correspond to the z-matrix used for conversion.

    :param gra: the graph
    """
    dummy_keys = sorted(atom_keys(gra, sym='X'))

    for dummy_key in reversed(dummy_keys):
        gra = _shift_remove_dummy_atom(gra, dummy_key)
    return gra


def _shift_remove_dummy_atom(gra, dummy_key):
    keys = sorted(automol.graph.atom_keys(gra))
    idx = keys.index(dummy_key)
    key_dct = {}
    key_dct.update({k: k for k in keys[:idx]})
    key_dct.update({k: k-1 for k in keys[(idx+1):]})
    gra = remove_atoms(gra, [dummy_key])
    gra = relabel(gra, key_dct)
    return gra


def without_bond_orders(gra):
    """ resonance graph with maximum spin (i.e. no pi bonds)
    """
    bnd_keys = list(bond_keys(gra))
    # don't set dummy bonds to one!
    bnd_ord_dct = bond_orders(gra)
    bnd_vals = [1 if v != 0 else 0
                for v in map(bnd_ord_dct.__getitem__, bnd_keys)]
    bnd_ord_dct = dict(zip(bnd_keys, bnd_vals))
    return set_bond_orders(gra, bnd_ord_dct)


def without_stereo_parities(gra):
    """ graph with stereo assignments wiped out
    """
    atm_ste_par_dct = dict_.by_key({}, atom_keys(gra), fill_val=None)
    bnd_ste_par_dct = dict_.by_key({}, bond_keys(gra), fill_val=None)
    gra = set_atom_stereo_parities(gra, atm_ste_par_dct)
    gra = set_bond_stereo_parities(gra, bnd_ste_par_dct)
    return gra


def frozen(gra):
    """ hashable, sortable, immutable container of graph data
    """
    atm_keys = sorted(atom_keys(gra))
    bnd_keys = sorted(bond_keys(gra), key=sorted)

    # make it sortable by replacing Nones with -infinity
    atm_vals = numpy.array(dict_.values_by_key(atoms(gra), atm_keys))
    bnd_vals = numpy.array(dict_.values_by_key(bonds(gra), bnd_keys))
    atm_vals[numpy.equal(atm_vals, None)] = -numpy.inf
    bnd_vals[numpy.equal(bnd_vals, None)] = -numpy.inf

    frz_atms = tuple(zip(atm_keys, map(tuple, atm_vals)))
    frz_bnds = tuple(zip(bnd_keys, map(tuple, bnd_vals)))
    return (frz_atms, frz_bnds)


# graph theory library
# # atom properties
def electron_count(gra, charge=0):
    """ the number of electrons in the molecule
    """
    atm_symb_dct = atom_symbols(explicit(gra))
    nelec = sum(map(ptab.to_number, atm_symb_dct.values())) - charge
    return nelec


def atom_count(gra, dummy=False, with_implicit=True):
    """ count the number of atoms in this molecule

    by default, this includes implicit hydrogens and excludes dummy atoms
    """
    if not dummy:
        gra = without_dummy_atoms(gra)
    natms = len(atoms(gra))
    if with_implicit:
        atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(gra)
        natms += sum(atm_imp_hyd_vlc_dct.values())
    return natms


def atom_count_by_type(gra, sym, keys=None):
    """ count the number of atoms with a given type (symbol)

    :param gra: the graph
    :param sym: the symbol
    :param keys: optionally, restrict the count to a subset of keys
    """
    keys = atom_keys(gra) if keys is None else keys
    symb_dct = atom_symbols(gra)
    symbs = list(map(symb_dct.__getitem__, keys))
    return symbs.count(sym)


def heavy_atom_count(gra, dummy=False):
    """ the number of heavy atoms
    """
    if not dummy:
        gra = without_dummy_atoms(gra)
    atm_symb_dct = atom_symbols(gra)
    nhvy_atms = sum(ptab.to_number(sym) != 1 for sym in atm_symb_dct.values())
    return nhvy_atms


def atom_neighborhood(gra, atm_key, bnd_keys=None, stereo=False):
    """ neighborhood subgraph for a specific atom
    """
    bnd_keys = bond_keys(gra) if bnd_keys is None else bnd_keys
    nbh_bnd_keys = set(k for k in bnd_keys if atm_key in k)
    nbh = bond_induced_subgraph(gra, nbh_bnd_keys, stereo=stereo)
    return nbh


def atom_neighborhoods(gra, stereo=False):
    """ neighborhood subgraphs, by atom
    """
    bnd_keys = bond_keys(gra)

    def _neighborhood(atm_key):
        return atom_neighborhood(gra, atm_key, bnd_keys=bnd_keys,
                                 stereo=stereo)

    atm_keys = list(atom_keys(gra))
    atm_nbh_dct = dict(zip(atm_keys, map(_neighborhood, atm_keys)))
    return atm_nbh_dct


def atom_sorted_neighbor_atom_keys(gra, atm_key, excl_atm_keys=(),
                                   incl_atm_keys=None, symbs_first=('C',),
                                   symbs_last=('H',)):
    """ get the next in a sorted list of neighbor keys, excluding some
    """
    atm_symb_dct = atom_symbols(gra)
    incl_atm_keys = atom_keys(gra) if incl_atm_keys is None else incl_atm_keys

    atm_nbh = atom_neighborhood(gra, atm_key)
    atm_keys = sorted(atom_keys(atm_nbh) - {atm_key} - set(excl_atm_keys))
    atm_keys = [k for k in atm_keys if k in incl_atm_keys]

    symbs = list(map(atm_symb_dct.__getitem__, atm_keys))
    srt = automol.formula.argsort_symbols(symbs, symbs_first, symbs_last)
    atm_keys = tuple(map(atm_keys.__getitem__, srt))
    return atm_keys


def atom_neighbor_atom_key(gra, atm_key, excl_atm_keys=(), incl_atm_keys=None,
                           symbs_first=('C',), symbs_last=('H',)):
    """ get the next in a sorted list of neighbor keys, excluding some
    """
    atm_keys = atom_sorted_neighbor_atom_keys(
        gra, atm_key, excl_atm_keys=excl_atm_keys, incl_atm_keys=incl_atm_keys,
        symbs_first=symbs_first, symbs_last=symbs_last)
    return atm_keys[0] if atm_keys else None


def atoms_neighbor_atom_keys(gra):
    """ keys of neighboring atoms, by atom
    """
    def _neighbor_keys(atm_key, atm_nbh):
        return frozenset(atom_keys(atm_nbh) - {atm_key})

    atm_ngb_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(gra), _neighbor_keys)
    return atm_ngb_keys_dct


def atoms_sorted_neighbor_atom_keys(gra, symbs_first=('C',), symbs_last=('H',),
                                    ords_last=(0.1,), prioritize_keys=()):
    """ keys of neighboring atoms, by atom

    :param gra: the graph
    :param symbs_first: atomic symbols to put put first in the sort order
    :param symbs_last: atomic symbols to put last in the sort order
    :param ords_last: neighors connected with a bond of this order will be put
        last in the sort order
    :param prioritize_keys: keys to put first no matter what
    """
    atm_symb_dct = atom_symbols(gra)
    bnd_ord_dct = bond_orders(gra)

    def _neighbor_keys(atm_key, atm_nbh):
        keys = sorted(atom_keys(atm_nbh) - {atm_key})
        bnd_keys = [frozenset({atm_key, k}) for k in keys]
        ords = list(map(bnd_ord_dct.__getitem__, bnd_keys))
        ords = [-1 if o not in ords_last else ords_last.index(o)
                for o in ords]
        symbs = list(map(atm_symb_dct.__getitem__, keys))
        pris = [0 if k in prioritize_keys else 1 for k in keys]
        srt_vals = list(zip(ords, pris, symbs))
        srt = automol.formula.argsort_symbols(
            srt_vals, symbs_first, symbs_last, idx=1)
        keys = tuple(map(keys.__getitem__, srt))
        return keys

    atm_ngb_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(gra), _neighbor_keys)
    return atm_ngb_keys_dct


def atoms_second_degree_neighbor_atom_keys(gra):
    """ keys of second-degree neighboring atoms, by atom

    That is, atoms that are connected through an intermediate atom
    """
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    atm_ngb2_keys_dct = {}
    for atm_key, atm_ngb_keys in atm_ngb_keys_dct.items():
        # Union of neighbors of neighbors
        atm_ngb2_keys = functools.reduce(
            operator.or_, map(atm_ngb_keys_dct.__getitem__, atm_ngb_keys))
        # Subtract of the atom itself and its neighbors (in case of 3-rings)
        atm_ngb2_keys -= {atm_key} | atm_ngb_keys

        atm_ngb2_keys_dct[atm_key] = frozenset(atm_ngb2_keys)
    return atm_ngb2_keys_dct


def dummy_atoms_neighbor_atom_key(gra):
    """ Atoms that are connected to dummy atoms, by dummy atom key

    (Requires that each dummy atom only be connected to one neighbor)
    """
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    dummy_atm_keys = atom_keys(gra, sym='X')

    dummy_ngb_key_dct = {}
    for key in dummy_atm_keys:
        ngb_keys = atm_ngb_keys_dct[key]
        assert len(ngb_keys) == 1, (
            "Dummy atoms should only be connected to one atom!")
        ngb_key, = ngb_keys
        dummy_ngb_key_dct[key] = ngb_key

    return dummy_ngb_key_dct


def atoms_bond_keys(gra):
    """ bond keys, by atom
    """
    return dict_.transform_values(atom_neighborhoods(gra), bond_keys)


def angle_keys(gra):
    """ triples of keys for pairs of adjacent bonds, with the central atom in
    the middle
    """
    bnd_keys = bond_keys(gra)

    ang_keys = []
    for bnd_key1, bnd_key2 in itertools.combinations(bnd_keys, r=2):
        if bnd_key1 != bnd_key2 and bnd_key1 & bnd_key2:
            atm2_key, = bnd_key1 & bnd_key2
            atm1_key, = bnd_key1 - {atm2_key}
            atm3_key, = bnd_key2 - {atm2_key}
            ang_keys.append((atm1_key, atm2_key, atm3_key))
            ang_keys.append((atm3_key, atm2_key, atm1_key))

    return tuple(ang_keys)


# # bond properties
def bonds_neighbor_atom_keys(gra):
    """ keys of neighboring atoms, by bond
    """
    def _neighbor_keys(bnd_key, bnd_nbh):
        return frozenset(atom_keys(bnd_nbh) - bnd_key)

    bnd_ngb_keys_dct = dict_.transform_items_to_values(
        bond_neighborhoods(gra), _neighbor_keys)
    return bnd_ngb_keys_dct


def bonds_neighbor_bond_keys(gra):
    """ keys of neighboring bonds, by bond
    """
    def _neighbor_keys(bnd_key, bnd_nbh):
        bnd_keys = bond_keys(bnd_nbh)
        bnd_keys -= {bnd_key}
        bnd_keys = frozenset(key for key in bnd_keys if key & bnd_key)
        return bnd_keys

    bnd_ngb_keys_dct = dict_.transform_items_to_values(
        bond_neighborhoods(gra), _neighbor_keys)
    return bnd_ngb_keys_dct


def bond_neighborhoods(gra, stereo=False):
    """ neighborhood subgraphs, by bond
    """
    bnd_keys = list(bond_keys(gra))

    def _neighborhood(bnd_key):
        nbh_bnd_keys = set(filter(lambda x: bnd_key & x, bnd_keys))
        return bond_induced_subgraph(gra, nbh_bnd_keys, stereo=stereo)

    bnd_nbh_dct = dict(zip(bnd_keys, map(_neighborhood, bnd_keys)))
    return bnd_nbh_dct


# # other properties
def terminal_heavy_atom_keys(gra):
    """ terminal heavy atoms, sorted by atom type and hydrogen count
    """
    gra = implicit(gra)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(gra)
    atm_keys = [key for key, ngb_keys in atoms_neighbor_atom_keys(gra).items()
                if len(ngb_keys) == 1]
    atm_keys = sorted(atm_keys, key=atm_imp_hyd_vlc_dct.__getitem__,
                      reverse=True)
    atm_symbs = dict_.values_by_key(atom_symbols(gra), atm_keys)
    srt = automol.formula.argsort_symbols(atm_symbs, symbs_first=('C',))
    atm_keys = tuple(map(atm_keys.__getitem__, srt))
    return atm_keys


def branch(gra, atm_key, bnd_key):
    """ branch extending along `bnd_key` away from `atm_key`
    """
    return bond_induced_subgraph(
        gra, branch_bond_keys(gra, atm_key, bnd_key), stereo=True)


def branch_atom_keys(gra, atm_key, bnd_key):
    """ atom keys for branch extending along `bnd_key` away from `atm_key`
    """
    bnch_atm_keys = atom_keys(branch(gra, atm_key, bnd_key))
    return bnch_atm_keys - {atm_key}


def branch_bond_keys(gra, atm_key, bnd_key):
    """ bond keys for branch extending along `bnd_key` away from `atm_key`
    """

    # bnd_key is the set of atom indices for the bond of interest
    # atm_bnd_keys_dct is a dictionary of atoms that are connected to each atom
    bnd_key = frozenset(bnd_key)
    assert atm_key in bnd_key

    atm_bnd_keys_dct = atoms_bond_keys(gra)

    bnch_bnd_keys = {bnd_key}
    seen_bnd_keys = set()
    excl_bnd_keys = atm_bnd_keys_dct[atm_key] - {bnd_key}

    new_bnd_keys = {bnd_key}

    bnd_ngb_keys_dct = bonds_neighbor_bond_keys(gra)

    while new_bnd_keys:
        new_bnd_ngb_keys = set(
            itertools.chain(
                *dict_.values_by_key(bnd_ngb_keys_dct, new_bnd_keys)))
        bnch_bnd_keys.update(new_bnd_ngb_keys - excl_bnd_keys)
        seen_bnd_keys.update(new_bnd_keys)
        new_bnd_keys = bnch_bnd_keys - seen_bnd_keys

    return frozenset(bnch_bnd_keys)


def is_connected(gra):
    """ is this a connected graph
    """
    return len(connected_components(gra)) == 1


def connected_components(gra):
    """ connected components in the graph
    """
    cmp_gra_atm_keys_lst = connected_components_atom_keys(gra)
    cmp_gras = tuple(subgraph(gra, cmp_gra_atm_keys, stereo=True)
                     for cmp_gra_atm_keys in cmp_gra_atm_keys_lst)
    return cmp_gras


def connected_components_atom_keys(gra):
    """ atom keys for each connected component in the graph
    """
    nxg = _networkx.from_graph(gra)
    cmp_gra_atm_keys_lst = _networkx.connected_component_atom_keys(nxg)
    return cmp_gra_atm_keys_lst


def atom_shortest_paths(gra):
    """ shortest paths between any two atoms in the graph

    :returns: a 2d dictionary keyed by pairs of atoms
    """
    nxg = _networkx.from_graph(gra)
    sp_dct = dict(_networkx.all_pairs_shortest_path(nxg))
    return sp_dct


def shortest_path_between_atoms(gra, key1, key2):
    """ shortest path between a pair of atoms
    """
    return shortest_path_between_groups(gra, [key1], [key2])


def shortest_path_between_groups(gra, keys1, keys2):
    """ shortest path between two groups of atoms

    Returns the atom pair from these groups that are nearest to each other and
    returns the path between them.
    """
    assert not set(keys1) & set(keys2), ("{:s} overlaps with {:s}"
                                         .format(*map(str, [keys1, keys2])))

    sp_dct = atom_shortest_paths(gra)
    keys = None
    for key1 in keys1:
        for key2 in keys2:
            if key2 in sp_dct[key1]:
                if keys is None or len(keys) > len(sp_dct[key1][key2]):
                    keys = sp_dct[key1][key2]

    return keys


def longest_chain(gra):
    """ longest chain in the graph
    """
    atm_keys = atom_keys(gra)

    max_chain = max((atom_longest_chain(gra, atm_key) for atm_key in atm_keys),
                    key=len)
    return max_chain


def atom_longest_chains(gra):
    """ longest chains, by atom
    """
    atm_keys = atom_keys(gra)

    long_chain_dct = {atm_key: atom_longest_chain(gra, atm_key)
                      for atm_key in atm_keys}
    return long_chain_dct


def atom_longest_chain(gra, atm_key):
    """ longest chain for a specific atom
    """
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    chains_lst = []
    if atm_ngb_keys:
        next_chains_lst = [
            [atm_key, atm_ngb_key] for atm_ngb_key in atm_ngb_keys]

        while True:
            chains_lst = next_chains_lst
            next_chains_lst = []
            for chain in chains_lst:
                atm_ngb_keys = atm_ngb_keys_dct[chain[-1]]
                next_atm_keys = sorted(atm_ngb_keys - set(chain))
                for next_atm_key in next_atm_keys:
                    next_chains_lst.append(chain + [next_atm_key])

            if not next_chains_lst:
                break

        max_chain = tuple(chains_lst[0])
    else:
        max_chain = tuple((atm_key,))
    return max_chain


def union(gra1, gra2, check=True):
    """ a union of two graphs
    """
    if check:
        assert not atom_keys(gra1) & atom_keys(gra2)
    atm_dct = {}
    atm_dct.update(atoms(gra1))
    atm_dct.update(atoms(gra2))

    bnd_dct = {}
    bnd_dct.update(bonds(gra1))
    bnd_dct.update(bonds(gra2))
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


def union_from_sequence(gras, check=True):
    """ a union of all parts of a sequence of graphs
    """
    def _union(gra1, gra2):
        return union(gra1, gra2, check=check)

    return tuple(functools.reduce(_union, gras))


def subgraph(gra, atm_keys, stereo=False):
    """ the subgraph induced by a subset of the atoms

    :param gra: the graph
    :param atm_keys: the atom keys to be included in the subgraph
    :param stereo: whether or not to include stereo in the subgraph
    :returns: the subgraph
    """
    atm_keys = set(atm_keys)
    assert atm_keys <= atom_keys(gra)
    bnd_keys = set(filter(lambda x: x <= atm_keys, bond_keys(gra)))
    atm_dct = dict_.by_key(atoms(gra), atm_keys)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    sub = _create.from_atoms_and_bonds(atm_dct, bnd_dct)
    if not stereo:
        sub = without_stereo_parities(sub)
    return sub


def bond_induced_subgraph(gra, bnd_keys, stereo=False):
    """ the subgraph induced by a subset of the bonds
    """
    atm_keys = set(itertools.chain(*bnd_keys))
    bnd_keys = set(bnd_keys)
    assert atm_keys <= atom_keys(gra)
    atm_dct = dict_.by_key(atoms(gra), atm_keys)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    sub = _create.from_atoms_and_bonds(atm_dct, bnd_dct)
    if not stereo:
        sub = without_stereo_parities(sub)
    return sub


# # transformations
def add_atom_implicit_hydrogen_valences(gra, inc_atm_imp_hyd_vlc_dct):
    """ add atom imlicit hydrogen valences

    (increments can be positive or negative)
    """
    atm_keys = list(inc_atm_imp_hyd_vlc_dct.keys())
    atm_imp_hyd_vlcs = numpy.add(
        dict_.values_by_key(atom_implicit_hydrogen_valences(gra), atm_keys),
        dict_.values_by_key(inc_atm_imp_hyd_vlc_dct, atm_keys))
    assert all(atm_imp_hyd_vlc >= 0 for atm_imp_hyd_vlc in atm_imp_hyd_vlcs)
    atm_imp_hyd_vlc_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_imp_hyd_vlcs)), int)
    return set_atom_implicit_hydrogen_valences(gra, atm_imp_hyd_vlc_dct)


def add_atoms(gra, symb_dct, imp_hyd_vlc_dct=None, ste_par_dct=None):
    """ add atoms to this molecular graph, setting their keys
    """
    atm_keys = atom_keys(gra)
    atm_symb_dct = atom_symbols(gra)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(gra)
    atm_ste_par_dct = atom_stereo_parities(gra)

    keys = set(symb_dct.keys())
    imp_hyd_vlc_dct = {} if imp_hyd_vlc_dct is None else imp_hyd_vlc_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct

    assert not keys & atm_keys
    assert set(imp_hyd_vlc_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys

    atm_symb_dct.update(symb_dct)
    atm_imp_hyd_vlc_dct.update(imp_hyd_vlc_dct)
    atm_ste_par_dct.update(ste_par_dct)

    atm_dct = _create.atoms_from_data(
        atom_symbols=atm_symb_dct,
        atom_implicit_hydrogen_valences=atm_imp_hyd_vlc_dct,
        atom_stereo_parities=atm_ste_par_dct)
    bnd_dct = bonds(gra)
    gra = _create.from_atoms_and_bonds(atoms=atm_dct, bonds=bnd_dct)
    return gra


def add_bonded_atom(gra, sym, atm_key, bnd_atm_key=None, imp_hyd_vlc=None,
                    atm_ste_par=None, bnd_ord=None, bnd_ste_par=None):
    """ add a single atom with a bond to an atom already in the graph
    """
    atm_keys = atom_keys(gra)

    bnd_atm_key = max(atm_keys) + 1 if bnd_atm_key is None else bnd_atm_key

    symb_dct = {bnd_atm_key: sym}
    imp_hyd_vlc_dct = ({bnd_atm_key: imp_hyd_vlc}
                       if imp_hyd_vlc is not None else None)
    atm_ste_par_dct = ({bnd_atm_key: atm_ste_par}
                       if atm_ste_par is not None else None)

    gra = add_atoms(gra, symb_dct, imp_hyd_vlc_dct=imp_hyd_vlc_dct,
                    ste_par_dct=atm_ste_par_dct)

    bnd_key = frozenset({bnd_atm_key, atm_key})
    bnd_ord_dct = {bnd_key: bnd_ord} if bnd_ord is not None else None
    bnd_ste_par_dct = ({bnd_key: bnd_ste_par}
                       if bnd_ste_par is not None else None)

    gra = add_bonds(gra, [bnd_key], ord_dct=bnd_ord_dct,
                    ste_par_dct=bnd_ste_par_dct)

    return gra


def insert_bonded_atom(gra, sym, atm_key, bnd_atm_key=None, imp_hyd_vlc=None,
                       atm_ste_par=None, bnd_ord=None, bnd_ste_par=None):
    """ insert a single atom with a bond to an atom already in the graph

    Keys will be standardized upon insertion
    """
    keys = sorted(atom_keys(gra))
    bnd_atm_key_ = max(keys) + 1

    gra = add_bonded_atom(gra, sym, atm_key, bnd_atm_key=bnd_atm_key_,
                          imp_hyd_vlc=imp_hyd_vlc, atm_ste_par=atm_ste_par,
                          bnd_ord=bnd_ord, bnd_ste_par=bnd_ste_par)
    if bnd_atm_key != bnd_atm_key_:
        assert bnd_atm_key in keys
        idx = keys.index(bnd_atm_key)
        key_dct = {}
        key_dct.update({k: k for k in keys[:idx]})
        key_dct[bnd_atm_key_] = bnd_atm_key
        key_dct.update({k: k+1 for k in keys[idx:]})
        gra = relabel(gra, key_dct)

    return gra


def add_dummy_atoms(gra, dummy_key_dct):
    """ add dummy atoms to the graph, with dummy bonds to particular atoms

    :param dummy_key_dct: keys are atoms in the graph on which to place a dummy
        atom; values are the desired keys of the dummy atoms themselves, which
        must not overlap with already existing atoms
    """
    atm_keys = atom_keys(gra)
    assert set(dummy_key_dct.keys()) <= atm_keys, (
        "Keys must be existing atoms in the graph.")
    assert not set(dummy_key_dct.values()) & atm_keys, (
        "Dummy atom keys cannot overlap with existing atoms.")

    for key, dummy_key in sorted(dummy_key_dct.items()):
        gra = add_bonded_atom(gra, 'X', key, bnd_atm_key=dummy_key, bnd_ord=0)

    return gra


def insert_dummy_atoms(gra, dummy_key_dct):
    """ add dummy atoms to the graph, with dummy bonds to particular atoms

    :param dummy_key_dct: keys are atoms in the graph on which to place a dummy
        atom; values are the desired keys of the dummy atoms themselves, which
        must not overlap with already existing atoms
    """
    for key, dummy_key in reversed(sorted(dummy_key_dct.items())):
        gra = insert_bonded_atom(
            gra, 'X', key, bnd_atm_key=dummy_key, bnd_ord=0)

    return gra


def standard_keys_without_dummy_atoms(gra):
    """ remove dummy atoms and standardize keys, returning the dummy key
    dictionary for converting back

    Requires that graph follows z-matrix ordering (this is checked)
    """
    dummy_ngb_key_dct = dummy_atoms_neighbor_atom_key(gra)

    dummy_keys_dct = {}
    last_dummy_key = -1
    decr = 0
    for dummy_key, key in sorted(dummy_ngb_key_dct.items()):
        assert last_dummy_key <= key < dummy_key, (
            "{:d} must follow previous dummy {:d} and preced next dummy {:d}"
            .format(key, last_dummy_key, dummy_key))

        dummy_keys_dct[key-decr] = dummy_key-decr
        gra = remove_atoms(gra, [dummy_key])

        decr += 1
        last_dummy_key = dummy_key-decr

    gra = standard_keys(gra)
    return gra, dummy_keys_dct


def remove_atoms(gra, atm_keys, check=True, stereo=False):
    """ remove atoms from the molecular graph
    """
    all_atm_keys = atom_keys(gra)
    atm_keys = set(atm_keys)

    if check:
        assert atm_keys <= all_atm_keys

    atm_keys_left = all_atm_keys - atm_keys
    return subgraph(gra, atm_keys_left, stereo=stereo)


def add_bonds(gra, keys, ord_dct=None, ste_par_dct=None, check=True):
    """ add bonds to this molecular graph
    """
    bnd_keys = set(bond_keys(gra))
    bnd_ord_dct = bond_orders(gra)
    bnd_ste_par_dct = bond_stereo_parities(gra)

    keys = set(map(frozenset, keys))
    ord_dct = {} if ord_dct is None else ord_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct

    ord_dct = dict_.transform_keys(ord_dct, frozenset)
    ste_par_dct = dict_.transform_keys(ste_par_dct, frozenset)

    if check:
        assert not keys & bnd_keys, (
            '{} and {} have a non-empty intersection'.format(keys, bnd_keys)
        )

    assert set(ord_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys

    bnd_keys.update(keys)
    bnd_ord_dct.update(ord_dct)
    bnd_ste_par_dct.update(ste_par_dct)

    atm_dct = atoms(gra)
    bnd_dct = _create.bonds_from_data(
        bond_keys=bnd_keys, bond_orders=bnd_ord_dct,
        bond_stereo_parities=bnd_ste_par_dct)

    gra = _create.from_atoms_and_bonds(atoms=atm_dct, bonds=bnd_dct)
    return gra


def atom_groups(gra, atm):
    """ return a list of groups off of one atom
    """
    adj_atms = atoms_neighbor_atom_keys(gra)
    keys = []
    for atmi in adj_atms[atm]:
        key = [atm, atmi]
        key.sort()
        key = frozenset(key)
        keys.append(key)
    gras = remove_bonds(gra, keys)
    return connected_components(gras)


def remove_bonds(gra, bnd_keys, check=True):
    """ remove bonds from the molecular graph
    """
    all_bnd_keys = bond_keys(gra)
    bnd_keys = set(map(frozenset, bnd_keys))

    if check:
        assert bnd_keys <= all_bnd_keys

    bnd_keys = all_bnd_keys - bnd_keys
    atm_dct = atoms(gra)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


def without_dummy_atoms(gra):
    """ remove dummy atoms from the graph
    """
    atm_symb_dct = atom_symbols(gra)
    atm_keys = [key for key, sym in atm_symb_dct.items()
                if ptab.to_number(sym)]
    return subgraph(gra, atm_keys, stereo=True)


def without_fractional_bonds(gra):
    """ rounds fractional bonds in the graph
    """
    ord_dct = dict_.transform_values(bond_orders(gra), func=round)
    gra = set_bond_orders(gra, ord_dct)
    return gra


def without_dummy_bonds(gra):
    """ remove 0-order bonds from the graph
    """
    ord_dct = dict_.filter_by_value(bond_orders(gra), func=lambda x: x == 0)
    gra = remove_bonds(gra, ord_dct.keys())
    return gra


def move_idx_to_top(gra, idx1, idx2):
    """ move indexing for atm at idx1 to idx2
    """
    atms, bnds = gra
    newatms = {}
    newbnds = {}
    for key in atms:
        if key == idx1:
            newatms[idx2] = atms[key]
        elif idx2 <= key < idx1:
            newatms[key + 1] = atms[key]
        else:
            newatms[key] = atms[key]
    for key in bnds:
        atm1, atm2 = list(key)
        if atm1 == idx1:
            key1 = idx2
        elif idx2 <= atm1 < idx1:
            key1 = atm1 + 1
        else:
            key1 = atm1
        if atm2 == idx1:
            key2 = idx2
        elif idx2 <= atm2 < idx1:
            key2 = atm2 + 1
        else:
            key2 = atm2
        newkey = [key1, key2]
        newkey.sort()
        newkey = frozenset(newkey)
        newbnds[newkey] = bnds[key]
    return (newatms, newbnds)


# implicit/explicit hydrogen functions
# # atom properties
def atom_explicit_hydrogen_valences(gra):
    """ explicit hydrogen valences, by atom
    """
    return dict_.transform_values(atom_explicit_hydrogen_keys(gra), len)


def atom_explicit_hydrogen_keys(gra):
    """ explicit hydrogen valences, by atom
    """
    exp_hyd_keys = explicit_hydrogen_keys(gra)
    atm_exp_hyd_keys_dct = dict_.transform_values(
        atoms_neighbor_atom_keys(gra), lambda x: x & exp_hyd_keys)
    return atm_exp_hyd_keys_dct


# # other properties
def backbone_keys(gra):
    """ backbone atom keys
    """
    bbn_keys = atom_keys(gra) - explicit_hydrogen_keys(gra)
    return bbn_keys


def explicit_hydrogen_keys(gra):
    """ explicit hydrogen keys (H types: explicit, implicit, backbone)
    """
    hyd_keys = dict_.keys_by_value(atom_symbols(gra), lambda x: x == 'H')
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    def _is_backbone(hyd_key):
        is_h2 = all(ngb_key in hyd_keys and hyd_key < ngb_key
                    for ngb_key in atm_ngb_keys_dct[hyd_key])
        is_multivalent = len(atm_ngb_keys_dct[hyd_key]) > 1
        return is_h2 or is_multivalent

    exp_hyd_keys = frozenset(fmit.filterfalse(_is_backbone, hyd_keys))
    return exp_hyd_keys


# # transformations
def add_atom_explicit_hydrogen_keys(gra, atm_exp_hyd_keys_dct):
    """ add explicit hydrogens by atom
    """
    assert set(atm_exp_hyd_keys_dct.keys()) <= atom_keys(gra), (
        '{} !<= {}'.format(
            set(atm_exp_hyd_keys_dct.keys()), atom_keys(gra))
    )
    for atm_key, atm_exp_hyd_keys in atm_exp_hyd_keys_dct.items():
        assert not set(atm_exp_hyd_keys) & atom_keys(gra)
        atm_exp_hyd_bnd_keys = {frozenset({atm_key, atm_exp_hyd_key})
                                for atm_exp_hyd_key in atm_exp_hyd_keys}
        atm_exp_hyd_symb_dct = dict_.by_key({}, atm_exp_hyd_keys, fill_val='H')
        gra = add_atoms(gra, atm_exp_hyd_symb_dct)
        gra = add_bonds(gra, atm_exp_hyd_bnd_keys)
    return gra


def implicit(gra, atm_keys=None):
    """ make the hydrogens at these atoms implicit
    """
    atm_keys = backbone_keys(gra) if atm_keys is None else atm_keys

    atm_exp_hyd_keys_dct = dict_.by_key(
        atom_explicit_hydrogen_keys(gra), atm_keys)

    inc_imp_hyd_keys_dct = dict_.transform_values(atm_exp_hyd_keys_dct, len)
    gra = add_atom_implicit_hydrogen_valences(gra, inc_imp_hyd_keys_dct)

    exp_hyd_keys = set(itertools.chain(*atm_exp_hyd_keys_dct.values()))
    gra = remove_atoms(gra, exp_hyd_keys)
    return gra


def explicit(gra, atm_keys=None):
    """ make the hydrogens at these atoms explicit
    """
    atm_keys = backbone_keys(gra) if atm_keys is None else atm_keys
    atm_keys = sorted(atm_keys)
    atm_imp_hyd_vlc_dct = dict_.by_key(
        atom_implicit_hydrogen_valences(gra), atm_keys)

    atm_exp_hyd_keys_dct = {}
    next_atm_key = max(atom_keys(gra)) + 1
    for atm_key in atm_keys:
        imp_hyd_vlc = atm_imp_hyd_vlc_dct[atm_key]
        atm_exp_hyd_keys_dct[atm_key] = set(
            range(next_atm_key, next_atm_key+imp_hyd_vlc))
        next_atm_key += imp_hyd_vlc

    gra = set_atom_implicit_hydrogen_valences(
        gra, dict_.by_key({}, atm_keys, fill_val=0))
    gra = add_atom_explicit_hydrogen_keys(gra, atm_exp_hyd_keys_dct)
    return gra


# # comparisons
# def isomorphisms(gra1, gra2):
#     """ isomorphisms between two graphs
#     """
#     assert gra1 == explicit(gra1) and gra2 == explicit(gra2)
#     igr1 = _igraph.from_graph(gra1)
#     igr2 = _igraph.from_graph(gra2)
#     iso_dcts = _igraph.isomorphisms(igr1, igr2)
#     return iso_dcts


def full_isomorphism(gra1, gra2, igraph=False): #True):
    """ full graph isomorphism
    """
    assert gra1 == explicit(gra1) and gra2 == explicit(gra2)
    if igraph:
        igr1 = _igraph.from_graph(gra1)
        igr2 = _igraph.from_graph(gra2)
        iso_dcts = _igraph.isomorphisms(igr1, igr2)
        iso_dct = iso_dcts[0] if iso_dcts else None
    else:
        nxg1 = _networkx.from_graph(gra1)
        nxg2 = _networkx.from_graph(gra2)
        iso_dct = _networkx.isomorphism(nxg1, nxg2)
    return iso_dct


def full_subgraph_isomorphism(gra1, gra2):
    """ gra2 is fully isomorphic to a subgraph of gra1
    """
    assert gra1 == explicit(gra1) and gra2 == explicit(gra2)
    nxg1 = _networkx.from_graph(gra1)
    nxg2 = _networkx.from_graph(gra2)
    iso_dct = _networkx.subgraph_isomorphism(nxg1, nxg2)
    return iso_dct


def backbone_isomorphic(gra1, gra2):
    """ are these molecular graphs backbone isomorphic?
    """
    return backbone_isomorphism(gra1, gra2) is not None


def backbone_isomorphism(gra1, gra2, igraph=False): #True):
    """ graph backbone isomorphism

    for implicit graphs, this is the relabeling of `gra1` to produce `gra2`
    for other graphs, it gives the correspondences between backbone atoms
    """
    gra1 = implicit(gra1)
    gra2 = implicit(gra2)
    if igraph:
        igr1 = _igraph.from_graph(gra1)
        igr2 = _igraph.from_graph(gra2)
        iso_dcts = _igraph.isomorphisms(igr1, igr2)
        iso_dct = iso_dcts[0] if iso_dcts else None
    else:
        nxg1 = _networkx.from_graph(gra1)
        nxg2 = _networkx.from_graph(gra2)
        iso_dct = _networkx.isomorphism(nxg1, nxg2)
    return iso_dct


def backbone_unique(gras):
    """ unique non-isomorphic graphs from a series
    """
    gras = _unique(gras, equiv=backbone_isomorphic)
    return gras


def _unique(itms, equiv):
    """ unique items from a list, according to binary comparison `equiv`
    """
    uniq_itms = []
    for itm in itms:
        if not any(map(functools.partial(equiv, itm), uniq_itms)):
            uniq_itms.append(itm)

    return tuple(uniq_itms)


# chemistry library
VALENCE_DCT = {
    None: 0,
    1: 1,   # H
    2: 2,   # Be
    13: 3,  # B
    14: 4,  # C
    15: 3,  # N
    16: 2,  # O
    17: 1,  # F
    18: 0,  # He
}

LONE_PAIR_COUNTS_DCT = {
    None: 0,
    1: 0,   # H
    2: 0,   # Be
    13: 0,  # B
    14: 0,  # C
    15: 1,  # N
    16: 2,  # O
    17: 3,  # F
    18: 4,  # He
}


# # atom properties
def atom_element_valences(gra):
    """ element valences (# possible single bonds), by atom
    """
    atm_symb_dct = atom_symbols(gra)
    atm_group_idx_dct = dict_.transform_values(atm_symb_dct, ptab.to_group)
    atm_elem_vlc_dct = dict_.transform_values(atm_group_idx_dct,
                                              VALENCE_DCT.__getitem__)
    return atm_elem_vlc_dct


def atom_lone_pair_counts(gra):
    """ lone pair counts, by atom
    """
    atm_symb_dct = atom_symbols(gra)
    atm_group_idx_dct = dict_.transform_values(atm_symb_dct, ptab.to_group)
    atm_lpc_dct = dict_.transform_values(atm_group_idx_dct,
                                         LONE_PAIR_COUNTS_DCT.__getitem__)
    atm_lpc_dct = dict_.transform_values(atm_lpc_dct, int)
    return atm_lpc_dct


def atom_bond_valences(gra, bond_order=True):
    """ bond count (bond valence), by atom
    """
    atm_keys = list(atom_keys(gra))
    gra = explicit(gra)
    if not bond_order:
        gra = without_bond_orders(gra)

    atm_nbhs = dict_.values_by_key(atom_neighborhoods(gra), atm_keys)
    atm_bnd_vlcs = [sum(bond_orders(nbh).values()) for nbh in atm_nbhs]
    atm_bnd_vlc_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_bnd_vlcs)), int)
    return atm_bnd_vlc_dct


def atom_unsaturated_valences(gra, bond_order=True):
    """ unsaturated valences, by atom

    element valences minus bonding valences = pi sites and radical electrons
    """
    atm_keys = list(atom_keys(gra))
    if not bond_order:
        gra = without_bond_orders(gra)

    atm_bnd_vlcs = dict_.values_by_key(atom_bond_valences(gra), atm_keys)
    atm_tot_vlcs = dict_.values_by_key(atom_element_valences(gra), atm_keys)
    atm_rad_vlcs = numpy.subtract(atm_tot_vlcs, atm_bnd_vlcs)
    atm_unsat_vlc_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_rad_vlcs)), int)
    return atm_unsat_vlc_dct


def unsaturated_atom_keys(gra):
    """ keys of unsaturated (radical or pi-bonded) atoms
    """
    atm_unsat_vlc_dct = atom_unsaturated_valences(gra, bond_order=False)
    unsat_atm_keys = frozenset(dict_.keys_by_value(atm_unsat_vlc_dct, bool))
    return unsat_atm_keys


# # other properties
def maximum_spin_multiplicity(gra, bond_order=True):
    """ the highest possible spin multiplicity for this molecular graph
    """
    atm_rad_vlc_dct = atom_unsaturated_valences(gra, bond_order=bond_order)
    return sum(atm_rad_vlc_dct.values()) + 1


def possible_spin_multiplicities(gra, bond_order=True):
    """ possible spin multiplicities for this molecular graph
    """
    mult_max = maximum_spin_multiplicity(gra, bond_order=bond_order)
    mult_min = 2 if mult_max % 2 == 0 else 1
    mults = tuple(range(mult_min, mult_max+1, 2))
    return mults


# miscellaneous
# # bond properties
def bond_symmetry_numbers(gra, frm_bnd_key=None, brk_bnd_key=None):
    """ symmetry numbers, by bond

    the (approximate) symmetry number of the torsional potential for this bond,
    based on the hydrogen counts for each atom
    It is reduced to 1 if one of the H atoms in the torsional bond is a
    neighbor to the special bonding atom (the atom that is being transferred)
    """
    imp_gra = implicit(gra)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(imp_gra)

    bnd_keys = bond_keys(imp_gra)

    tfr_atm = None
    if frm_bnd_key and brk_bnd_key:
        for atm_f in list(frm_bnd_key):
            for atm_b in list(brk_bnd_key):
                if atm_f == atm_b:
                    tfr_atm = atm_f

        if tfr_atm:
            neighbor_dct = atoms_neighbor_atom_keys(gra)
            nei_tfr = neighbor_dct[tfr_atm]

            atms = gra[0]
            all_hyds = []
            for atm in atms:
                if atms[atm][0] == 'H':
                    all_hyds.append(atm)
        else:
            nei_tfr = {}

    bnd_symb_num_dct = {}
    bnd_symb_nums = []
    for bnd_key in bnd_keys:
        bnd_sym = 1
        vlc = max(map(atm_imp_hyd_vlc_dct.__getitem__, bnd_key))
        if vlc == 3:
            bnd_sym = 3
            if tfr_atm:
                for atm in nei_tfr:
                    nei_s = neighbor_dct[atm]
                    h_nei = 0
                    for nei in nei_s:
                        if nei in all_hyds:
                            h_nei += 1
                    if h_nei == 3:
                        bnd_sym = 1
        bnd_symb_nums.append(bnd_sym)

    bnd_symb_num_dct = dict(zip(bnd_keys, bnd_symb_nums))

    # fill in the rest of the bonds for completeness
    bnd_symb_num_dct = dict_.by_key(
        bnd_symb_num_dct, bond_keys(gra), fill_val=1)

    return bnd_symb_num_dct
