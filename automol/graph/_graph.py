""" molecular graph
"""

import operator
import itertools
import functools
import collections.abc
import numpy
from phydat import ptab
import automol.formula
from automol import util
from automol.util import dict_
from automol.graph._graph_dep import atoms
from automol.graph._graph_dep import bonds
from automol.graph._graph_dep import atom_keys
from automol.graph._graph_dep import bond_keys
from automol.graph._graph_dep import atom_symbols
from automol.graph._graph_dep import atom_implicit_hydrogen_valences
from automol.graph._graph_dep import set_atom_symbols
from automol.graph._graph_dep import relabel
from automol.graph._graph_dep import without_stereo_parities
from automol.graph._graph_dep import atom_explicit_hydrogen_keys
from automol.graph._graph_dep import without_dummy_atoms
from automol.graph._graph_dep import add_bonds
from automol.graph._graph_dep import remove_atoms
from automol.graph._graph_dep import remove_bonds
from automol.graph._graph_dep import implicit
from automol.graph._graph_dep import explicit
from automol.graph._graph_dep import add_atoms
from automol.graph._graph_dep import atom_unsaturated_valences
from automol.graph._graph_dep import maximum_spin_multiplicity
from automol.graph._graph_dep import atom_neighborhoods
from automol.graph._graph_dep import bond_neighborhoods
from automol.graph._graph_dep import atoms_neighbor_atom_keys
from automol.graph._graph_dep import subgraph
from automol.graph._graph_dep import bond_induced_subgraph
from automol.graph._graph_dep import dummy_atoms_neighbor_atom_key
from automol.graph._embed_dep import atom_shortest_paths
from automol.graph._embed_dep import union
from automol.graph._embed_dep import backbone_isomorphic
# graphbase
from automol.graph import _networkx
from automol.graph import _igraph


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

    gras = tuple(relabel(gra, atm_key_dct)
                 for gra, atm_key_dct in zip(gras, atm_key_dcts))
    atm_key_dcts = tuple(atm_key_dcts)

    return gras, atm_key_dcts


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


def frozen(gra):
    """ hashable, sortable, immutable container of graph data
    """
    atm_keys = sorted(atom_keys(gra))
    bnd_keys = sorted(bond_keys(gra), key=sorted)

    # make it sortable by replacing Nones with -infinity
    atm_vals = numpy.array(dict_.values_by_key(atoms(gra), atm_keys),
                           dtype=numpy.object)
    bnd_vals = numpy.array(dict_.values_by_key(bonds(gra), bnd_keys),
                           dtype=numpy.object)
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


def is_branched(gra):
    """ determine is the molecule has a branched chain
    """
    _is_branched = False
    gra = implicit(gra)
    chain_length = len(longest_chain(gra))
    natoms = atom_count(gra, with_implicit=False)

    if natoms != chain_length:
        _is_branched = True
    return _is_branched


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


def union_from_sequence(gras, check=True):
    """ a union of all parts of a sequence of graphs
    """
    def _union(gra1, gra2):
        return union(gra1, gra2, check=check)

    return tuple(functools.reduce(_union, gras))


# # transformations
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
        # If the dummy key comes first, the anchor key will be offset by 1
        if dummy_key < key:
            key -= 1

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
        assert last_dummy_key <= key, (
            "{:d} must follow previous dummy {:d}"
            .format(key, last_dummy_key))

        dummy_keys_dct[key-decr] = dummy_key-decr
        gra = remove_atoms(gra, [dummy_key])

        decr += 1
        last_dummy_key = dummy_key-decr

    gra = standard_keys(gra)
    return gra, dummy_keys_dct


def atom_groups(gra, atm, stereo=False):
    """ return a list of groups off of one atom
    """
    if not stereo:
        gra = without_stereo_parities(gra)
    adj_atms = atoms_neighbor_atom_keys(gra)
    keys = []
    for atmi in adj_atms[atm]:
        key = [atm, atmi]
        key.sort()
        key = frozenset(key)
        keys.append(key)
    gras = remove_bonds(gra, keys)
    return connected_components(gras)


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


# # other properties
# # transformations
# # comparisons
# def isomorphisms(gra1, gra2):
#     """ isomorphisms between two graphs
#     """
#     assert gra1 == explicit(gra1) and gra2 == explicit(gra2)
#     igr1 = _igraph.from_graph(gra1)
#     igr2 = _igraph.from_graph(gra2)
#     iso_dcts = _igraph.isomorphisms(igr1, igr2)
#     return iso_dcts


def isomorphism(gra1, gra2, backbone_only=False, stereo=True, dummy=True):
    """ Obtain an isomorphism between two graphs

    This should eventually replace the other isomorphism functions.

    :param backbone_only: Compare backbone atoms only?
    :type backbone_only: bool
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: The isomorphism mapping `gra1` onto `gra2`
    :rtype: dict
    """
    if backbone_only:
        gra1 = implicit(gra1)
        gra2 = implicit(gra2)

    if not stereo:
        gra1 = without_stereo_parities(gra1)
        gra2 = without_stereo_parities(gra2)

    if not dummy:
        gra1 = without_dummy_atoms(gra1)
        gra2 = without_dummy_atoms(gra2)

    return _isomorphism(gra1, gra2)


def _isomorphism(gra1, gra2, igraph=False):
    """
    """
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


def equivalent_atoms(gra, atm_key, stereo=True, dummy=True):
    """ Identify sets of isomorphically equivalent atoms

    Two atoms are equivalent if they transform into each other under an
    automorphism

    :param gra: A graph
    :param atm_key: An atom key for the graph
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: Keys to equivalent atoms
    :rtype: frozenset
    """
    assert atm_key in atom_keys(gra), (
        "{} not in {}".format(atm_key, atom_keys(gra)))

    atm_symb_dct = atom_symbols(gra)
    atm_ngbs_dct = atoms_neighbor_atom_keys(gra)

    def _neighbor_symbols(key):
        return sorted(map(atm_symb_dct.__getitem__, atm_ngbs_dct[key]))

    # 1. Find atoms with the same symbols
    atm_symb = atm_symb_dct[atm_key]
    cand_keys = atom_keys(gra, sym=atm_symb)

    # 2. Of those, find atoms with the same neighboring atom types
    atm_ngb_symbs = _neighbor_symbols(atm_key)
    cand_keys = [k for k in cand_keys
                 if _neighbor_symbols(k) == atm_ngb_symbs]

    # 3. Find the equivalent atoms from the list of candidates.
    # Strategy: Change the atom symbol to 'Ts' and check for isomorphism.
    # Assumes none of the compounds have element 117.
    atm_keys = []
    for key in cand_keys:
        if are_equivalent_atoms(gra, atm_key, key, stereo=stereo, dummy=dummy):
            atm_keys.append(key)

    return frozenset(atm_keys)


def equivalent_bonds(gra, bnd_key, stereo=True, dummy=True):
    """ Identify sets of isomorphically equivalent bonds

    Two bonds are equivalent if they transform into each other under an
    automorphism

    :param gra: A graph
    :param bnd_key: An bond key for the graph, which may be sorted or unsorted
    :param backbone_only: Compare backbone atoms only?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: Keys to equivalent bonds
    :rtype: frozenset
    """
    bnd_key = tuple(bnd_key)
    bnd_keys = list(map(tuple, map(sorted, bond_keys(gra))))
    bnd_keys += list(map(tuple, map(reversed, bnd_keys)))
    assert bnd_key in bnd_keys, "{} not in {}".format(bnd_key, bnd_keys)

    atm_symb_dct = atom_symbols(gra)
    atm_ngbs_dct = atoms_neighbor_atom_keys(gra)

    def _symbols(bnd_key):
        return list(map(atm_symb_dct.__getitem__, bnd_key))

    def _neighbor_symbols(bnd_key):
        key1, key2 = bnd_key
        nsymbs1 = sorted(map(atm_symb_dct.__getitem__, atm_ngbs_dct[key1]))
        nsymbs2 = sorted(map(atm_symb_dct.__getitem__, atm_ngbs_dct[key2]))
        return nsymbs1, nsymbs2

    # 1. Find bonds with the same atom types
    bnd_symbs = _symbols(bnd_key)
    cand_keys = [k for k in bnd_keys if _symbols(k) == bnd_symbs]

    # 2. Of those, find bonds with the same neighboring atom types
    bnd_ngb_symbs = _neighbor_symbols(bnd_key)
    cand_keys = [k for k in cand_keys
                 if _neighbor_symbols(k) == bnd_ngb_symbs]

    # 3. Find the equivalent bonds from the list of candidates.
    # Strategy: Change the atom symbols to 'Lv' and 'Ts' and check for
    # isomorphism.  Assumes none of the compounds have element 116 or 117.
    bnd_keys = []
    for key in cand_keys:
        if are_equivalent_bonds(gra, bnd_key, key, stereo=stereo, dummy=dummy):
            bnd_keys.append(key)

    return frozenset(bnd_keys)


def are_equivalent_atoms(gra, atm1_key, atm2_key, stereo=True, dummy=True):
    """ Determine whether two atoms are isomorphically equivalent.

    Two atoms are equivalent if they transform into each other under an
    automorphism

    :param gra: A graph
    :param atm1_key: The first atom
    :type atm1_key: int
    :param atm2_key: The first atom
    :type atm2_key: int
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: True if the atoms are equivalent, False otherwise
    :rtype: bool
    """
    gra1 = set_atom_symbols(gra, {atm1_key: 'Ts'})
    gra2 = set_atom_symbols(gra, {atm2_key: 'Ts'})
    are_equiv = bool(isomorphism(gra1, gra2, stereo=stereo, dummy=dummy))
    return are_equiv


def are_equivalent_bonds(gra, bnd1_key, bnd2_key, stereo=True, dummy=True):
    """ Determine whether two bonds are isomorphically equivalent.

    Two bonds are equivalent if they transform into each other under an
    automorphism

    :param gra: A graph
    :param bnd1_key: The first atom
    :type bnd1_key: int
    :param bnd2_key: The first atom
    :type bnd2_key: int
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: True if the atoms are equivalent, False otherwise
    :rtype: bool
    """
    order_matters = (isinstance(bnd1_key, collections.abc.Sequence) and
                     isinstance(bnd2_key, collections.abc.Sequence))

    bnd1_key = list(bnd1_key)
    bnd2_key = list(bnd2_key)
    gra1 = set_atom_symbols(gra, {bnd1_key[0]: 'Lv', bnd1_key[1]: 'Ts'})
    gra2 = set_atom_symbols(gra, {bnd2_key[0]: 'Lv', bnd2_key[1]: 'Ts'})
    are_equiv = bool(isomorphism(gra1, gra2, stereo=stereo, dummy=dummy))

    # If order doesn't matter, check swap atoms and check again
    if not order_matters:
        gra2 = set_atom_symbols(gra, {bnd2_key[1]: 'Lv', bnd2_key[0]: 'Ts'})
        are_equiv |= bool(isomorphism(gra1, gra2, stereo=stereo, dummy=dummy))

    return are_equiv


def atom_equivalence_class_reps(gra, atm_keys=None, stereo=True, dummy=True):
    """ Find equivalence class representatives for atoms in the class

    This function identifies isomorphically unique atoms, which do not
    transform into each other under an automorphism.

    Optionally, a subset of atoms can be passed in to consider class
    representatives from within that list.

    :param gra: A graph
    :param atm_keys: An optional list of atom keys from which to determine
        equivalence class representatives. If None, the full set of atom keys
        will be used.
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: The list of equivalence class reprentatives/unique atoms.
    :rtype: frozenset[int]
    """
    atm_keys = atom_keys(gra) if atm_keys is None else atm_keys

    def _equiv(atm1_key, atm2_key):
        return are_equivalent_atoms(gra, atm1_key, atm2_key, stereo=stereo,
                                    dummy=dummy)

    eq_classes = util.equivalence_partition(atm_keys, _equiv)
    class_reps = frozenset(next(iter(c)) for c in eq_classes)
    return class_reps


def bond_equivalence_class_reps(gra, bnd_keys=None, stereo=True, dummy=True):
    """ Find equivalence class representatives for bonds in the class

    This function identifies isomorphically unique bonds, which do not
    transform into each other under an automorphism.

    Optionally, a subset of bonds can be passed in to consider class
    representatives from within that list.

    :param gra: A graph
    :param bnd_keys: An optional list of bond keys from which to determine
        equivalence class representatives. If None, the full set of bond keys
        will be used.
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: The list of equivalence class reprentatives/unique bonds.
    :rtype: frozenset[int]
    """
    bnd_keys = bond_keys(gra) if bnd_keys is None else bnd_keys

    def _equiv(bnd1_key, bnd2_key):
        return are_equivalent_bonds(gra, bnd1_key, bnd2_key, stereo=stereo,
                                    dummy=dummy)

    eq_classes = util.equivalence_partition(bnd_keys, _equiv)
    class_reps = frozenset(next(iter(c)) for c in eq_classes)
    return class_reps


# def full_isomorphism(gra1, gra2, igraph=True):
def full_isomorphism(gra1, gra2, igraph=False):
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


# # atom properties
def unsaturated_atom_keys(gra):
    """ keys of unsaturated (radical or pi-bonded) atoms
    """
    atm_unsat_vlc_dct = atom_unsaturated_valences(gra, bond_order=False)
    unsat_atm_keys = frozenset(dict_.keys_by_value(atm_unsat_vlc_dct, bool))
    return unsat_atm_keys


# # other properties
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


if __name__ == '__main__':
    GRA = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
            3: ('C', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
            6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None),
            9: ('H', 0, None), 10: ('H', 0, None), 11: ('H', 0, None),
            12: ('H', 0, None), 13: ('H', 0, None)},
           {frozenset({0, 3}): (1, None), frozenset({0, 4}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({0, 6}): (1, None),
            frozenset({1, 3}): (1, None), frozenset({1, 7}): (1, None),
            frozenset({8, 1}): (1, None), frozenset({1, 9}): (1, None),
            frozenset({2, 3}): (1, None), frozenset({2, 10}): (1, None),
            frozenset({2, 11}): (1, None), frozenset({2, 12}): (1, None),
            frozenset({3, 13}): (1, None)})
    print(equivalent_bonds(GRA, [0, 3]))
    # GRA1 = ({0: ('Z', 3, None), 1: ('C', 3, None), 2: ('C', 2, None),
    #          3: ('C', 2, None)},
    #         {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
    #          frozenset({2, 3}): (1, None)})
    # GRA2 = ({0: ('C', 3, None), 1: ('Z', 3, None), 2: ('C', 2, None),
    #          3: ('C', 2, None)},
    #         {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
    #          frozenset({2, 3}): (1, None)})
    # print(isomorphism(GRA1, GRA2))
