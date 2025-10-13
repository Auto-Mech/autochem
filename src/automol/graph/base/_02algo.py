"""algorithm functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import collections.abc
import functools
import itertools
import numbers
import operator
from collections import defaultdict
from typing import Any, Dict, List, Tuple

import more_itertools as mit
import networkx

from ... import util
from . import _01networkx
from ._00core import (
    AtomKeys,
    add_bonds,
    atom_bond_keys,
    atom_count,
    atom_keys,
    atom_neighbor_atom_keys,
    atom_symbols,
    atoms_bond_keys,
    atoms_neighbor_atom_keys,
    bond_induced_subgraph,
    bond_keys,
    bond_orders,
    clear_unstable_bond_orders,
    frozen,
    implicit,
    remove_bonds,
    set_atom_symbols,
    string,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reactants_graph_without_stereo,
    ts_reacting_bond_keys,
    union,
    union_from_sequence,
    without_dummy_atoms,
    without_stereo,
)
from ._00core import (
    subgraph as subgraph_,
)


# # isomorphisms and equivalence
def isomorphism(
    gra1,
    gra2,
    backbone_only=False,
    stereo=True,
    dummy=True,
    subgraph=False,
    unstable=True,
):
    """Obtain an isomorphism between two graphs.

    :param backbone_only: Compare backbone atoms only?
    :type backbone_only: bool
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :param subgraph: Test whether `gra2` is isomorphic to a subgraph of `gra1`?
    :type subgraph: bool
    :param unstable: Whether to include unstable product bonds in the isomorphism
    :returns: The isomorphism mapping `gra1` onto `gra2`
    :rtype: dict
    """
    if backbone_only:
        gra1 = implicit(gra1)
        gra2 = implicit(gra2)

    if not stereo:
        gra1 = without_stereo(gra1)
        gra2 = without_stereo(gra2)

    if not dummy:
        gra1 = without_dummy_atoms(gra1)
        gra2 = without_dummy_atoms(gra2)

    if not unstable:
        gra1 = clear_unstable_bond_orders(gra1)
        gra2 = clear_unstable_bond_orders(gra2)

    nxg1 = _01networkx.from_graph(gra1)
    nxg2 = _01networkx.from_graph(gra2)
    if subgraph:
        iso_dct = _01networkx.subgraph_isomorphism(nxg1, nxg2)
    else:
        iso_dct = _01networkx.isomorphism(nxg1, nxg2)
    return iso_dct


def isomorphic(gra1, gra2, backbone_only=False, stereo=True, dummy=True):
    """Determine whether two graphs are isomorphic

    :param backbone_only: Compare backbone atoms only?
    :type backbone_only: bool
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: The isomorphism mapping `gra1` onto `gra2`
    :rtype: dict
    """
    return (
        isomorphism(gra1, gra2, backbone_only=backbone_only, stereo=stereo, dummy=dummy)
        is not None
    )


def unique(gras, backbone_only=False, stereo=True, dummy=True) -> tuple[Any, ...]:
    """Get the subset of non-isomorphic graphs from a list with redundancies.

    :param backbone_only: Compare backbone atoms only?
    :param stereo: Consider stereo?
    :param dummy: Consider dummy atoms?
    :returns: Unique graphs
    """
    gras_with_counts = unique_with_counts(
        gras, backbone_only=backbone_only, stereo=stereo, dummy=dummy
    )
    if gras_with_counts:
        gras, *_ = zip(*gras_with_counts, strict=True)
        return gras

    return ()


def unique_with_counts(
    gras, backbone_only=False, stereo=True, dummy=True
) -> list[tuple[Any, int]]:
    """Get the subset of non-isomorphic graphs from a list with redundancies.

    :param backbone_only: Compare backbone atoms only?
    :param stereo: Consider stereo?
    :param dummy: Consider dummy atoms?
    :returns: Unique graphs along with their counts.
    """

    def _equiv(gra1, gra2):
        return isomorphic(
            gra1, gra2, backbone_only=backbone_only, stereo=stereo, dummy=dummy
        )

    gras_with_counts = _unique_with_counts(gras, equiv=_equiv)
    return gras_with_counts


def _unique_with_counts(items, equiv) -> list[tuple[Any, int]]:
    """Get unique items from a list, according to binary comparison `equiv`."""
    items0 = items
    items = []
    counts = []
    for item in items0:
        idx = next((ix for ix, it in enumerate(items) if equiv(item, it)), None)
        if idx is not None:
            counts[idx] += 1
        else:
            items.append(item)
            counts.append(1)

    return list(zip(items, counts, strict=True))


def sequence_isomorphism(gras1, gras2, backbone_only=False, stereo=True, dummy=True):
    """Obtain an isomorphism between two sequences of graphs

    :param backbone_only: Compare backbone atoms only?
    :type backbone_only: bool
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: First, a sequence of integers mapping each graph in `gras1` to an
        isomorph in `gras2`, followed by a sequence of dictionaries mapping
        their atoms onto each other.
    :rtype: tuple[int], tuple[dict]
    """
    gras1_pool = dict(enumerate(gras1))
    order = []
    iso_dcts = []
    for gra2 in gras2:
        found_match = False
        for idx, gra1 in gras1_pool.items():
            iso_dct = isomorphism(
                gra2, gra1, backbone_only=backbone_only, stereo=stereo, dummy=dummy
            )
            if iso_dct is not None:
                found_match = True
                order.append(idx)
                iso_dcts.append(iso_dct)
                gras1_pool.pop(idx)
                break

        if not found_match:
            break

    if not found_match:
        order = None
        iso_dcts = None
    else:
        order = tuple(order)
        iso_dcts = tuple(iso_dcts)

    return order, iso_dcts


def subgraph_isomorphism(gra1, gra2, backbone_only=False, stereo=True, dummy=True):
    """Test whether a graph contains a subgraph which is isomorphic to another

    :param backbone_only: Compare backbone atoms only?
    :type backbone_only: bool
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: The isomorphism mapping a subgraph of `gra1` onto `gra2`
    :rtype: dict
    """
    return isomorphism(
        gra1,
        gra2,
        backbone_only=backbone_only,
        stereo=stereo,
        dummy=dummy,
        subgraph=True,
    )


def equivalent_atoms(gra, atm_key, stereo=True, dummy=True):
    """Identify sets of isomorphically equivalent atoms

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
    assert atm_key in atom_keys(gra), f"{atm_key} not in {atom_keys(gra)}"

    atm_symb_dct = atom_symbols(gra)
    atm_ngbs_dct = atoms_neighbor_atom_keys(gra)

    def _neighbor_symbols(key):
        return sorted(map(atm_symb_dct.__getitem__, atm_ngbs_dct[key]))

    # 1. Find atoms with the same symbols
    atm_symb = atm_symb_dct[atm_key]
    cand_keys = atom_keys(gra, symb=atm_symb)

    # 2. Of those, find atoms with the same neighboring atom types
    atm_ngb_symbs = _neighbor_symbols(atm_key)
    cand_keys = [k for k in cand_keys if _neighbor_symbols(k) == atm_ngb_symbs]

    # 3. Find the equivalent atoms from the list of candidates.
    # Strategy: Change the atom symbol to 'Ts' and check for isomorphism.
    # Assumes none of the compounds have element 117.
    atm_keys = []
    for key in cand_keys:
        if are_equivalent_atoms(gra, atm_key, key, stereo=stereo, dummy=dummy):
            atm_keys.append(key)

    return frozenset(atm_keys)


def equivalent_bonds(gra, bnd_key, stereo=True, dummy=True):
    """Identify sets of isomorphically equivalent bonds

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
    assert bnd_key in bnd_keys, f"{bnd_key} not in {bnd_keys}"

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
    cand_keys = [k for k in cand_keys if _neighbor_symbols(k) == bnd_ngb_symbs]

    # 3. Find the equivalent bonds from the list of candidates.
    # Strategy: Change the atom symbols to 'Lv' and 'Ts' and check for
    # isomorphism.  Assumes none of the compounds have element 116 or 117.
    bnd_keys = []
    for key in cand_keys:
        if are_equivalent_bonds(gra, bnd_key, key, stereo=stereo, dummy=dummy):
            bnd_keys.append(key)

    return frozenset(bnd_keys)


def are_equivalent_atoms(gra, atm1_key, atm2_key, stereo=True, dummy=True):
    """Determine whether two atoms are isomorphically equivalent.

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
    gra1 = set_atom_symbols(gra, {atm1_key: "Ts"})
    gra2 = set_atom_symbols(gra, {atm2_key: "Ts"})
    are_equiv = bool(isomorphism(gra1, gra2, stereo=stereo, dummy=dummy))
    return are_equiv


def are_equivalent_bonds(gra, bnd1_key, bnd2_key, stereo=True, dummy=True):
    """Determine whether two bonds are isomorphically equivalent.

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
    order_matters = isinstance(bnd1_key, collections.abc.Sequence) and isinstance(
        bnd2_key, collections.abc.Sequence
    )

    bnd1_key = list(bnd1_key)
    bnd2_key = list(bnd2_key)
    gra1 = set_atom_symbols(gra, {bnd1_key[0]: "Lv", bnd1_key[1]: "Ts"})
    gra2 = set_atom_symbols(gra, {bnd2_key[0]: "Lv", bnd2_key[1]: "Ts"})
    are_equiv = bool(isomorphism(gra1, gra2, stereo=stereo, dummy=dummy))

    # If order doesn't matter, check swap atoms and check again
    if not order_matters:
        gra2 = set_atom_symbols(gra, {bnd2_key[1]: "Lv", bnd2_key[0]: "Ts"})
        are_equiv |= bool(isomorphism(gra1, gra2, stereo=stereo, dummy=dummy))

    return are_equiv


def atom_equivalence_class_reps(gra, atm_keys=None, stereo=True, dummy=True, symb=None):
    """Identify isomorphically unique atoms, which do not transform into each
    other by an automorphism

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
    :param symb: An atomic symbol, restricting the function to atoms of a
        particular type.
    :type symb: str
    :returns: The list of equivalence class reprentatives/unique atoms.
    :rtype: frozenset[int]
    """
    atm_keys = atom_keys(gra, symb=symb) if atm_keys is None else atm_keys

    def _equiv(atm1_key, atm2_key):
        return are_equivalent_atoms(gra, atm1_key, atm2_key, stereo=stereo, dummy=dummy)

    eq_classes = util.equivalence_partition(atm_keys, _equiv)
    class_reps = frozenset(next(iter(c)) for c in eq_classes)
    return class_reps


def bond_equivalence_class_reps(gra, bnd_keys=None, stereo=True, dummy=True):
    """Identify isomorphically unique bonds, which do not transform into each
    other by an automorphism

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
        return are_equivalent_bonds(gra, bnd1_key, bnd2_key, stereo=stereo, dummy=dummy)

    eq_classes = util.equivalence_partition(bnd_keys, _equiv)
    class_reps = frozenset(next(iter(c)) for c in eq_classes)
    return class_reps


# # algorithms
def connected_components(gra, stereo=True):
    """connected components in the graph"""
    cmp_gra_atm_keys_lst = connected_components_atom_keys(gra)
    cmp_gras = tuple(
        subgraph_(gra, cmp_gra_atm_keys, stereo=stereo)
        for cmp_gra_atm_keys in cmp_gra_atm_keys_lst
    )
    return cmp_gras


def connected_components_atom_keys(gra):
    """atom keys for each connected component in the graph"""
    nxg = _01networkx.from_graph(gra)
    cmp_gra_atm_keys_lst = _01networkx.connected_component_atom_keys(nxg)
    return cmp_gra_atm_keys_lst


def is_connected(gra):
    """is this a connected graph"""
    return len(connected_components(gra)) == 1


def dfs_(gra, key) -> list[tuple[int, int]]:
    """Get the result of a depth-first search.

    :param gra: A graph
    :param key: The starting atom key
    :return: The depth-first search, encoded as edge tuples in visitation order
    """
    nxg = _01networkx.from_graph(gra)
    dfs = list(networkx.dfs_edges(nxg, source=key, sort_neighbors=sorted))
    return dfs


def dfs_atom_keys(dfs: list[tuple[int, int]]) -> list[int]:
    """Get the atom keys in order from a depth-first search.

    :param dfs: The depth-first search, encoded as edge tuples in visitation order
    :return: The atom keys in visitation order
    """
    return list(mit.unique_justseen(itertools.chain(*dfs)))


def dfs_bond_keys(dfs: list[tuple[int, int]]) -> list[frozenset[int]]:
    """Get the bond keys in order from a depth-first search.

    :param dfs: The depth-first search, encoded as edge tuples in visitation order
    :return: The bond keys in visitation order
    """
    return list(map(frozenset, dfs))


def dfs_children(dfs: list[tuple[int, int]]) -> dict[int, list[int]]:
    """Get a dictionary of children from a depth-first search.

    :param dfs: The edges visited by the DFS search
    :return: A dictionary mapping parents onto their children, in order
    """
    child_dct = defaultdict(list)
    for pkey, ckey in dfs:
        child_dct[pkey].append(ckey)
    return dict(child_dct)


def dfs_parents(dfs: list[tuple[int, int]]) -> dict[int, list[int]]:
    """Get a dictionary of parents from a depth-first search.

    :param dfs: The edges visited by the DFS search
    :return: A dictionary mapping children onto their parents, in order
    """
    return {c: p for p, c in dfs}


def dfs_missing_bond_keys(gra, dfs: list[tuple[int, int]]) -> frozenset[frozenset[int]]:
    """Get the (ring) bond keys missed by a depth-first search.

    :param gra: A graph
    :param dfs: The depth-first search, encoded as edge tuples in visitation order
    :return: The bond keys missed by the search
    """
    return bond_keys(gra) - frozenset(dfs_bond_keys(dfs))


def dfs_missing_children(gra, dfs: list[tuple[int, int]]) -> dict[int, list[int]]:
    """Get a dictionary of children which are missing from a depth-first search.

    These are the children that are necessary to reconnect missing ring bonds.

    :param gra: A graph
    :param dfs: The depth-first search, encoded as edge tuples in visitation order
    :return: A dictionary mapping parents onto their missing children, in order
    """
    dfs_keys = dfs_atom_keys(dfs)
    miss_bkeys = dfs_missing_bond_keys(gra, dfs)
    miss_child_dct = defaultdict(list)
    for bkey in miss_bkeys:
        pkey, ckey = sorted(bkey, key=dfs_keys.index)
        miss_child_dct[ckey].append(pkey)
    miss_child_dct = util.dict_.transform_values(miss_child_dct, sorted)
    return miss_child_dct


def atom_shortest_paths(gra):
    """shortest paths between any two atoms in the graph

    :returns: a 2d dictionary keyed by pairs of atoms
    """
    nxg = _01networkx.from_graph(gra)
    sp_dct = dict(_01networkx.all_pairs_shortest_path(nxg))
    return sp_dct


def shortest_path_between_atoms(gra, key1, key2):
    """shortest path between a pair of atoms"""
    return shortest_path_between_groups(gra, [key1], [key2])


def simple_paths_between_atoms(gra, key1, key2, overlap: bool = False):
    """shortest paths between any two atoms in the graph

    :param overlap: Whether to include overlapping paths. If not, the shortest one is
        retained.
    :returns: a 2d dictionary keyed by pairs of atoms
    """
    nxg = _01networkx.from_graph(gra)
    return _01networkx.simple_paths(nxg, key1, key2, overlap=overlap)


def shortest_path_between_groups(gra, keys1, keys2):
    """shortest path between two groups of atoms

    Returns the atom pair from these groups that are nearest to each other
    and returns the path between them.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param keys1: Atom keys
    :param keys2: Atom keys to determine path to
    """
    sp_dct = atom_shortest_paths(gra)
    keys = None
    if set(keys1) & set(keys2):
        keys = ()
    else:
        for key1 in keys1:
            for key2 in keys2:
                if key2 in sp_dct[key1]:
                    if keys is None or len(keys) > len(sp_dct[key1][key2]):
                        keys = sp_dct[key1][key2]

    return keys


def atom_longest_chains(gra):
    """longest chains, by atom"""
    atm_keys = atom_keys(gra)

    long_chain_dct = {atm_key: atom_longest_chain(gra, atm_key) for atm_key in atm_keys}
    return long_chain_dct


def atom_longest_chain(gra, atm_key):
    """longest chain for a specific atom"""
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    chains_lst = []
    if atm_ngb_keys:
        next_chains_lst = [[atm_key, atm_ngb_key] for atm_ngb_key in atm_ngb_keys]

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


def longest_chain(gra):
    """longest chain in the graph"""
    atm_keys = atom_keys(gra)

    max_chain = max((atom_longest_chain(gra, atm_key) for atm_key in atm_keys), key=len)
    return max_chain


def weighted_maximal_matching(gra, bnd_weight_dct=None):
    """Calculates a weighted maximal matching of a graph

    That is, a set of edges that covers the graph as much as possible
    """
    edge_attrib_dct = None if bnd_weight_dct is None else {"weight": bnd_weight_dct}
    nxg = _01networkx.from_graph(gra, edge_attrib_dct=edge_attrib_dct)
    bnd_keys = _01networkx.weighted_maximal_matching(nxg, edge_attrib_name="weight")
    return bnd_keys


# # branches and groups
def ring_atom_chirality(gra, atm, ring_atms, stereo=False):
    """is this ring atom a chiral center?"""
    if not stereo:
        gra = without_stereo(gra)
    adj_atms = atoms_neighbor_atom_keys(gra)
    keys = []
    for atmi in adj_atms[atm]:
        key = [atm, atmi]
        key.sort()
        key = frozenset(key)
        keys.append(key)
        if atmi in ring_atms:
            for atmj in adj_atms[atmi]:
                if atmj in ring_atms:
                    key = [atmj, atmi]
                    key.sort()
                    key = frozenset(key)
                    keys.append(key)
    gras = remove_bonds(gra, keys)
    cgras = connected_components(gras)
    ret_gras = []
    for gra_i in cgras:
        atms_i = atom_keys(gra_i)
        if [x for x in atms_i if x in adj_atms[atm] or x == atm]:
            ret_gras.append(gra_i)
    return ret_gras


def branches(gra, atm_key, keep_root=False, stereo=False):
    """Get the branches extending from an atom, as a list

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_key: the root atom from which the branch originates
    :type atm_key: int
    :param keep_root: Keep root atom as part of branch? defaults to False
    :type keep_root: bool, optional
    :param stereo: Keep stereochemistry in the subgraph? defaults to False
    :type stereo: bool, optional
    :return: A dictionary associating each neighbor of `atm_key` with the
        subgraph for its branch
    :rtype: dict(int: automol graph data structure)
    """
    bnch_dct = branch_dict(gra, atm_key, keep_root=keep_root, stereo=stereo)
    return list(bnch_dct.values())


def branch_dict(gra, atm_key, keep_root=False, stereo=False):
    """Get the branches extending from an atom, by neighboring key

    If `keep_root` is set to `True`, then `atm_key` will be included in each
    branch. When this is done, note that if `atm_key` is part of a ring, then
    the branch for each neighboring key will only include that specific
    neighbor's bond to `atm_key`, not that of any other neighbors. This ensures
    that each branch has a distinct subgraph, even when `atm_key` is part of
    one or more rings. The latter will not be the case if `keep_root` is set to
    False.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_key: the root atom from which the branch originates
    :type atm_key: int
    :param keep_root: Keep root atom as part of branch? defaults to False
    :type keep_root: bool, optional
    :param stereo: Keep stereochemistry in the subgraph? defaults to False
    :type stereo: bool, optional
    :return: A dictionary associating each neighbor of `atm_key` with the
        subgraph for its branch
    :rtype: dict(int: automol graph data structure)
    """
    island_gra = remove_bonds(gra, atom_bond_keys(gra, atm_key), stereo=stereo)
    comps = connected_components(island_gra, stereo=stereo)
    root_atm = subgraph_(gra, {atm_key}, stereo=stereo)

    bnch_dct = {}
    nkeys = atom_neighbor_atom_keys(gra, atm_key)
    for nkey in nkeys:
        bnch = next(c for c in comps if nkey in atom_keys(c))
        if keep_root:
            bnch = union(root_atm, bnch)
            bnch = add_bonds(bnch, [{atm_key, nkey}])
        bnch_dct[nkey] = bnch
    return bnch_dct


def branch(
    gra,
    atm_key,
    branch_key,
    keep_root: bool = False,
    stereo: bool = False,
    return_missing_neighbor: bool = False,
):
    """Get the atom keys for a branch off of an atom extending a long a
    specific neighboring atom or bond.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_key: the root atom from which the branch originates
    :type atm_key: int
    :param branch_key: the adjacent atom or bond key that begins the branch
    :type branch_key: int or frozenset[int]
    :param keep_root: Keep root atom as part of branch? defaults to False
    :param stereo: Keep stereochemistry in the subgraph? defaults to False
    :param return_missing_neighbor: If the branch key is not a neighbor, return it as a
        one-atom branch?
    :return: The branch atom keys
    :rtype: frozenset[int]
    """
    # If `branch_key` was given as a bond key, get the atom
    if not isinstance(branch_key, numbers.Number):
        assert len(branch_key) == 2, f"Invalid branch key: {branch_key}"
        (branch_key,) = set(branch_key) - {atm_key}

    assert branch_key in atoms_neighbor_atom_keys(gra, atm_key), (
        f"No branch from {atm_key} along non-adjacent atom {branch_key}.\n"
        f"Input graph:\n{string(gra)}"
    )
    bnch_dct = branch_dict(gra, atm_key, keep_root=keep_root, stereo=stereo)
    if return_missing_neighbor:
        return bnch_dct.get(branch_key, subgraph_(gra, {branch_key}, stereo=stereo))

    return bnch_dct[branch_key]


def branch_atom_keys(
    gra,
    atm_key,
    branch_key,
    keep_root: bool = False,
    return_missing_neighbor: bool = False,
):
    """Get the atom keys for a branch off of an atom extending a long a
    specific neighboring atom or bond

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_key: the root atom from which the branch originates
    :type atm_key: int
    :param branch_key: the adjacent atom or bond key that begins the branch
    :type branch_key: int or frozenset[int]
    :param keep_root: Keep root atom as part of branch? defaults to False
    :param return_missing_neighbor: If the branch key is not a neighbor, return it as a
        one-atom branch?
    :return: The branch atom keys
    :rtype: frozenset[int]
    """
    return atom_keys(
        branch(
            gra,
            atm_key,
            branch_key,
            keep_root=keep_root,
            return_missing_neighbor=return_missing_neighbor,
        )
    )


def is_branched(gra):
    """determine is the molecule has a branched chain"""
    _is_branched = False
    gra = implicit(gra)
    chain_length = len(longest_chain(gra))
    natoms = atom_count(gra, with_implicit=False)

    if natoms != chain_length:
        _is_branched = True
    return _is_branched


# # rings
def rings(gra, ts_=True):
    """rings in the graph (minimal basis)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such
    :type ts_: bool
    :returns: A set of automol graph data structures for each ring
    """
    rng_bnd_keys_lst = rings_bond_keys(gra, ts_=ts_)
    gras = [bond_induced_subgraph(gra, bnd_keys) for bnd_keys in rng_bnd_keys_lst]
    return tuple(sorted(gras, key=frozen))


def rings_atom_keys(gra, ts_=True):
    """atom keys for each ring in the graph sorted by connectivity (minimal
    basis)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such
    :type ts_: bool
    :returns: A set of tuples of atom keys for each ring.
    """
    gra = gra if ts_ else ts_reactants_graph_without_stereo(gra)

    # Define weights, giving the breaking and forming bonds a higher value so that they
    # are only included if necessary
    ord_dct = bond_orders(gra)
    wgt_dct = util.dict_.transform_values(
        ord_dct, lambda o: int(o * 10000 if o < 1 else o)
    )

    nxg = _01networkx.from_graph(gra, edge_attrib_dct={"weight": wgt_dct})
    return frozenset(_01networkx.minimum_cycle_basis(nxg, weight="weight"))


def rings_bond_keys(gra, ts_=True):
    """bond keys for each ring in the graph (minimal basis)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such
    :type ts_: bool
    :returns: A set of sets of bond keys for each ring.
    """
    rng_atm_keys_lst = rings_atom_keys(gra, ts_=ts_)
    return frozenset(
        frozenset(map(frozenset, util.ring.edges(ks))) for ks in rng_atm_keys_lst
    )


def sorted_ring_atom_keys(rng):
    """get a ring's atom keys, sorted in order of connectivity"""
    return sorted_ring_atom_keys_from_bond_keys(bond_keys(rng))


def sorted_ring_atom_keys_from_bond_keys(rng_bnd_keys):
    """get a ring's atom keys, sorted in order of connectivity, from its bond
    keys
    """
    rng_bnd_keys = list(rng_bnd_keys)
    bnd_key = min(rng_bnd_keys, key=sorted)
    first_atm_key, atm_key = sorted(bnd_key)
    rng_bnd_keys.remove(bnd_key)
    rng_atm_keys = [first_atm_key, atm_key]
    while rng_bnd_keys:
        bnd_key = next(filter(lambda x: atm_key in x, rng_bnd_keys))
        rng_bnd_keys.remove(bnd_key)
        bnd_key = set(bnd_key)
        bnd_key.remove(atm_key)
        atm_key = next(iter(bnd_key))
        rng_atm_keys.append(atm_key)
    rng_atm_keys.pop(-1)
    rng_atm_keys = tuple(rng_atm_keys)
    return rng_atm_keys


def is_ring_key_sequence(gra, keys):
    """does this sequence of keys share a ring?"""
    keys = set(keys)
    return any(keys <= rng_keys for rng_keys in rings_atom_keys(gra))


def ring_arc_complement_atom_keys(gra, rng):
    """non-intersecting arcs from a ring that shares segments with a graph"""
    gra_atm_bnd_dct = atoms_bond_keys(gra)
    rng_atm_bnd_dct = atoms_bond_keys(rng)

    # 1. find divergence points, given by the atom at which the divergence
    # occurs and the bond followed by the ring as it diverges
    div_dct = {}

    for atm_key in atom_keys(gra) & atom_keys(rng):
        div = rng_atm_bnd_dct[atm_key] - gra_atm_bnd_dct[atm_key]
        if div:
            (bnd_key,) = div
            div_dct[atm_key] = bnd_key

    # 2. cycle through the ring atoms; if you meet a starting divergence, start
    # an arc; extend the arc until you meet an ending divergence; repeat until
    # all divergences are accounted for
    atm_keys = sorted_ring_atom_keys_from_bond_keys(bond_keys(rng))

    arcs = []
    arc = []
    for atm_key, next_atm_key in mit.windowed(itertools.cycle(atm_keys), 2):
        bnd_key = frozenset({atm_key, next_atm_key})

        # if we haven't started an arc, see if we are at a starting divergence;
        # if so, start the arc now and cross the divergence from our list
        if not arc:
            if atm_key in div_dct and div_dct[atm_key] == bnd_key:
                div_dct.pop(atm_key)

                arc.append(atm_key)
        # if we've started an arc, extend it; then, check if we are at an
        # ending divergence; if so, end the arc and cross the divergence from
        # our list; add it to our list of arcs
        else:
            arc.append(atm_key)

            if next_atm_key in div_dct and div_dct[next_atm_key] == bnd_key:
                div_dct.pop(next_atm_key)

                arc.append(next_atm_key)
                arcs.append(arc)
                arc = []

        # if no divergences are left, break out of the loop
        if not div_dct:
            break

    arcs = tuple(map(tuple, arcs))
    return arcs


def ring_systems(gra, lump_spiro=True):
    """polycyclic ring systems in the graph"""
    gras = [
        bond_induced_subgraph(gra, bnd_keys, stereo=True)
        for bnd_keys in ring_systems_bond_keys(gra, lump_spiro=lump_spiro)
    ]
    return tuple(sorted(gras, key=frozen))


def ring_systems_atom_keys(gra, lump_spiro=True):
    """bond keys for polycyclic ring systems in the graph"""
    atm_keys_lst = tuple(map(atom_keys, ring_systems(gra, lump_spiro=lump_spiro)))
    return atm_keys_lst


def ring_systems_bond_keys(gra, lump_spiro=True):
    """bond keys for polycyclic ring systems in the graph

    :param lump_spiro: Lump spiro rings into the same system?, default True
    :type lump_spiro: bool, optional
    """

    def _are_connected(bnd_keys1, bnd_keys2):
        """see if two rings are connected based on their bond keys"""
        atm_keys1 = functools.reduce(operator.or_, bnd_keys1)
        atm_keys2 = functools.reduce(operator.or_, bnd_keys2)
        ret = set(bnd_keys1) & set(bnd_keys2)
        if lump_spiro:
            ret |= set(atm_keys1) & set(atm_keys2)
        return ret

    rng_bnd_keys_lst = rings_bond_keys(gra)
    rsy_bnd_keys_lsts = util.equivalence_partition(rng_bnd_keys_lst, _are_connected)
    rsy_bnd_keys_lst = [
        frozenset(functools.reduce(operator.or_, bnd_keys_lst))
        for bnd_keys_lst in rsy_bnd_keys_lsts
    ]
    return rsy_bnd_keys_lst


def spiros(gra: Any, repeat_end: bool = False) -> Dict[int, Tuple[AtomKeys, AtomKeys]]:
    """Identify spiro-joined rings in the graph

    :param gra: Molecular graph
    :param repeat_end: Repeat the spiro atom at the end of the ring keys?
    :returns: A dictionary of pairs of spiro-joined rings, by joining spiro atom
    """
    rkeys_lst = rings_atom_keys(gra)
    spi_dct = {}
    for rkeys1, rkeys2 in itertools.combinations(rkeys_lst, r=2):
        shared_keys = set(rkeys1) & set(rkeys2)
        if len(shared_keys) == 1:
            (spi_key,) = shared_keys
            rkeys1 = list(util.ring.cycle_item_to_front(rkeys1, spi_key))
            rkeys2 = list(util.ring.cycle_item_to_front(rkeys2, spi_key))
            if repeat_end:
                rkeys1.append(spi_key)
                rkeys2.append(spi_key)
            spi_dct[spi_key] = (tuple(rkeys1), tuple(rkeys2))
    return spi_dct


def spiro_atom_keys(gra):
    """Identify spiro atoms in the graph

    Spiro atoms are the only atoms shared by one or more ring systems

    :param gra: Molecular graph
    :type gra: automol graph data structure
    """
    return frozenset(spiros(gra))


def is_ring_system(gra):
    """is this graph a ring system?"""
    gra = without_stereo(gra)
    return union_from_sequence(rings(gra), check=False) == gra


def ring_system_decomposed_atom_keys(rsy, rng_keys=None, check=True):
    """decomposed atom keys for a polycyclic ring system in a graph

    The ring system is decomposed into a ring and a series of arcs that can
    be used to successively construct the system

    :param rsy: the ring system
    :param rng_keys: keys for the first ring in the decomposition; if None, the
        smallest ring in the system will be chosen
    """
    if rng_keys is None:
        rng = sorted(rings(rsy), key=atom_count)[0]
        rng_keys = sorted_ring_atom_keys(rng)

    # check the arguments, if requested
    if check:
        # check that the graph is connected
        assert is_connected(rsy), "Ring system can't be disconnected."

        # check that the graph is actually a ring system
        assert is_ring_system(rsy), f"This is not a ring system graph:\n{string(rsy):s}"

        # check that rng is a subgraph of rsy
        assert set(rng_keys) <= atom_keys(rsy), (
            f"{string(rsy, one_indexed=False)}\n^ "
            "Rings system doesn't contain ring as subgraph:\n"
            f"{str(rng_keys)}"
        )

    bnd_keys = list(mit.windowed(rng_keys + rng_keys[:1], 2))

    # Save some information to handle spiros
    spi_dct = spiros(rsy, repeat_end=True)

    # Remove bonds for the ring
    rsy = remove_bonds(rsy, bnd_keys)
    keys_lst = [rng_keys]
    done_keys = set(rng_keys)

    while bond_keys(rsy):
        sp_dct = atom_shortest_paths(rsy)
        arc_pool = [
            sp_dct[i][j]
            for i, j in itertools.combinations(done_keys, 2)
            if j in sp_dct[i]
        ]

        # Include spiro rings as arcs
        rng_keys_pool = list(map(set, rings_atom_keys(rsy)))
        arc_pool.extend(
            rk
            for k, rks in spi_dct.items()
            for rk in rks
            if k in done_keys and set(rk) in rng_keys_pool
        )

        # Identify the next shortest arc
        arc_keys = min(arc_pool, key=len)
        # # Determine shortest paths for the graph with one more ring/arc deleted
        # sp_dct = atom_shortest_paths(rsy)

        # # The shortest path will be the next shortest arc in the system
        # # fused
        # all_arcs = list(
        #     sp_dct[i][j] for i, j in itertools.combinations(done_keys, 2)
        #     if j in sp_dct[i])
        # # spiro
        # spiro_key = None
        # if not all_arcs:
        #     for bnd_key in bond_keys(rsy):
        #         keya, keyb = bnd_key
        #         if keya in done_keys:
        #             spiro_key = keya
        #         elif keyb in done_keys:
        #             spiro_key = keyb
        #             keyb = keya
        #         if spiro_key is not None:
        #             sp_dct = atom_shortest_paths(remove_bonds(rsy, (bnd_key,)))
        #             all_arcs = (sp_dct[spiro_key][keyb],)
        #             break
        # arc_keys = min(all_arcs, key=len)

        # Add this arc to the list
        keys_lst.append(arc_keys)

        # Add these keys to the list of done keys
        done_keys |= set(arc_keys)

        # Delete the bond keys for the new arc and continue to the next iteration
        bnd_keys = list(map(frozenset, mit.windowed(arc_keys, 2)))
        # if spiro_key is not None:
        #     bnd_keys += (bnd_key,)
        rsy = remove_bonds(rsy, bnd_keys)
    keys_lst = tuple(map(tuple, keys_lst))
    return keys_lst


def ring_systems_decomposed_atom_keys(gra):
    """decomposed atom keys for polycyclic ring systems in the graph

    each ring system is decomposed into a ring and a series of arcs that can be
    used to successively construct the system
    """
    rsys = ring_systems(gra)
    decomps = tuple(map(_decompose_ring_system_atom_keys, rsys))
    return decomps


def _decompose_ring_system_atom_keys(rsy):
    """decompose a ring system into a ring and a series of arcs"""
    # sort from smallest to largest
    rngs_pool = sorted(rings(rsy), key=lambda x: atom_count(x, with_implicit=False))

    decomp = ()
    decomp_bnd_keys = set({})

    rng = rngs_pool.pop(0)
    bnd_keys = bond_keys(rng)
    atm_keys = sorted_ring_atom_keys_from_bond_keys(bnd_keys)

    decomp += (atm_keys,)
    decomp_bnd_keys.update(bnd_keys)

    while rngs_pool:
        decomp_rsy = bond_induced_subgraph(rsy, decomp_bnd_keys)
        for idx, rng in enumerate(rngs_pool):
            arcs = ring_arc_complement_atom_keys(decomp_rsy, rng)
            if arcs:
                rngs_pool.pop(idx)
                decomp += arcs
                decomp_bnd_keys.update(bond_keys(rng))

    return decomp


def is_bimolecular(tsg) -> bool:
    """Is this a TS for a bimolecular reaction?

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: `True` if it is, `False` if it isn't
    :rtype: bool
    """
    return len(connected_components(ts_reactants_graph_without_stereo(tsg))) == 2


def reacting_rings_bond_keys(tsg):
    """get the bond keys to rings forming in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    bnd_keys = ts_reacting_bond_keys(tsg)
    rngs_bnd_keys = tuple(
        bks for bks in rings_bond_keys(tsg, ts_=True) if bnd_keys & bks
    )
    return rngs_bnd_keys


def reacting_rings_atom_keys(tsg) -> List[List[int]]:
    """Get the atom keys to rings containing breaking or forming bonds

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    rngs_bnd_keys = reacting_rings_bond_keys(tsg)
    rngs_atm_keys = tuple(map(sorted_ring_atom_keys_from_bond_keys, rngs_bnd_keys))
    return rngs_atm_keys


def forming_rings_bond_keys(tsg):
    """get the bond keys to rings forming in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    frm_bnd_keys = ts_forming_bond_keys(tsg)
    frm_rngs_bnd_keys = tuple(
        bks for bks in rings_bond_keys(tsg, ts_=True) if frm_bnd_keys & bks
    )
    return frm_rngs_bnd_keys


def forming_rings_atom_keys(tsg) -> List[List[int]]:
    """Get the atom keys to rings forming in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    frm_rngs_bnd_keys = forming_rings_bond_keys(tsg)
    frm_rngs_atm_keys = tuple(
        map(sorted_ring_atom_keys_from_bond_keys, frm_rngs_bnd_keys)
    )
    return frm_rngs_atm_keys


def breaking_rings_bond_keys(tsg):
    """get the bond keys to rings breaking in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    brk_bnd_keys = ts_breaking_bond_keys(tsg)
    brk_rngs_bnd_keys = tuple(
        bks for bks in rings_bond_keys(tsg, ts_=True) if brk_bnd_keys & bks
    )
    return brk_rngs_bnd_keys


def breaking_rings_atom_keys(tsg):
    """get the atom keys to rings breaking in the TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    """
    brk_rngs_bnd_keys = breaking_rings_bond_keys(tsg)
    brk_rngs_atm_keys = tuple(
        map(sorted_ring_atom_keys_from_bond_keys, brk_rngs_bnd_keys)
    )
    return brk_rngs_atm_keys


def has_reacting_ring(tsg) -> bool:
    """Does this TS graph have a ring which is involved in the reaction?

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return bool(forming_rings_atom_keys(tsg) or breaking_rings_atom_keys(tsg))
