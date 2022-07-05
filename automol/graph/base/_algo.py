""" algorithm functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import operator
import collections.abc
import functools
import itertools
import more_itertools as mit
from automol import util
from automol.util import dict_
from automol.graph.base import _networkx
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bond_keys
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import set_atom_symbols
from automol.graph.base._core import string
from automol.graph.base._core import frozen
from automol.graph.base._core import atom_count
from automol.graph.base._core import atom_symbol_keys
from automol.graph.base._core import remove_atoms
from automol.graph.base._core import remove_bonds
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import implicit
from automol.graph.base._core import explicit
from automol.graph.base._core import union_from_sequence
from automol.graph.base._core import subgraph
from automol.graph.base._core import bond_induced_subgraph
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atoms_bond_keys
from automol.graph.base._core import bonds_neighbor_bond_keys


# # isomorphisms and equivalence
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


def _isomorphism(gra1, gra2):
    """
    """
    nxg1 = _networkx.from_graph(gra1)
    nxg2 = _networkx.from_graph(gra2)
    iso_dct = _networkx.isomorphism(nxg1, nxg2)
    return iso_dct


def sequence_isomorphism(gras1, gras2, backbone_only=False, stereo=True,
                         dummy=True):
    """ Obtain an isomorphism between two sequences of graphs

    :param backbone_only: Compare backbone atoms only?
    :type backbone_only: bool
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: First, a sequence of integers mapping each graph in `gras1` to an
        isomorph in `gras2`, followes by a sequence of dictionaries mapping
        their atoms onto each other.
    :rtype: tuple[int], tuple[dict]
    """
    gras1_pool = dict(enumerate(gras1))
    order = []
    iso_dcts = []
    for gra2 in gras2:
        found_match = False
        for idx, gra1 in gras1_pool.items():
            iso_dct = isomorphism(gra2, gra1, backbone_only=backbone_only,
                                  stereo=stereo, dummy=dummy)
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


def full_isomorphism(gra1, gra2):
    """ full graph isomorphism

    TODO: DEPRECATE
    """
    assert gra1 == explicit(gra1) and gra2 == explicit(gra2)
    nxg1 = _networkx.from_graph(gra1)
    nxg2 = _networkx.from_graph(gra2)
    iso_dct = _networkx.isomorphism(nxg1, nxg2)
    return iso_dct


def full_subgraph_isomorphism(gra1, gra2):
    """ gra2 is fully isomorphic to a subgraph of gra1

    TODO: DEPRECATE
    """
    assert gra1 == explicit(gra1) and gra2 == explicit(gra2)
    nxg1 = _networkx.from_graph(gra1)
    nxg2 = _networkx.from_graph(gra2)
    iso_dct = _networkx.subgraph_isomorphism(nxg1, nxg2)
    return iso_dct


def backbone_isomorphism(gra1, gra2):
    """ graph backbone isomorphism

    TODO: DEPRECATE

    for implicit graphs, this is the relabeling of `gra1` to produce `gra2`
    for other graphs, it gives the correspondences between backbone atoms
    """
    gra1 = implicit(gra1)
    gra2 = implicit(gra2)
    nxg1 = _networkx.from_graph(gra1)
    nxg2 = _networkx.from_graph(gra2)
    iso_dct = _networkx.isomorphism(nxg1, nxg2)
    return iso_dct


def backbone_isomorphic(gra1, gra2):
    """ are these molecular graphs backbone isomorphic?

    TODO: DEPRECATE
    """
    return backbone_isomorphism(gra1, gra2) is not None


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
        f"{atm_key} not in {atom_keys(gra)}")

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
    """ Identify isomorphically unique atoms, which do not transform into each
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
    """ Identify isomorphically unique bonds, which do not transform into each
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
        return are_equivalent_bonds(gra, bnd1_key, bnd2_key, stereo=stereo,
                                    dummy=dummy)

    eq_classes = util.equivalence_partition(bnd_keys, _equiv)
    class_reps = frozenset(next(iter(c)) for c in eq_classes)
    return class_reps


def chem_unique_atoms_of_type(gra, symb):
    """ For the given atom type, determine the idxs of all the
         chemically unique atoms.

    TODO: DEPRECATE this can be done with the more general
    atom_equivalence_class_reps function. This approach to determining
    equivalence seems risky, because deleting an atom also deletes its bonds.
    In particular, if this is a disconnected graph with multiple radical
    molecules, then this approach breaks down completely and will treat atoms
    as equivalent when they are not. The approach used above is robust in these
    cases.

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :param symb: atomic symbol to determine symbols for
        :type symb: str
        :rtype: tuple(int)
    """

    # Get the indices for the atom type
    symb_idx_dct = atom_symbol_keys(gra)
    atom_idxs = symb_idx_dct[symb]

    # Loop over each idx
    uni_idxs = tuple()
    uni_del_gras = []
    for idx in atom_idxs:

        # Remove the atom from the graph
        del_gra = remove_atoms(gra, [idx])

        # Test if the del_gra is isomorphic to any of the uni_del_gras
        new_uni = True
        for uni_del_gra in uni_del_gras:
            iso_dct = full_isomorphism(del_gra, uni_del_gra)
            if iso_dct:
                new_uni = False
                break

        # Add graph and idx to lst if del gra is unique
        if new_uni:
            uni_del_gras.append(del_gra)
            uni_idxs += (idx,)

    return uni_idxs


# # algorithms
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


def is_connected(gra):
    """ is this a connected graph
    """
    return len(connected_components(gra)) == 1


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
    assert not set(keys1) & set(keys2), (
        f"{str(keys1):s} overlaps with {str(keys2):s}"
    )

    sp_dct = atom_shortest_paths(gra)
    keys = None
    for key1 in keys1:
        for key2 in keys2:
            if key2 in sp_dct[key1]:
                if keys is None or len(keys) > len(sp_dct[key1][key2]):
                    keys = sp_dct[key1][key2]

    return keys


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


def longest_chain(gra):
    """ longest chain in the graph
    """
    atm_keys = atom_keys(gra)

    max_chain = max((atom_longest_chain(gra, atm_key) for atm_key in atm_keys),
                    key=len)
    return max_chain


def weighted_maximal_matching(gra, bnd_weight_dct=None):
    """ Calculates a weighted maximal matching of a graph

        That is, a set of edges that covers the graph as much as possible
    """
    edge_attrib_dct = (None if bnd_weight_dct is None else
                       {'weight': bnd_weight_dct})
    nxg = _networkx.from_graph(gra, edge_attrib_dct=edge_attrib_dct)
    bnd_keys = _networkx.weighted_maximal_matching(
        nxg, edge_attrib_name='weight')
    return bnd_keys


# # branches and groups
def ring_atom_chirality(gra, atm, ring_atms, stereo=False):
    """is this ring atom a chiral center?
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


def atom_groups(gra, atm, stereo=False):
    """ return a list of groups off of one atom

    TODO: MERGE WITH BRANCH FUNCTIONS OR MAKE NAMING CONSISTENT SOMEHOW
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


# # rings
def rings(gra):
    """ rings in the graph (minimal basis)
    """
    gras = [bond_induced_subgraph(gra, bnd_keys)
            for bnd_keys in rings_bond_keys(gra)]
    return tuple(sorted(gras, key=frozen))


def rings_atom_keys(gra):
    """ atom keys for each ring in the graph sorted by connectivity (minimal basis)
    """
    rng_atm_keys_lst = frozenset(
        map(sorted_ring_atom_keys_from_bond_keys, rings_bond_keys(gra)))
    return rng_atm_keys_lst


def rings_bond_keys(gra):
    """ bond keys for each ring in the graph (minimal basis)
    """
    bnd_keys = bond_keys(gra)

    def _ring_bond_keys(rng_atm_keys):
        return frozenset(filter(lambda x: x <= rng_atm_keys, bnd_keys))

    nxg = _networkx.from_graph(gra)
    rng_atm_keys_lst = _networkx.minimum_cycle_basis(nxg)
    rng_bnd_keys_lst = frozenset(map(_ring_bond_keys, rng_atm_keys_lst))
    return rng_bnd_keys_lst


def sorted_ring_atom_keys(rng):
    """ get a ring's atom keys, sorted in order of connectivity
    """
    return sorted_ring_atom_keys_from_bond_keys(bond_keys(rng))


def sorted_ring_atom_keys_from_bond_keys(rng_bnd_keys):
    """ get a ring's atom keys, sorted in order of connectivity, from its bond
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
    """ does this sequence of keys share a ring?
    """
    keys = set(keys)
    return any(keys <= rng_keys for rng_keys in rings_atom_keys(gra))


def cycle_ring_atom_key_to_front(keys, key, end_key=None):
    """ helper function to cycle ring atom keys until one is in front

    :param keys: ring keys
    :parm key: the key to cycle to the font
    :param end_key: optionally, ensure that another key is the last key in the
        ring; note that this is only possible if key and end_key are adjacent
    """
    assert key in keys, (f"{key:d} is not in {str(keys):s}")
    keys = tuple(itertools.islice(
        itertools.dropwhile(lambda x: x != key, itertools.cycle(keys)),
        len(keys)))

    if end_key is not None and keys[-1] != end_key:
        assert keys[1] == end_key, (
            f"end_key {key:d} is not adjacent to {end_key:d} in the ring")
        keys = list(reversed(keys))
        keys = cycle_ring_atom_key_to_front(keys, key)

    return keys


def ring_arc_complement_atom_keys(gra, rng):
    """ non-intersecting arcs from a ring that shares segments with a graph
    """
    gra_atm_bnd_dct = atoms_bond_keys(gra)
    rng_atm_bnd_dct = atoms_bond_keys(rng)

    # 1. find divergence points, given by the atom at which the divergence
    # occurs and the bond followed by the ring as it diverges
    div_dct = {}

    for atm_key in atom_keys(gra) & atom_keys(rng):
        div = rng_atm_bnd_dct[atm_key] - gra_atm_bnd_dct[atm_key]
        if div:
            bnd_key, = div
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


def ring_systems(gra):
    """ polycyclic ring systems in the graph
    """
    gras = [bond_induced_subgraph(gra, bnd_keys, stereo=True)
            for bnd_keys in ring_systems_bond_keys(gra)]
    return tuple(sorted(gras, key=frozen))


def ring_systems_atom_keys(gra):
    """ bond keys for polycyclic ring systems in the graph
    """
    atm_keys_lst = tuple(map(atom_keys, ring_systems(gra)))
    return atm_keys_lst


def ring_systems_bond_keys(gra):
    """ bond keys for polycyclic ring systems in the graph
    """

    def _are_connected(bnd_keys1, bnd_keys2):
        """ see if two rings are connected based on their bond keys
        """
        atm_keys1 = functools.reduce(operator.or_, bnd_keys1)
        atm_keys2 = functools.reduce(operator.or_, bnd_keys2)
        common_bonds = set(bnd_keys1) & set(bnd_keys2)
        common_atoms = set(atm_keys1) & set(atm_keys2)
        return bool(common_bonds) or bool(common_atoms)

    rng_bnd_keys_lst = rings_bond_keys(gra)
    rsy_bnd_keys_lsts = util.equivalence_partition(rng_bnd_keys_lst,
                                                   _are_connected)
    rsy_bnd_keys_lst = [frozenset(functools.reduce(operator.or_, bnd_keys_lst))
                        for bnd_keys_lst in rsy_bnd_keys_lsts]
    return rsy_bnd_keys_lst


def is_ring_system(gra):
    """ is this graph a ring system?
    """
    gra = without_stereo_parities(gra)
    return union_from_sequence(rings(gra), check=False) == gra


def ring_system_decomposed_atom_keys(rsy, rng_keys=None, check=True):
    """ decomposed atom keys for a polycyclic ring system in a graph

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
        assert is_ring_system(rsy), (
            f"This is not a ring system graph:\n{string(rsy):s}")

        # check that rng is a subgraph of rsy
        assert set(rng_keys) <= atom_keys(rsy), (
            f"{string(rsy, one_indexed=False)}\n^ "
            "Rings system doesn't contain ring as subgraph:\n"
            f"{str(rng_keys)}")

    bnd_keys = list(mit.windowed(rng_keys + rng_keys[:1], 2))

    # Remove bonds for the ring
    rsy = remove_bonds(rsy, bnd_keys)
    keys_lst = [rng_keys]
    done_keys = set(rng_keys)

    while bond_keys(rsy):

        # Determine shortest paths for the graph with one more ring/arc deleted
        sp_dct = atom_shortest_paths(rsy)

        # The shortest path will be the next shortest arc in the system
        arc_keys = min(
            (sp_dct[i][j] for i, j in itertools.combinations(done_keys, 2)
             if j in sp_dct[i]), key=len)

        # Add this arc to the list
        keys_lst.append(arc_keys)

        # Add these keys to the list of done keys
        done_keys |= set(arc_keys)

        # Delete tbond keys for the new arc and continue to the next iteration
        bnd_keys = list(map(frozenset, mit.windowed(arc_keys, 2)))
        rsy = remove_bonds(rsy, bnd_keys)

    keys_lst = tuple(map(tuple, keys_lst))
    return keys_lst


def ring_systems_decomposed_atom_keys(gra):
    """ decomposed atom keys for polycyclic ring systems in the graph

    each ring system is decomposed into a ring and a series of arcs that can be
    used to successively construct the system
    """
    rsys = ring_systems(gra)
    decomps = tuple(map(_decompose_ring_system_atom_keys, rsys))
    return decomps


def _decompose_ring_system_atom_keys(rsy):
    """ decompose a ring system into a ring and a series of arcs
    """
    # sort from smallest to largest
    rngs_pool = sorted(
        rings(rsy), key=lambda x: atom_count(x, with_implicit=False))

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
