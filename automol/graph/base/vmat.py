""" graph-based z-matrix builder
"""

import itertools
from typing import Any, List

from ... import vmat
from ...util import dict_
from ...util import ring as ring_
from ._00core import (
    BondKeys,
    atom_count,
    atom_keys,
    atom_symbols,
    atoms_neighbor_atom_keys,
    atoms_zmat_sorted_neighbor_atom_keys,
    bond_keys,
    bonds_neighbor_atom_keys,
    is_ts_graph,
    remove_bonds,
    sort_by_size,
    string,
    subgraph,
    terminal_atom_keys,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reactants_graph_without_stereo,
)
from ._02algo import (
    connected_components,
    has_reacting_ring,
    is_bimolecular,
    is_connected,
    reacting_rings_atom_keys,
    reacting_rings_bond_keys,
    ring_system_decomposed_atom_keys,
    ring_systems,
    rings,
    rings_atom_keys,
    shortest_path_between_groups,
    sorted_ring_atom_keys,
)
from ._04class import atom_transfers, ring_forming_scissions


def ts_zmatrix_sorted_reactants_keys(tsg) -> List[List[int]]:
    """For bimolecular reactions without a TS ring, return keys for the reactants in the
    order they should appear in the z-matrix

    For unimolecular reactions or bimolecular reactions with a TS ring, returns None

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: For bimolecular reactions without a TS ring, the keys for each reactant,
        in z-matrix order, i.e. atom donor first or, if there isn't one, larger reactant
        first, as measured by (heavy atoms, total atoms, electrons)
    :rtype: List[int]
    """
    if not is_ts_graph(tsg) or has_reacting_ring(tsg) or not is_bimolecular(tsg):
        return None

    rcts_gra = ts_reactants_graph_without_stereo(tsg, dummy=True)
    rct_gras = connected_components(rcts_gra)

    # 1. If there is an atom transfer, put the donor reagent first
    tra_keys = set(atom_transfers(tsg))
    if tra_keys:
        rct_gras = sorted(rct_gras, key=lambda g: atom_keys(g) & tra_keys, reverse=True)
    # 2. Otherwise, put the larger reagent first
    else:
        rct_gras = sort_by_size(rct_gras)

    rcts_keys = tuple(map(tuple, map(sorted, map(atom_keys, rct_gras))))
    return rcts_keys


def ts_zmatrix_starting_ring_keys(tsg) -> List[int]:
    """Return keys for a TS ring to start from, sorted in z-matrix order

    If there isn't a TS ring, this returns `None`

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: The ring keys, sorted to exclude breaking bonds and include forming bonds
        as late as possible
    :rtype: List[int]
    """
    if not has_reacting_ring(tsg):
        return None

    rngs_keys = reacting_rings_atom_keys(tsg)
    rngs_bkeys = reacting_rings_bond_keys(tsg)

    if len(rngs_keys) > 1:
        raise NotImplementedError(f"Not implemented for multiple reacting rings: {tsg}")

    (rng_keys,) = rngs_keys
    (rng_bkeys,) = rngs_bkeys
    frm_bkeys = [bk for bk in ts_forming_bond_keys(tsg) if bk in rng_bkeys]
    brk_bkeys = [bk for bk in ts_breaking_bond_keys(tsg) if bk in rng_bkeys]

    # Determine the neighbor atom counts for each bond
    bnkeys_dct = bonds_neighbor_atom_keys(tsg, group=False)
    bncount_dct = dict_.transform_values(bnkeys_dct, len)
    rem_bkeys = sorted(
        set(rng_bkeys) - (set(frm_bkeys) | set(brk_bkeys)), key=bncount_dct.get
    )

    # If there are forming keys in the ring, drop the breaking keys, if possible
    if frm_bkeys:
        drop_bkeys = brk_bkeys if brk_bkeys else rem_bkeys
        back_keys = list(itertools.chain(*frm_bkeys))
    else:
        drop_bkeys = rem_bkeys
        back_keys = list(itertools.chain(*brk_bkeys))

    # Special handling for ring-forming scissions
    rsciss_dct = ring_forming_scissions(tsg)
    if rsciss_dct:
        tra_key, (_, frm_nkey) = next(iter(rsciss_dct.items()))
        drop_bkeys = [frozenset({tra_key, frm_nkey})]
        back_keys = [tra_key]

    return ring_.cycle_to_optimal_split(rng_keys, drop_bkeys, back_keys)


def vmatrix(gra):
    """V-matrix for a connected graph or TS graph

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :returns: The v-matrix
    :rtype: automol vmat data structure
    """
    # Generate a v-matrix for the graph and get the z-matrix reordering
    rng_keys = ts_zmatrix_starting_ring_keys(gra)
    rcts_keys = ts_zmatrix_sorted_reactants_keys(gra)
    if rcts_keys is not None:
        rct1_keys, rct2_keys = rcts_keys
        vma, zma_keys = connected_vmatrix(gra, keys=rct1_keys)
        vma, zma_keys = continue_vmatrix(gra, rct2_keys, vma, zma_keys)
    else:
        vma, zma_keys = connected_vmatrix(gra, rng_keys=rng_keys)

    return vma, tuple(zma_keys)


def connected_vmatrix(gra, keys=None, rng_keys=None):
    """v-matrix for a connected graph

    :param gra: the graph
    :param keys: restrict the v-matrix to a subset of keys, which must span a
        connected graph
    :param rng_keys: keys for a ring to start from
    """
    if keys is not None:
        gra = subgraph(gra, keys)

    assert is_connected(gra), "Graph must be connected!"

    # Start with the ring systems and their connections. If there aren't any,
    # start with the first terminal atom
    if ring_systems(gra):
        vma, zma_keys = connected_ring_systems(gra, rng_keys=rng_keys)
    else:
        term_keys = sorted(terminal_atom_keys(gra, backbone=True))
        if term_keys:
            start_key = term_keys[0]
        else:
            start_key = sorted(atom_keys(gra))[0]

        vma, zma_keys = start_at(gra, start_key)

    rem_keys = atom_keys(gra) - set(zma_keys)
    vma, zma_keys = continue_vmatrix(gra, rem_keys, vma, zma_keys)
    return vma, zma_keys


def missing_bond_keys(gra: Any, vma: Any) -> BondKeys:
    """Identify the bonds that are missing in a V-Matrix

    :param gra: A graph
    :param vma: A V-Matrix
    :returns: The missing bond keys
    """
    dist_names = vmat.distance_names(vma)
    dist_coos = dict_.values_by_key(vmat.coordinates(vma, multi=False), dist_names)
    miss_bkeys = bond_keys(gra) - set(map(frozenset, dist_coos))
    return miss_bkeys


def continue_vmatrix(gra, keys, vma, zma_keys):
    """continue a v-matrix for a subset of keys, starting from a partial
    v-matrix
    """
    gra = subgraph(gra, set(keys) | set(zma_keys))

    vma, zma_keys = continue_connected_ring_systems(gra, keys, vma, zma_keys)

    # Complete any incomplete branches
    branch_keys = _atoms_missing_neighbors(gra, zma_keys)
    for key in branch_keys:
        vma, zma_keys = complete_branch(gra, key, vma, zma_keys)

    return vma, zma_keys


def connected_ring_systems(gra, rng_keys=None, check=True):
    """generate a v-matrix covering a graph's ring systems and the connections
    between them
    """
    if check:
        assert is_connected(gra), "Graph must be connected!"

    rsys = sorted(ring_systems(gra), key=atom_count)
    # Construct the v-matrix for the first ring system, choosing which ring
    # to start from
    if rng_keys is None:
        rsy = rsys.pop(0)
        rngs = sorted(rings(rsy), key=atom_count)
        rng_keys = sorted_ring_atom_keys(rngs.pop(0))
    else:
        idx = next(
            (i for i, ks in enumerate(map(atom_keys, rsys)) if set(rng_keys) <= ks),
            None,
        )
        assert (
            idx is not None
        ), f"The ring {str(rng_keys)} is not in this graph:\n{string(gra)}"
        rsy = rsys.pop(idx)

    keys_lst = list(ring_system_decomposed_atom_keys(rsy, rng_keys=rng_keys))

    vma, zma_keys = ring_system(gra, keys_lst)

    keys = atom_keys(gra) - set(zma_keys)

    vma, zma_keys = continue_connected_ring_systems(
        gra, keys, vma, zma_keys, rsys=rsys, check=False
    )

    return vma, zma_keys


def continue_connected_ring_systems(gra, keys, vma, zma_keys, rsys=None, check=True):
    """generate the connected ring systems for a subset of keys, continuing on
    from a partial v-matrix

    The subset must have at least one neighbor that already exists in the
    v-matrix

    :param gra: the graph for which the v-matrix will be constructed
    :param keys: the subset of keys to be added to the v-matrix
    :param vma: a partial v-matrix from which to continue
    :param zma_keys: row keys for the partial v-matrix, identifying the atom
        specified by each row of `vma` in order
    :param rsys: optionally, pass the ring systems in to avoid recalculating
    """
    gra = subgraph(gra, set(keys) | set(zma_keys))
    sub = subgraph(gra, keys)
    if check:
        assert is_connected(gra), "Graph must be connected!"

    if rsys is None:
        rsys = sorted(ring_systems(sub), key=atom_count)

    rsys = list(rsys)

    while rsys:
        # Find the next ring system with a connection to the current
        # v-vmatrix and connect them
        conn = False
        for idx, rsy_keys in enumerate(map(atom_keys, rsys)):
            if set(zma_keys) & rsy_keys:
                # ring systems are connected by one bond -- no chain needed
                keys = set(zma_keys) & rsy_keys
                assert len(keys) == 1, (
                    "Attempting to add redundant keys to v-matrix: " f"{str(keys)}"
                )
                (key,) = keys

                conn = True
            else:
                # see if the ring systems are connected by a chain
                keys = shortest_path_between_groups(gra, zma_keys, rsy_keys)

                # if so, build a bridge from the current v-matrix to this next
                # ring system
                vma, zma_keys = continue_chain(
                    gra, keys[:-1], vma, zma_keys, term_hydrogens=False
                )
                key = keys[-1]

                conn = bool(keys is not None)

            if conn:
                rsy = rsys.pop(idx)
                break

        assert keys is not None, "This is a disconnected graph!"

        # 2. Decompose the ring system with the connecting ring first
        rng_keys = next(rks for rks in rings_atom_keys(rsy) if key in rks)
        keys_lst = ring_system_decomposed_atom_keys(rsy, rng_keys=rng_keys)

        # 3. Build the next ring system
        vma, zma_keys = continue_ring_system(gra, keys_lst, vma, zma_keys)

    return vma, zma_keys


def ring_system(gra, keys_lst):
    """generate a v-matrix for a ring system

    :param gra: the graph
    :param keys_lst: the first entry contains keys for a ring and each next one
        contains keys for an arc that starts and ends on atoms in the preceding
        entries
    """
    # First, get the ring keys
    keys_lst = list(keys_lst)
    keys = keys_lst.pop(0)
    incl_keys = keys
    # Break the bonds joining the last pair of atoms in each arc
    gra = remove_bonds(gra, [(k[-1], k[-2]) for k in keys_lst])

    # Start by constructing the v-matrix for the first ring
    vma, zma_keys = ring(gra, keys)
    # Now, complete the ring system by continuing the v-matrix along each arc
    for keys in keys_lst:
        start_key = 1
        if keys[-1] in incl_keys:
            end_key = -1
        else:
            end_key = None
        # Note that the atoms on each end of the arc are already in the
        # v-matrix, so we ignore those
        vma, zma_keys = continue_chain(gra, keys[start_key:end_key], vma, zma_keys)
        incl_keys += keys[start_key:end_key]
    return vma, zma_keys


def continue_ring_system(gra, keys_lst, vma, zma_keys):
    """continue constructing a v-matrix for a ring system

    Exactly one atom in the ring system must already be in the v-matrix, and
    this atom must be in the starting ring of the decomposed ring system key
    list.
    """
    # First, get the ring keys
    keys_lst = list(keys_lst)
    keys = keys_lst.pop(0)

    # Break the bonds joining the last pair of atoms in each arc
    gra = remove_bonds(gra, [(k[-1], k[-2]) for k in keys_lst])

    # Start by constructing the v-matrix for the first ring
    vma, zma_keys = continue_ring(gra, keys, vma, zma_keys)

    # Now, complete the ring system by continuing the v-matrix along each arc
    for keys in keys_lst:
        # Note that the atoms on each end of the arc are already in the
        # v-matrix, so we ignore those
        vma, zma_keys = continue_chain(
            gra, keys[1:-1], vma, zma_keys, term_hydrogens=False
        )

    return vma, zma_keys


def ring(gra, keys):
    """generate a v-matrix for a ring

    All neighboring atoms along the ring will be included

    :param gra: the graph
    :param keys: ring keys, in the order they should appear in the z-matrix
    """
    # Break the bond between the first and last atoms to make this a chain
    gra = remove_bonds(gra, [(keys[0], keys[-1])])

    # Now, construct a v-matrix for the chain
    vma, zma_keys = chain(gra, keys, term_hydrogens=True)
    return vma, zma_keys


def continue_ring(gra, keys, vma, zma_keys):
    """continue constructing a v-matrix around a ring

    All neighboring atoms along the ring will be included

    Exactly one atom in the ring must already be in the v-matrix.
    """
    # Find the connecting key
    key = next((k for k in keys if k in zma_keys), None)
    assert key is not None, "There must be a ring atom already in the v-matrix"

    # Cycle the connecting key to the front of the ring
    keys = ring_.cycle_item_to_front(keys, key)

    # Break the bond between the first and last atoms to make this a chain
    gra = remove_bonds(gra, [(keys[0], keys[-1])])

    # Now, construct a v-matrix for the chain
    vma, zma_keys = continue_chain(gra, keys, vma, zma_keys)
    return vma, zma_keys


def chain(gra, keys, term_hydrogens=True):
    """generate a v-matrix for a chain

    All neighboring atoms along the chain will be included

    :param gra: the graph
    :param keys: a list of keys for the chain
    :param term_hydrogens: whether or not to extend the chain to include
        terminal hydrogens, if present
    """
    if term_hydrogens:
        keys = _extend_chain_to_include_terminal_hydrogens(gra, keys)

    # 1. Start the chain
    vma, zma_keys = start_at(gra, keys[0])

    # If the chain isn't complete, continue it
    if set(keys) - set(zma_keys):
        miss_keys = set(_atoms_missing_neighbors(gra, zma_keys))
        (start_key,) = set(keys) & miss_keys
        keys = keys[keys.index(start_key) :]

        # 2. Continue the chain
        if keys:
            vma, zma_keys = continue_chain(gra, keys, vma=vma, zma_keys=zma_keys)

    return vma, zma_keys


def continue_chain(gra, keys, vma, zma_keys, term_hydrogens=True):
    """continue constructing a v-matrix along a chain

    All neighboring atoms along the chain will be included

    Exactly one atom in the chain must already be in the v-matrix

    :param gra: the graph for which the v-matrix will be constructed
    :param keys: the keys for atoms along the chain, which must be contiguous;
        the first atom must already appear in the v-matrix
    :param vma: a partial v-matrix from which to continue
    :param zma_keys: row keys for the partial v-matrix, identifying the atom
        specified by each row of `vma` in order
    :param extend: whether to extend the chain's start to include the three
        anchoring atoms
    :param term_hydrogens: whether to extend the chain's end to include
        terminal hydrogens
    """
    if term_hydrogens:
        keys = _extend_chain_to_include_terminal_hydrogens(gra, keys, start=False)

    vma, zma_keys = complete_branch(gra, keys[0], vma, zma_keys, branch_keys=keys)

    return vma, zma_keys


def start_at(gra, key):
    """start a v-matrix at a specific atom

    Returns the started vmatrix, along with keys to atoms whose neighbors are
    missing from it
    """
    symb_dct = atom_symbols(gra)
    nkeys_dct = atoms_zmat_sorted_neighbor_atom_keys(gra)

    nkeys = nkeys_dct[key]
    if not nkeys:
        # There are no neighbors, so this will be monatomic
        zma_keys = [key]
    elif len(nkeys) == 1:
        # Need special handling for atoms with only one neighbor
        if symb_dct[key] == "X" or (
            symb_dct[key] == "H" and symb_dct[nkeys[0]] not in ("H", "X")
        ):
            key2 = nkeys[0]
            zma_keys = (key2,) + nkeys_dct[key2]
        else:
            key2 = nkeys[0]
            nkeys = tuple(k for k in nkeys_dct[key2] if k != key)
            zma_keys = (key, key2) + nkeys
    else:
        zma_keys = (key,) + nkeys_dct[key]

    vma = ()
    for row, key_ in enumerate(zma_keys):
        idx1 = idx2 = idx3 = None
        if row > 0:
            key1 = next(k for k in nkeys_dct[key_] if k in zma_keys[:row])
            idx1 = zma_keys.index(key1)
        if row > 1:
            key2 = next(k for k in nkeys_dct[key1] if k in zma_keys[:row] and k != key_)
            idx2 = zma_keys.index(key2)
        if row > 2:
            key3 = next(k for k in zma_keys[:row] if k not in (key_, key1, key2))
            idx3 = zma_keys.index(key3)

        sym = symb_dct[key_]
        key_row = [idx1, idx2, idx3]
        vma = vmat.add_atom(vma, sym, key_row)

    return vma, zma_keys


def complete_branch(gra, key, vma, zma_keys, branch_keys=None):
    """continue constructing a v-matrix along a chain

    All neighboring atoms along the chain will be included

    Exactly one atom in the chain must already be in the v-matrix

    :param gra: the graph for which the v-matrix will be constructed
    :param keys: the keys for atoms along the chain, which must be contiguous;
        the first atom must already appear in the v-matrix
    :param vma: a partial v-matrix from which to continue
    :param zma_keys: row keys for the partial v-matrix, identifying the atom
        specified by each row of `vma` in order
    :param branch_keys: optionally, restrict the v-matrix to these keys and
        their neighbors; if `None`, the entire branch will be included
    """
    branch_keys = atom_keys(gra) if branch_keys is None else branch_keys
    keys = _extend_chain_to_include_anchoring_atoms(gra, [key], zma_keys)

    zma_keys = list(zma_keys)
    symb_dct = atom_symbols(gra)
    nkeys_dct = atoms_zmat_sorted_neighbor_atom_keys(gra, prioritize_keys=branch_keys)

    # If this z-matrix is being continued from a partial z-matrix, the leading
    # atom for a torsion may have already be defined. To handle this case, I
    # make a dictionary of these leading atoms and use them below where needed.
    lead_key_dct = {}
    for idx, key_row in enumerate(vmat.key_matrix(vma)):
        axis = key_row[:2]
        if None not in axis:
            axis = tuple(map(zma_keys.__getitem__, axis))
            if axis not in lead_key_dct:
                lead_key_dct[axis] = zma_keys[idx]

    def _continue(key1, key2, key3, vma, zma_keys):
        k3ns = list(nkeys_dct[key3])
        for k3n in set(k3ns) & set(zma_keys):
            k3ns.remove(k3n)

        if k3ns:
            key4 = k3ns.pop(0)

            lead_key = None
            if (key3, key2) in lead_key_dct:
                lead_key = lead_key_dct[(key3, key2)]

            # Add the leading atom to the v-matrix
            symb = symb_dct[key4]
            dkey = key1 if lead_key is None else lead_key
            key_row = [
                zma_keys.index(k) if k is not None else None for k in (key3, key2, dkey)
            ]
            vma = vmat.add_atom(vma, symb, key_row)
            assert key4 not in zma_keys, f"Atom {key4:d} already in v-matrix."
            zma_keys.append(key4)

            dkey = key4 if lead_key is None else lead_key
            # Add the neighbors of atom 3 (if any) to the v-matrix, decoupled
            # from atom 1 for properly decopuled torsions
            for k3n in k3ns:
                sym = symb_dct[k3n]

                if symb_dct[dkey] == "X":
                    key_row = list(map(zma_keys.index, (key3, dkey, key2)))
                else:
                    key_row = list(map(zma_keys.index, (key3, key2, dkey)))

                vma = vmat.add_atom(vma, sym, key_row)
                assert k3n not in zma_keys, f"Atom {k3n:d} already in v-matrix."
                zma_keys.append(k3n)

            # Recursion
            if key4 in branch_keys:
                vma, zma_keys = _continue(key2, key3, key4, vma, zma_keys)

            if symb_dct[key4] == "X":
                key2 = key4

            for k3n in k3ns:
                if k3n in branch_keys:
                    vma, zma_keys = _continue(key2, key3, k3n, vma, zma_keys)

        return vma, zma_keys

    key1, key2, key3 = keys[0], keys[1], keys[2]
    vma, zma_keys = _continue(key1, key2, key3, vma, zma_keys)

    return vma, zma_keys


# helpers
def _extend_chain_to_include_anchoring_atoms(gra, keys, zma_keys):
    """extend chain to include three atoms already specified in v-matrix

    :param gra: the graph
    :param keys: keys in the chain; the first atom should already be specified
    :param zma_keys: keys currently in the v-matrix
    """
    sub_gra = subgraph(gra, zma_keys)
    nkeys_dct = atoms_zmat_sorted_neighbor_atom_keys(sub_gra)

    symb_dct = atom_symbols(gra)

    key3 = keys[0]
    assert key3 in zma_keys

    # 1. Try to find a chain key3-key2-key1
    def _find_key1_and_key2():
        for key2 in nkeys_dct[key3]:
            key1s = nkeys_dct[key2] if symb_dct[key2] != "X" else nkeys_dct[key3][1:]
            for key1 in key1s:
                if key1 != key3:
                    return key1, key2
        return None, None

    key1, key2 = _find_key1_and_key2()
    # 2. If no complete chain can be found, allow key1 to be discontinuous
    if key1 is None and key2 is None:
        # Unless the starting z-matrix is monatomic, there *must* be a neighbor for key3
        key2 = next((k for k in nkeys_dct[key3] if k in zma_keys), None)
        # Since _find_key1_and_key2() failed, there must not be a neighbor for key2
        # aside from key3, so use another, discontinous z-matrix key, if there is one
        key1 = next((k for k in zma_keys if k not in (key2, key3)), None)

    keys = (
        key1,
        key2,
    ) + tuple(keys)

    return keys


def _extend_chain_to_include_terminal_hydrogens(gra, keys, start=True, end=True):
    """extend each end of a chain to include terminal hydrogens, if any"""
    symb_dct = atom_symbols(gra)
    nkeys_dct = atoms_neighbor_atom_keys(gra)

    sta_nkeys = nkeys_dct[keys[0]] - set(keys)
    end_nkeys = nkeys_dct[keys[-1]] - set(keys)

    sta_nkey = min((k for k in sta_nkeys if symb_dct[k] == "H"), default=None)
    end_nkey = min((k for k in end_nkeys if symb_dct[k] == "H"), default=None)

    keys = tuple(keys)

    if start and sta_nkey is not None:
        keys = (sta_nkey,) + keys

    if end and end_nkey is not None:
        keys = keys + (end_nkey,)

    return keys


def _atoms_missing_neighbors(gra, zma_keys):
    """get atoms from the list currently in the v-matrix with neighbors that
    are not in the v-matrix
    """
    nkeys_dct = atoms_neighbor_atom_keys(gra)
    keys = []
    for key in zma_keys:
        if any(k not in zma_keys for k in nkeys_dct[key]):
            keys.append(key)
    keys = tuple(keys)
    return keys
