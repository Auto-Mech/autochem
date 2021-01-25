""" graph-based z-matrix builder
"""
import automol.vmat
from automol.graph._graph_base import string
from automol.graph._graph import atom_count
from automol.graph._graph import atom_keys
from automol.graph._graph import atom_symbols
from automol.graph._graph import remove_bonds
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import sorted_atom_neighbor_keys
from automol.graph._graph import is_connected
from automol.graph._graph import terminal_heavy_atom_keys
from automol.graph._graph import shortest_path_between_groups
from automol.graph._ring import rings
from automol.graph._ring import rings_atom_keys
from automol.graph._ring import sorted_ring_atom_keys
from automol.graph._ring import ring_systems
from automol.graph._ring import ring_system_decomposed_atom_keys
from automol.graph._ring import cycle_ring_atom_key_to_front


def vmatrix(gra, rng_keys=None):
    """ v-matrix for a connected graph

    :param gra: the graph
    :param rng_keys: keys for a ring to start from
    """
    assert is_connected(gra), "Graph must be connected!"

    rsys = sorted(ring_systems(gra), key=atom_count)

    # Start with the ring systems and their connections. If there aren't any,
    # start with the first terminal atom
    if rsys:
        vma, zma_keys = connected_ring_systems(gra, rng_keys=rng_keys)
    else:
        term_keys = sorted(terminal_heavy_atom_keys(gra))
        if term_keys:
            start_key = term_keys[0]
        else:
            start_key = sorted(atom_keys(gra))[0]

        vma, zma_keys = start_at(gra, start_key)

    # Complete any incomplete branches
    branch_keys = _atoms_missing_neighbors(gra, zma_keys)
    for key in branch_keys:
        vma, zma_keys = complete_branch(gra, key, vma, zma_keys)

    return vma, zma_keys


def connected_ring_systems(gra, rng_keys=None, check=True):
    """ generate a v-matrix covering a graph's ring systems and the connections
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
        idx = next((i for i, ks in enumerate(map(atom_keys, rsys))
                    if set(rng_keys) <= ks), None)
        assert idx is not None, (
            "The ring {} is not in this graph:\n{}"
            .format(str(rng_keys), string(gra, one_indexed=False)))
        rsy = rsys.pop(idx)

    keys_lst = list(ring_system_decomposed_atom_keys(rsy, rng_keys=rng_keys))

    vma, zma_keys = ring_system(gra, keys_lst)

    while rsys:
        # Find the next ring system with a connection to the current
        # v-vmatrix and connect them
        conn = False
        for idx, rsy_keys in enumerate(map(atom_keys, rsys)):
            if set(zma_keys) & rsy_keys:
                # ring systems are connected by one bond -- no chain needed
                keys = set(zma_keys) & rsy_keys
                assert len(keys) == 1, (
                    "Attempting to add redundant keys to v-matrix: {}"
                    .format(str(keys)))
                key, = keys

                conn = True
            else:
                # see if the ring systems are connected by a chain
                keys = shortest_path_between_groups(
                    gra, zma_keys, atom_keys(rsy))

                # if so, build a bridge from the current v-matrix to this next
                # ring system
                vma, zma_keys = continue_chain(gra, keys[:-1], vma, zma_keys,
                                               term_hydrogens=False)
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
    """ generate a v-matrix for a ring system

    :param gra: the graph
    :param keys_lst: the first entry contains keys for a ring and each next one
        contains keys for an arc that starts and ends on atoms in the preceding
        entries
    """
    # First, get the ring keys
    keys_lst = list(keys_lst)
    keys = keys_lst.pop(0)

    # Break the bonds joining the last pair of atoms in each arc
    gra = remove_bonds(gra, [(k[-1], k[-2]) for k in keys_lst])

    # Start by constructing the v-matrix for the first ring
    vma, zma_keys = ring(gra, keys)

    # Now, complete the ring system by continuing the v-matrix along each arc
    for keys in keys_lst:
        # Note that the atoms on each end of the arc are already in the
        # v-matrix, so we ignore those
        vma, zma_keys = continue_chain(gra, keys[1:-1], vma, zma_keys)

    return vma, zma_keys


def continue_ring_system(gra, keys_lst, vma, zma_keys):
    """ continue constructing a v-matrix for a ring system

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
        vma, zma_keys = continue_chain(gra, keys[1:-1], vma, zma_keys)

    return vma, zma_keys


def ring(gra, keys):
    """ generate a v-matrix for a ring

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
    """ continue constructing a v-matrix around a ring

    All neighboring atoms along the ring will be included

    Exactly one atom in the ring must already be in the v-matrix.
    """
    # Find the connecting key
    key = next((k for k in keys if k in zma_keys), None)
    assert key is not None, (
        "There must be a ring atom already in the v-matrix")

    # Cycle the connecting key to the front of the ring
    keys = cycle_ring_atom_key_to_front(keys, key)

    # Break the bond between the first and last atoms to make this a chain
    gra = remove_bonds(gra, [(keys[0], keys[-1])])

    # Now, construct a v-matrix for the chain
    vma, zma_keys = continue_chain(gra, keys, vma, zma_keys)
    return vma, zma_keys


def chain(gra, keys, term_hydrogens=True):
    """ generate a v-matrix for a chain

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

    start_key, = set(keys) & set(_atoms_missing_neighbors(gra, zma_keys))
    keys = keys[keys.index(start_key):]

    # 2. Continue the chain
    if keys:
        vma, zma_keys = continue_chain(gra, keys, vma=vma, zma_keys=zma_keys)

    return vma, zma_keys


def continue_chain(gra, keys, vma, zma_keys, term_hydrogens=True):
    """ continue constructing a v-matrix along a chain

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
        keys = _extend_chain_to_include_terminal_hydrogens(gra, keys,
                                                           start=False)

    vma, zma_keys = complete_branch(gra, keys[0], vma, zma_keys,
                                    branch_keys=keys)

    return vma, zma_keys


def start_at(gra, key):
    """ start a v-matrix at a specific atom

    Returns the started vmatrix, along with keys to atoms whose neighbors are
    missing from it
    """
    symb_dct = atom_symbols(gra)
    ngb_keys_dct = sorted_atom_neighbor_keys(
        gra, symbs_first=('X', 'C',), symbs_last=('H',), ords_last=(0.1,))

    ngb_keys = ngb_keys_dct[key]
    if not ngb_keys:
        zma_keys = []
    elif len(ngb_keys) == 1:
        # Need special handling for atoms with only one neighbor
        if symb_dct[key] in ('H', 'X'):
            key2 = ngb_keys[0]
            zma_keys = (key2,) + ngb_keys_dct[key2]
        else:
            key2 = ngb_keys[0]
            ngb_keys = tuple(k for k in ngb_keys_dct[key2] if k != key)
            zma_keys = (key, key2) + ngb_keys
    else:
        zma_keys = (key,) + ngb_keys_dct[key]

    vma = ()
    for row, key_ in enumerate(zma_keys):
        idx1 = idx2 = idx3 = None
        if row > 0:
            key1 = next(k for k in ngb_keys_dct[key_] if k in zma_keys[:row])
            idx1 = zma_keys.index(key1)
        if row > 1:
            key2 = next(k for k in ngb_keys_dct[key1] if k in zma_keys[:row]
                        and k != key_)
            idx2 = zma_keys.index(key2)
        if row > 2:
            key3 = next(k for k in zma_keys[:row]
                        if k not in (key_, key1, key2))
            idx3 = zma_keys.index(key3)

        sym = symb_dct[key_]
        key_row = [idx1, idx2, idx3]
        vma = automol.vmat.add_atom(vma, sym, key_row)

    return vma, zma_keys


def complete_branch(gra, key, vma, zma_keys, branch_keys=None):
    """ continue constructing a v-matrix along a chain

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
    ngb_keys_dct = sorted_atom_neighbor_keys(
        gra, symbs_first=('X', 'C',), symbs_last=('H',), ords_last=(0.1,))

    def _continue(key1, key2, key3, vma, zma_keys):
        k3ns = list(ngb_keys_dct[key3])
        for k3n in set(k3ns) & set(zma_keys):
            k3ns.remove(k3n)

        if k3ns:
            key4 = k3ns.pop(0)

            # Add the leading atom to the v-matrix
            vma = automol.vmat.add_atom(vma, sym, key_row)
            assert key4 not in zma_keys, ("Atom {:d} already in v-matrix."
                                          .format(key4))
            zma_keys.append(key4)

            # Add the neighbors of atom 3 (if any) to the v-matrix, decoupled
            # from atom 1 for properly decopuled torsions
            for k3n in k3ns:
                sym = symb_dct[k3n]

                if symb_dct[key4] == 'X':
                    key_row = list(map(zma_keys.index, (key3, key4, key2)))
                else:
                    key_row = list(map(zma_keys.index, (key3, key2, key4)))

                vma = automol.vmat.add_atom(vma, sym, key_row)
                assert k3n not in zma_keys, ("Atom {:d} already in v-matrix."
                                             .format(k3n))
                zma_keys.append(k3n)

            # Recursion
            if key4 in branch_keys:
                vma, zma_keys = _continue(key2, key3, key4, vma, zma_keys)

            if symb_dct[key4] == 'X':
                key2 = key4

            for k3n in k3ns:
                if k3n in branch_keys:
                    vma, zma_keys = _continue(key2, key3, k3n, vma, zma_keys)

        return vma, zma_keys

    key1, key2, key3 = keys[:3]
    vma, zma_keys = _continue(key1, key2, key3, vma, zma_keys)

    return vma, zma_keys


# helpers
def _extend_chain_to_include_anchoring_atoms(gra, keys, zma_keys):
    """ extend chain to include three atoms already specified in v-matrix

    :param gra: the graph
    :param keys: keys in the chain; the first atom should already be specified
    :param zma_keys: keys currently in the v-matrix
    """
    ngb_keys_dct = sorted_atom_neighbor_keys(
        gra, symbs_first=('X', 'C',), symbs_last=('H',), ords_last=(0.1,))

    key3 = keys[0]
    assert key3 in zma_keys
    key2 = next(k for k in ngb_keys_dct[key3] if k in zma_keys)
    key1 = next(k for k in ngb_keys_dct[key2] if k in zma_keys and k != key3)
    keys = (key1, key2,) + tuple(keys)

    return keys


def _extend_chain_to_include_terminal_hydrogens(gra, keys,
                                                start=True, end=True):
    """ extend each end of a chain to include terminal hydrogens, if any
    """
    symb_dct = atom_symbols(gra)
    atm_ngb_dct = atom_neighbor_keys(gra)

    sta_ngbs = atm_ngb_dct[keys[0]] - {keys[1]}
    end_ngbs = atm_ngb_dct[keys[-1]] - {keys[-2]}

    sta_ngb = min((k for k in sta_ngbs if symb_dct[k] == 'H'), default=None)
    end_ngb = min((k for k in end_ngbs if symb_dct[k] == 'H'), default=None)

    keys = tuple(keys)

    if start and sta_ngb is not None:
        keys = (sta_ngb,) + keys

    if end and end_ngb is not None:
        keys = keys + (end_ngb,)

    return keys


def _atoms_missing_neighbors(gra, zma_keys):
    """ get atoms from the list currently in the v-matrix with neighbors that
    are not in the v-matrix
    """
    ngb_keys_dct = atom_neighbor_keys(gra)
    keys = []
    for key in zma_keys:
        if any(k not in zma_keys for k in ngb_keys_dct[key]):
            keys.append(key)
    keys = tuple(keys)
    return keys
