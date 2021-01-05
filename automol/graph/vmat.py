""" graph-based z-matrix builder
"""
import more_itertools as mit
import automol.vmat
from automol.graph._graph import atom_count
from automol.graph._graph import atom_keys
from automol.graph._graph import atom_symbols
from automol.graph._graph import remove_bonds
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import sorted_atom_neighbor_keys
from automol.graph._graph import is_connected
from automol.graph._graph import longest_chain
from automol.graph._graph import shortest_path_between_groups
from automol.graph._ring import rings
from automol.graph._ring import ring_systems
from automol.graph._ring import ring_system_decomposed_atom_keys
from automol.graph._ring import cycle_ring_atom_key_to_front


def vmatrix(gra):
    """ v-matrix for a connected graph

    :param gra: the graph
    :param lin_keys: keys for linear atoms in the graph (dummy atoms will be
        inserted)
    """
    assert is_connected(gra), "Graph must be connected!"

    rsys = sorted(ring_systems(gra), key=atom_count)

    # Start with the ring systems and their connections. If there aren't any,
    # start with the longest chain.
    vma, row_keys = (_connected_ring_systems(gra, check=False) if rsys else
                     _start_chain(gra, longest_chain(gra)))

    # Finally, complete any incomplete branches
    ngb_keys_dct = sorted_atom_neighbor_keys(gra)
    for key in row_keys:
        # If any neighbor keys are not in the v-matrix, this is an incomplete
        # branch -- complete it.
        if any(k not in row_keys for k in ngb_keys_dct[key]):
            vma, row_keys = _complete_branch(gra, key, vma, row_keys)

    return vma, row_keys


def _connected_ring_systems(gra, check=True):
    """ generate a v-matrix covering a graph's ring systems and the connections
    between them
    """
    if check:
        assert is_connected(gra), "Graph must be connected!"

    rsys = sorted(ring_systems(gra), key=atom_count)

    # Construct the v-matrix for the first ring system, choosing which ring
    # to start from
    rsy = rsys.pop(0)
    rngs = sorted(rings(rsy), key=atom_count)
    rng = rngs.pop(0)
    keys_lst = list(ring_system_decomposed_atom_keys(rsy, rng=rng))

    vma, row_keys = _ring_system(gra, keys_lst)

    while rsys:
        # Find the next ring system with a connection to the current
        # v-vmatrix
        for idx, rsy in enumerate(rsys):
            keys = shortest_path_between_groups(
                gra, row_keys, atom_keys(rsy))
            if keys is not None:
                rsys.pop(idx)
                break

        assert keys is not None, "This is a disconnected graph!"

        # 1. Build a bridge from the current v-matrix to the next ring
        # system
        vma, row_keys = _continue_chain(gra, keys, vma, row_keys,
                                        term_hydrogens=False)

        # 2. Decompose the ring system with the connecting ring first
        rng = next(r for r in rings(rsy) if keys[-1] in atom_keys(r))
        keys_lst = ring_system_decomposed_atom_keys(rsy, rng=rng)

        # 3. Build the next ring system
        vma, row_keys = _continue_ring_system(gra, keys_lst, vma, row_keys)

    return vma, row_keys


def _ring_system(gra, keys_lst):
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
    vma, row_keys = _ring(gra, keys)

    # Now, complete the ring system by continuing the v-matrix along each arc
    for keys in keys_lst:
        # Note that the atoms on each end of the arc are already in the
        # v-matrix, so we ignore those
        vma, row_keys = _continue_chain(gra, keys[1:-1], vma, row_keys)

    return vma, row_keys


def _continue_ring_system(gra, keys_lst, vma, row_keys):
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
    vma, row_keys = _continue_ring(gra, keys, vma, row_keys)

    # Now, complete the ring system by continuing the v-matrix along each arc
    for keys in keys_lst:
        # Note that the atoms on each end of the arc are already in the
        # v-matrix, so we ignore those
        vma, row_keys = _continue_chain(gra, keys[1:-1], vma, row_keys)

    return vma, row_keys


def _ring(gra, keys):
    """ generate a v-matrix for a ring

    All neighboring atoms along the ring will be included

    :param gra: the graph
    :param keys: ring keys, in the order they should appear in the z-matrix
    """
    # Break the bond between the first and last atoms to make this a chain
    gra = remove_bonds(gra, [(keys[0], keys[-1])])

    # Now, construct a v-matrix for the chain
    vma, row_keys = _chain(gra, keys, term_hydrogens=True)
    return vma, row_keys


def _continue_ring(gra, keys, vma, row_keys):
    """ continue constructing a v-matrix around a ring

    All neighboring atoms along the ring will be included

    Exactly one atom in the ring must already be in the v-matrix.
    """
    # Find the connecting key
    key = next((k for k in keys if k in row_keys), None)
    assert key is not None, (
        "There must be a ring atom already in the v-matrix")

    # Cycle the connecting key to the front of the ring
    keys = cycle_ring_atom_key_to_front(keys, key)

    # Break the bond between the first and last atoms to make this a chain
    gra = remove_bonds(gra, [(keys[0], keys[-1])])

    # Now, construct a v-matrix for the chain
    vma, row_keys = _continue_chain(gra, keys, vma, row_keys)
    return vma, row_keys


def _chain(gra, keys, term_hydrogens=True):
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
    vma, row_keys = _start_chain(gra, keys)

    # 2. Continue the chain
    if len(keys) > 3:
        vma, row_keys = _continue_chain(gra, keys, vma=vma, row_keys=row_keys,
                                        extend=False)

    return vma, row_keys


def _complete_branch(gra, key, vma, row_keys):
    """ finish constructing the v-matrix for a branch
    """
    keys = _extend_chain_to_include_anchoring_atoms(gra, [key], row_keys)

    sym_dct = atom_symbols(gra)
    ngb_keys_dct = sorted_atom_neighbor_keys(gra, syms_first=('C',),
                                             syms_last=('H',))

    def _recurse(key1, key2, key3, vma, row_keys):
        k3ns = list(ngb_keys_dct[key3])
        k3ns.remove(key2)

        if k3ns:
            # Add the leading atom to the v-matrix
            key4 = k3ns.pop(0)
            sym = sym_dct[key4]
            key_row = list(map(row_keys.index, (key3, key2, key1)))
            vma = automol.vmat.add_atom(vma, sym, key_row)
            assert key4 not in row_keys, ("Atom {:d} already in v-matrix."
                                          .format(key4))
            row_keys.append(key4)

            for k3n in k3ns:
                sym = sym_dct[k3n]
                key_row = list(map(row_keys.index, (key3, key2, key4)))
                vma = automol.vmat.add_atom(vma, sym, key_row)
                assert k3n not in row_keys, ("Atom {:d} already in v-matrix."
                                             .format(k3n))
                row_keys.append(k3n)

            # Recursion
            for k3n in [key4] + k3ns:
                vma, row_keys = _recurse(key2, key3, k3n, vma, row_keys)

        return vma, row_keys

    key1, key2, key3 = keys
    vma, row_keys = _recurse(key1, key2, key3, vma, row_keys)

    return vma, row_keys


def _continue_chain(gra, keys, vma, row_keys, extend=True,
                    term_hydrogens=True):
    """ continue constructing a v-matrix along a chain

    All neighboring atoms along the chain will be included

    Exactly one atom in the chain must already be in the v-matrix

    :param gra: the graph for which the v-matrix will be constructed
    :param keys: the keys for atoms along the chain, which must be contiguous;
        if `extend` is False, the first three atoms in this list must already
        be specified in the v-matrix
    :param vma: a partial v-matrix from which to continue
    :param row_keys: row keys for the partial v-matrix, identifying the atom
        specified by each row of `vma` in order (None indicates an inserted
        dummy atom)
    :param extend: whether to extend the chain's start to include the three
        anchoring atoms
    :param term_hydrogens: whether to extend the chain's end to include
        terminal hydrogens
    """
    if extend:
        keys = _extend_chain_to_include_anchoring_atoms(gra, keys, row_keys)
    if term_hydrogens:
        keys = _extend_chain_to_include_terminal_hydrogens(gra, keys,
                                                           start=False)

    sym_dct = atom_symbols(gra)
    ngb_keys_dct = sorted_atom_neighbor_keys(gra, syms_first=('C',),
                                             syms_last=('H',))

    assert set(keys[:3]) <= set(row_keys), (
        "The first three atoms must already be specified.")
    assert not set(keys[3:]) & set(row_keys), (
        "Cannot add atoms {:s} that are already in the z-matrix."
        .format(str(set(keys[3:]) & set(row_keys))))

    for key1, key2, key3, key4 in mit.windowed(keys, 4):
        sym = sym_dct[key4]

        # Add the atom 4 to the v-matrix
        key_row = list(map(row_keys.index, (key3, key2, key1)))
        vma = automol.vmat.add_atom(vma, sym, key_row)
        assert key4 not in row_keys, ("Atom {:d} already in v-matrix."
                                      .format(key4))
        row_keys.append(key4)

        # Add the neighbors of atom 3 (if any) to the v-matrix, decoupled from
        # atom 1 for properly decopuled torsions
        k3ns = list(ngb_keys_dct[key3])
        k3ns.remove(key2)
        k3ns.remove(key4)
        for k3n in k3ns:
            sym = sym_dct[k3n]
            key_row = list(map(row_keys.index, (key3, key2, key4)))
            vma = automol.vmat.add_atom(vma, sym, key_row)
            assert k3n not in row_keys, ("Atom {:d} already in v-matrix."
                                         .format(k3n))
            row_keys.append(k3n)

    return vma, row_keys


def _start_chain(gra, keys, term_hydrogens=True):
    """ start v-matrix for a chain

    (Seems to be relatively robust)
    """
    if term_hydrogens:
        keys = _extend_chain_to_include_terminal_hydrogens(gra, keys,
                                                           end=False)

    sym_dct = atom_symbols(gra)
    ngb_keys_dct = sorted_atom_neighbor_keys(gra, syms_first=('C',),
                                             syms_last=('H',))

    if len(keys) > 2:
        ngb_keys = list(ngb_keys_dct[keys[1]])
        ngb_keys.remove(keys[0])
        ngb_keys.remove(keys[2])

        if sym_dct[keys[0]] == 'H':
            row_keys = [keys[1], keys[2], keys[0]] + ngb_keys
        else:
            row_keys = [keys[0], keys[1], keys[2]] + ngb_keys
    else:
        row_keys = keys[:3]

    syms = tuple(map(sym_dct.__getitem__, row_keys))

    vma = ()
    for row, sym in enumerate(syms):
        key_row = [0, 1, 2][:row]
        vma = automol.vmat.add_atom(vma, sym, key_row)

    return vma, row_keys


def _extend_chain_to_include_anchoring_atoms(gra, keys, row_keys):
    """ extend chain to include three atoms already specified in v-matrix

    :param gra: the graph
    :param keys: keys in the chain; the first atom should already be specified
    :param row_keys: keys currently in the v-matrix
    """
    atm_ngb_dct = sorted_atom_neighbor_keys(gra)

    key3 = keys[0]
    assert key3 in row_keys
    key2 = next(k for k in atm_ngb_dct[key3] if k in row_keys)
    key1 = next(k for k in atm_ngb_dct[key2] if k in row_keys and k != key3)
    keys = (key1, key2,) + tuple(keys)

    return keys


def _extend_chain_to_include_terminal_hydrogens(gra, keys,
                                                start=True, end=True):
    """ extend each end of a chain to include terminal hydrogens, if any
    """
    sym_dct = atom_symbols(gra)
    atm_ngb_dct = atom_neighbor_keys(gra)

    sta_ngbs = atm_ngb_dct[keys[0]] - {keys[1]}
    end_ngbs = atm_ngb_dct[keys[-1]] - {keys[-2]}

    sta_ngb = min((k for k in sta_ngbs if sym_dct[k] == 'H'), default=None)
    end_ngb = min((k for k in end_ngbs if sym_dct[k] == 'H'), default=None)

    keys = tuple(keys)

    if start and sta_ngb is not None:
        keys = (sta_ngb,) + keys

    if end and end_ngb is not None:
        keys = keys + (end_ngb,)

    return keys


if __name__ == '__main__':
    import automol
    # ICH = automol.smiles.inchi('CC(C)C')
    # ICH = automol.smiles.inchi('CCCC(OO)CC(CC(N)(CC)CC)CCCCC')
    # ICH = automol.smiles.inchi('C1CCC(CCC2CCCC2)CC1')
    ICH = automol.smiles.inchi('C12C(OON)C3C(CC2)CC1'
                               '.C3C(C(C)C)C4'
                               '.C45C(CC6)CC(CCO)C56')
    # ICH = automol.smiles.inchi('C#CCCCC#CCCCC#C')
    GEO = automol.inchi.geometry(ICH)
    print(automol.geom.string(GEO))
    GRA = automol.geom.connectivity_graph(GEO)
    # Yuri's code:
    ZMA = automol.geom.zmatrix(GEO)
    print(automol.zmatrix.string(ZMA, one_indexed=False))
    print(automol.geom.zmatrix_torsion_coordinate_names(GEO))
    # My code:
    VMA, ROW_KEYS = vmatrix(GRA)
    SUBGEO = automol.geom.from_subset(GEO, ROW_KEYS)
    SUBZMA = automol.zmat.from_geometry(VMA, SUBGEO)
    SUBGEO = automol.zmat.geometry(SUBZMA)
    SUBGEO = automol.geom.mass_centered(SUBGEO)
    print(automol.zmat.string(SUBZMA, one_indexed=False))
    print(automol.geom.string(SUBGEO))
