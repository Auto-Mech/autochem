""" graph-based z-matrix builder
"""
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
                     _start(gra, longest_chain(gra)[0]))

    # Finally, complete any incomplete branches
    branch_keys = _atoms_missing_neighbors(gra, row_keys)
    for key in branch_keys:
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
        vma, row_keys = _continue_chain(gra, keys[:-1], vma, row_keys,
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
    vma, row_keys = _start(gra, keys[0])

    start_key, = set(keys) & set(_atoms_missing_neighbors(gra, row_keys))
    keys = keys[keys.index(start_key):]

    # 2. Continue the chain
    if keys:
        vma, row_keys = _continue_chain(gra, keys, vma=vma,
                                        row_keys=row_keys)

    return vma, row_keys


def _continue_chain(gra, keys, vma, row_keys, term_hydrogens=True):
    """ continue constructing a v-matrix along a chain

    All neighboring atoms along the chain will be included

    Exactly one atom in the chain must already be in the v-matrix

    :param gra: the graph for which the v-matrix will be constructed
    :param keys: the keys for atoms along the chain, which must be contiguous;
        the first atom must already appear in the v-matrix
    :param vma: a partial v-matrix from which to continue
    :param row_keys: row keys for the partial v-matrix, identifying the atom
        specified by each row of `vma` in order
    :param extend: whether to extend the chain's start to include the three
        anchoring atoms
    :param term_hydrogens: whether to extend the chain's end to include
        terminal hydrogens
    """
    if term_hydrogens:
        keys = _extend_chain_to_include_terminal_hydrogens(gra, keys,
                                                           start=False)

    vma, row_keys = _complete_branch(gra, keys[0], vma, row_keys,
                                     branch_keys=keys)

    return vma, row_keys


def _extend_chain_to_include_anchoring_atoms(gra, keys, row_keys):
    """ extend chain to include three atoms already specified in v-matrix

    :param gra: the graph
    :param keys: keys in the chain; the first atom should already be specified
    :param row_keys: keys currently in the v-matrix
    """
    ngb_keys_dct = sorted_atom_neighbor_keys(gra, syms_first=('X', 'C',),
                                             syms_last=('H',))

    key3 = keys[0]
    assert key3 in row_keys
    key2 = next(k for k in ngb_keys_dct[key3] if k in row_keys)
    key1 = next(k for k in ngb_keys_dct[key2] if k in row_keys and k != key3)
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


def _start(gra, key):
    """ start a v-matrix at a specific atom

    Returns the started vmatrix, along with keys to atoms whose neighbors are
    missing from it
    """
    sym_dct = atom_symbols(gra)
    ngb_keys_dct = sorted_atom_neighbor_keys(gra, syms_first=('X', 'C',),
                                             syms_last=('H',))

    ngb_keys = ngb_keys_dct[key]
    if not ngb_keys:
        row_keys = []
    elif len(ngb_keys) == 1:
        # Need special handling for atoms with only one neighbor
        if sym_dct[key] in ('H', 'X'):
            key2 = ngb_keys[0]
            row_keys = (key2,) + ngb_keys_dct[key2]
        else:
            key2 = ngb_keys[0]
            ngb_keys = tuple(k for k in ngb_keys_dct[key2] if k != key)
            row_keys = (key, key2) + ngb_keys
    else:
        row_keys = (key,) + ngb_keys_dct[key]

    vma = ()
    for row, key_ in enumerate(row_keys):
        idx1 = idx2 = idx3 = None
        if row > 0:
            key1 = next(k for k in ngb_keys_dct[key_] if k in row_keys[:row])
            idx1 = row_keys.index(key1)
        if row > 1:
            key2 = next(k for k in ngb_keys_dct[key1] if k in row_keys[:row]
                        and k != key_)
            idx2 = row_keys.index(key2)
        if row > 2:
            key3 = next(k for k in row_keys[:row]
                        if k not in (key_, key1, key2))
            idx3 = row_keys.index(key3)

        sym = sym_dct[key_]
        key_row = [idx1, idx2, idx3]
        vma = automol.vmat.add_atom(vma, sym, key_row)

    return vma, row_keys


def _atoms_missing_neighbors(gra, row_keys):
    """ get atoms from the list currently in the v-matrix with neighbors that
    are not in the v-matrix
    """
    ngb_keys_dct = atom_neighbor_keys(gra)
    keys = []
    for key in row_keys:
        if any(k not in row_keys for k in ngb_keys_dct[key]):
            keys.append(key)
    keys = tuple(keys)
    return keys


def _complete_branch(gra, key, vma, row_keys, branch_keys=None):
    """ continue constructing a v-matrix along a chain

    All neighboring atoms along the chain will be included

    Exactly one atom in the chain must already be in the v-matrix

    :param gra: the graph for which the v-matrix will be constructed
    :param keys: the keys for atoms along the chain, which must be contiguous;
        the first atom must already appear in the v-matrix
    :param vma: a partial v-matrix from which to continue
    :param row_keys: row keys for the partial v-matrix, identifying the atom
        specified by each row of `vma` in order
    :param branch_keys: optionally, restrict the v-matrix to these keys and
        their neighbors; if `None`, the entire branch will be included
    """
    branch_keys = atom_keys(gra) if branch_keys is None else branch_keys
    keys = _extend_chain_to_include_anchoring_atoms(gra, [key], row_keys)

    row_keys = list(row_keys)
    sym_dct = atom_symbols(gra)
    ngb_keys_dct = sorted_atom_neighbor_keys(gra, syms_first=('X', 'C',),
                                             syms_last=('H',))

    def _continue(key1, key2, key3, vma, row_keys):
        k3ns = list(ngb_keys_dct[key3])
        for k3n in set(k3ns) & set(row_keys):
            k3ns.remove(k3n)

        if k3ns:
            key4 = k3ns.pop(0)

            # Add the leading atom to the v-matrix
            sym = sym_dct[key4]
            key_row = list(map(row_keys.index, (key3, key2, key1)))
            vma = automol.vmat.add_atom(vma, sym, key_row)
            assert key4 not in row_keys, ("Atom {:d} already in v-matrix."
                                          .format(key4))
            row_keys.append(key4)

            # Add the neighbors of atom 3 (if any) to the v-matrix, decoupled
            # from atom 1 for properly decopuled torsions
            for k3n in k3ns:
                sym = sym_dct[k3n]

                if sym_dct[key4] == 'X':
                    key_row = list(map(row_keys.index, (key3, key4, key2)))
                else:
                    key_row = list(map(row_keys.index, (key3, key2, key4)))

                vma = automol.vmat.add_atom(vma, sym, key_row)
                assert k3n not in row_keys, ("Atom {:d} already in v-matrix."
                                             .format(k3n))
                row_keys.append(k3n)

            # Recursion
            if key4 in branch_keys:
                vma, row_keys = _continue(key2, key3, key4, vma, row_keys)

            if sym_dct[key4] == 'X':
                key2 = key4

            for k3n in k3ns:
                if k3n in branch_keys:
                    vma, row_keys = _continue(key2, key3, k3n, vma, row_keys)

        return vma, row_keys

    key1, key2, key3 = keys[:3]
    vma, row_keys = _continue(key1, key2, key3, vma, row_keys)

    return vma, row_keys


if __name__ == '__main__':
    import automol
    # ICH = automol.smiles.inchi('CC(C)C#C')
    # ICH = automol.smiles.inchi('CCCC(OO)CC(CC(N)(CC)CC)C=C=CC#C')
    # ICH = automol.smiles.inchi('C1CCCC2C1.C2C3.C4C3CCC4')
    # ICH = automol.smiles.inchi('C1CCC(CCC2CCCC2)CC1')
    # ICH = automol.smiles.inchi('C12C(OON)C3C(CC2)CC1'
    #                            '.C3C#CC(C(C)C)C4'
    #                            '.C45C(CC6)CC(CCO)C56')
    # ICH = automol.smiles.inchi('C1CCCCC1')
    # ICH = automol.smiles.inchi('C#CCCCC#CCCCC#C')
    # ICH = automol.smiles.inchi('C=C=C')
    # ICH = automol.smiles.inchi('C#C')
    ICH = 'InChI=1S/C3H7O4/c1-3(7-5)2-6-4/h3-4H,2H2,1H3/t3-/m0/s1'
    GEO = automol.inchi.geometry(ICH)
    # # Yuri's code:
    # ZMA = automol.geom.zmatrix(GEO)
    # print(automol.zmat.string(ZMA, one_indexed=False))
    # print(automol.geom.zmatrix_torsion_coordinate_names(GEO))
    # GEO = automol.zmat.geometry(ZMA)
    # My code:
    GEO = automol.geom.insert_dummies_on_linear_atoms(GEO)
    GRA = automol.geom.connectivity_graph(GEO, dummy_bonds=True)
    print(automol.geom.string(GEO))
    print(automol.graph.string(GRA, one_indexed=False))
    # KEYS = longest_chain(GRA)
    # VMA, ROW_KEYS = _start(GRA, KEYS[0])
    VMA, ROW_KEYS = vmatrix(GRA)
    print(automol.vmat.string(VMA, one_indexed=False))
    SUBGEO = automol.geom.from_subset(GEO, ROW_KEYS)
    SUBZMA = automol.zmat.from_geometry(VMA, SUBGEO)
    print(automol.zmat.string(SUBZMA, one_indexed=False))
    SUBGEO = automol.zmat.geometry(SUBZMA)
    SUBGEO = automol.geom.mass_centered(SUBGEO)
    SUBGEO = automol.geom.without_dummy_atoms(SUBGEO)
    print(automol.geom.string(SUBGEO))
    ICH_OUT = automol.geom.inchi(SUBGEO)
    print(ICH_OUT)
    assert ICH == ICH_OUT
