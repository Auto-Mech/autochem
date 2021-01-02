""" graph-based z-matrix builder
"""
import more_itertools as mit
import automol.vmat
from automol.graph._graph import atom_symbols
from automol.graph._graph import remove_bonds
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import sorted_atom_neighbor_keys
from automol.graph._graph import longest_chain
from automol.graph._graph import shortest_path_between_groups
from automol.graph._ring import rings_atom_keys
from automol.graph._ring import cycle_ring_atom_key_to_front


def main_ring(geo):
    """ v-matrix for a chain of heavy atoms

    works for rings connected to each other
    """
    gra = automol.geom.connectivity_graph(geo)
    rng_keys_lst = sorted(rings_atom_keys(gra), key=len, reverse=True)

    keys = rng_keys_lst.pop(0)

    auxg = remove_bonds(gra, [(keys[0], keys[-1])])

    vma, row_keys = _chain(auxg, keys)

    print('ring 1 keys', keys)
    print('row keys', row_keys)

    # Connect all rings:
    # (Eventually this will need to be done for ring *systems* instead)
    while rng_keys_lst:
        keys = None

        # Find the next ring with a connection to the current v-matrix and
        # remove it from the list
        rng_keys = None
        for idx, rng_keys in enumerate(rng_keys_lst):
            keys = shortest_path_between_groups(gra, row_keys, rng_keys)
            if keys is not None:
                rng_keys_lst.pop(idx)
                break

        assert keys is not None, "This is a disconnected graph!"

        # 1. Build the bridge between the current v-matrix and the next ring
        keys = _extend_chain_to_include_anchoring_atoms(auxg, keys, row_keys)
        vma, row_keys = _continue_chain(auxg, keys, vma, row_keys,
                                        term_hydrogens=False)
        print('connecting keys', keys)
        print('row keys', row_keys)

        # 2. Build the next ring
        conn_key = next(k for k in rng_keys if k in row_keys)
        print('ring 2 keys', rng_keys)
        keys = cycle_ring_atom_key_to_front(rng_keys, conn_key)
        print('ring 2 keys', keys)
        auxg = remove_bonds(auxg, [(keys[0], keys[-1])])
        keys = _extend_chain_to_include_anchoring_atoms(auxg, keys, row_keys)
        print('ring 2 keys', keys)
        print('row keys', row_keys)
        vma, row_keys = _continue_chain(auxg, keys, vma, row_keys)

    subgeo = automol.geom.from_subset(geo, row_keys)
    subzma = automol.zmat.from_geometry(vma, subgeo)
    subgeo = automol.zmat.geometry(subzma)
    print(automol.zmat.string(subzma))
    print(automol.geom.string(subgeo))


def main_linear(geo):
    """ v-matrix for a chain of heavy atoms

    works for linear chains
    """
    gra = automol.geom.connectivity_graph(geo)
    keys = list(longest_chain(gra))

    # 1. Start a chain
    vma, row_keys = _chain(gra, keys)

    subgeo = automol.geom.from_subset(geo, row_keys)
    subzma = automol.zmat.from_geometry(vma, subgeo)
    subgeo = automol.zmat.geometry(subzma)
    print(automol.zmat.string(subzma))
    print(automol.geom.string(subgeo))


def _chain(gra, keys, term_hydrogens=True):
    """ generate a v-matrix for a chain

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
        vma, row_keys = _continue_chain(gra, keys, vma=vma, row_keys=row_keys)

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


def _continue_chain(gra, keys, vma, row_keys, term_hydrogens=True):
    """ continue constructing a v-matrix along a chain

    All neighboring atoms along the chain will be included

    :param gra: the graph for which the v-matrix will be constructed
    :param keys: the keys for atoms along the chain, which must be contiguous;
        the first three atoms in this list must already be specified in the
        v-matrix
    :param vma: a partial v-matrix from which to continue
    :param vma_keys: row keys for the partial v-matrix, identifying the atom
        specified by each row of `vma` in order (None indicates an inserted
        dummy atom)
    """
    if term_hydrogens:
        keys = _extend_chain_to_include_terminal_hydrogens(gra, keys,
                                                           start=False)

    sym_dct = atom_symbols(gra)
    ngb_keys_dct = sorted_atom_neighbor_keys(gra, syms_first=('C',),
                                             syms_last=('H',))

    assert set(keys[:3]) <= set(row_keys), (
        "The first three atoms must already be specified.")
    print('KEYS', keys)
    print('ROW KEYS', row_keys)
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
    ICH = automol.smiles.inchi('C1CCC(CCC2CCCC2)CC1')
    GEO = automol.inchi.geometry(ICH)
    ZMA = automol.geom.zmatrix(GEO)
    print(automol.zmatrix.string(ZMA, one_indexed=False))
    print(automol.geom.zmatrix_torsion_coordinate_names(GEO))
    # main_linear(GEO)
    main_ring(GEO)
