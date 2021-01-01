""" graph-based z-matrix builder
"""
import more_itertools as mit
import automol.vmat
from automol.graph._graph import atom_symbols
from automol.graph._graph import sorted_atom_neighbor_keys
from automol.graph._graph import longest_chain


def main(geo):
    """ v-matrix for a chain of heavy atoms
    """
    gra = automol.geom.connectivity_graph(geo)
    keys = list(longest_chain(gra))

    # 1. Start a chain
    vma, row_keys = _start_chain(gra, keys)

    # 2. Continue the chain
    if len(keys) > 3:
        vma, row_keys = _continue_chain(gra, keys, vma=vma, row_keys=row_keys)

    subgeo = automol.geom.from_subset(geo, row_keys)
    subzma = automol.zmat.from_geometry(vma, subgeo)
    subgeo = automol.zmat.geometry(subzma)
    print(automol.zmat.string(subzma))
    print(automol.geom.string(subgeo))


def _continue_chain(gra, keys, vma, row_keys):
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
    sym_dct = atom_symbols(gra)
    ngb_keys_dct = sorted_atom_neighbor_keys(gra, syms_first=('C',),
                                             syms_last=('H',))

    for key1, key2, key3, key4 in mit.windowed(keys, 4):
        sym = sym_dct[key4]

        # Add the atom 4 to the v-matrix
        key_row = list(map(row_keys.index, (key3, key2, key1)))
        vma = automol.vmat.add_atom(vma, sym, key_row)
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
            row_keys.append(k3n)

    return vma, row_keys


def _start_chain(gra, keys):
    """ start v-matrix for a chain

    (Seems to be relatively robust)
    """
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


if __name__ == '__main__':
    import automol
    ICH = automol.smiles.inchi('CCC(CC)CCC')
    GEO = automol.inchi.geometry(ICH)
    ZMA = automol.geom.zmatrix(GEO)
    print(automol.zmatrix.string(ZMA, one_indexed=False))
    print(automol.geom.zmatrix_torsion_coordinate_names(GEO))
    main(GEO)
