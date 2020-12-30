""" graph-based z-matrix builder
"""
import automol.vmat
from automol.graph._graph import atom_symbols
from automol.graph._graph import sorted_atom_neighbor_keys
from automol.graph._graph import longest_chain


def main(geo):
    """ v-matrix for a chain of heavy atoms
    """
    gra = automol.geom.connectivity_graph(geo)
    chain_keys = list(longest_chain(gra))

    # 1. start from a given atom
    vma, row_keys = _start_from_atom(gra, chain_keys[1])

    # 2. Fill in a branch
    _chain(vma, row_keys, gra, chain_keys[3:])

    # sub_geo = automol.geom.from_subset(geo, row_keys)
    # zma = automol.zmat.from_geometry(vma, sub_geo)
    # geo2 = automol.zmat.geometry(zma)
    # print(automol.zmat.string(zma, one_indexed=False))
    # print(automol.geom.string(geo2))


def _start_from_atom(gra, start_key):
    """ start a z-matrix from a specific atom
    """
    sym_dct = atom_symbols(gra)
    ngb_keys_dct = sorted_atom_neighbor_keys(gra, syms_first=('C',),
                                             syms_last=('H',))
    ngb_keys = ngb_keys_dct[start_key]

    row_keys = (start_key,) + ngb_keys
    syms = tuple(map(sym_dct.__getitem__, row_keys))

    vma = ()
    for row, sym in enumerate(syms):
        key_row = [0, 1, 2][:row]
        vma = automol.vmat.add_atom(vma, sym, key_row)

    return vma, row_keys


def _key_row(row_keys, key, ngb_keys_dct):
    """ find three connected anchoring atoms
    """
    key1 = next(k for k in ngb_keys_dct[key] if k in row_keys)
    key2 = next(k for k in ngb_keys_dct[key1] if k in row_keys)
    key3 = next(k for k in ngb_keys_dct[key2] if k in row_keys and k != key1)
    keys = (key1, key2, key3)
    key_row = tuple(map(row_keys.index, keys))
    return key_row


def _chain(vma, row_keys, gra, chain_keys):
    """ continue a z-matrix along a chain

    (fills in hydrogens, but does not extend branches)
    """
    sym_dct = atom_symbols(gra)
    ngb_keys_dct = sorted_atom_neighbor_keys(gra, syms_first=('C',),
                                             syms_last=('H',))

    for key in chain_keys:
        vma = automol.vmat.add_atom(
            vma, sym=sym_dct[key],
            key_row=_key_row(row_keys, key, ngb_keys_dct))
        row_keys += (key,)
        hkeys = tuple(k for k in ngb_keys_dct[key] if sym_dct[k] == 'H')

        for hkey in hkeys:
            vma = automol.vmat.add_atom(
                vma, sym=sym_dct[hkey],
                key_row=_key_row(row_keys, hkey, ngb_keys_dct))
            row_keys += (hkey,)

        print(automol.vmat.string(vma, one_indexed=False))
        import sys
        sys.exit()
    # sym_dct = atom_symbols(gra)
    # for key in chain_keys:
    #     import sys
    #     sys.exit()


if __name__ == '__main__':
    import automol
    ICH = automol.smiles.inchi('CCC')
    GEO = automol.inchi.geometry(ICH)
    ZMA = automol.geom.zmatrix(GEO)
    print(automol.zmatrix.string(ZMA, one_indexed=False))
    print(automol.geom.zmatrix_torsion_coordinate_names(GEO))
    main(GEO)
