""" functions for working with rotational bonds and groups
"""
import itertools
from automol.graph._graph_base import atom_symbols
from automol.graph._graph_base import bond_keys
from automol.graph._graph import explicit
from automol.graph._graph import branch_atom_keys
from automol.graph._graph import atom_neighbor_keys
from automol.graph._ring import rings_bond_keys
from automol.graph._res import resonance_dominant_bond_orders


def rotational_bond_keys(gra, with_h_rotors=True):
    """ get all rotational bonds for a graph
    """
    gra = explicit(gra)
    sym_dct = atom_symbols(gra)
    ngb_keys_dct = atom_neighbor_keys(gra)
    bnd_ord_dct = resonance_dominant_bond_orders(gra)
    rng_bnd_keys = list(itertools.chain(*rings_bond_keys(gra)))

    def _is_rotational_bond(bnd_key):
        ngb_keys_lst = [ngb_keys_dct[k] - bnd_key for k in bnd_key]

        is_single = max(bnd_ord_dct[bnd_key]) <= 1
        has_neighbors = all(ngb_keys_lst)
        not_in_ring = bnd_key not in rng_bnd_keys

        is_h_rotor = any(set(map(sym_dct.__getitem__, ks)) == {'H'}
                         for ks in ngb_keys_lst)
        return is_single and has_neighbors and not_in_ring and (
            not is_h_rotor or with_h_rotors)

    rot_bnd_keys = frozenset(filter(_is_rotational_bond, bond_keys(gra)))
    return rot_bnd_keys


def rotational_groups(gra, key1, key2):
    """ get the rotational groups for a given rotational axis

    :param gra: the graph
    :param key1: the first atom key
    :param key2: the second atom key
    """
    bnd_key = frozenset({key1, key2})
    assert bnd_key in bond_keys(gra)
    grp1 = branch_atom_keys(gra, key2, bnd_key) - {key1}
    grp2 = branch_atom_keys(gra, key1, bnd_key) - {key2}
    grp1 = tuple(sorted(grp1))
    grp2 = tuple(sorted(grp2))
    return grp1, grp2


if __name__ == '__main__':
    import automol
    ICH = automol.smiles.inchi('CCCC')
    # ICH = automol.smiles.inchi('C1CCCC1')
    GEO = automol.inchi.geometry(ICH)
    ZMA = automol.convert.geom.zmatrix_x2z(GEO)
    GEO = automol.zmat.geometry(ZMA)
    ZMA = automol.convert.geom.zmatrix_x2z(GEO)
    NAMES = automol.geom.zmatrix_torsion_coordinate_names(GEO)
    COO_DCT = automol.vmat.coordinates(ZMA)
    COOS = [x[0][1:-1] for x in map(COO_DCT.__getitem__, NAMES)]
    COOS = set(map(frozenset, COOS))
    print(automol.zmat.string(ZMA))
    print(NAMES)
    print(COOS)
    print()

    GRA = automol.geom.graph(GEO)
    ROT_BND_KEYS = rotational_bond_keys(GRA)
    print(len(COOS))
    print(len(ROT_BND_KEYS))
    print(ROT_BND_KEYS)
    print(COOS == ROT_BND_KEYS)
