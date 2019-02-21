""" graph theory library
"""
from itertools import chain as _chain
from ._networkx import from_graph as _nxg_from_graph
from ._networkx import minimal_rings_atom_keys as _nxg_minimal_rings_atom_keys
from ._dict import by_key as _by_key
from ._dict import values_by_key as _values_by_key
from ._dict import transform_values as _transform_values
from ._dict import transform_items_to_values as _transform_items_to_values
from ._core import from_atoms_and_bonds as _from_atoms_and_bonds
from ._core import atoms as _atoms
from ._core import bonds as _bonds
from ._core import atom_keys as _atom_keys
from ._core import bond_keys as _bond_keys
from ._core import frozen as _frozen


# atom properties
def atom_neighbor_keys(xgr):
    """ keys of neighboring atoms, by atom
    """
    def _neighbor_keys(atm_key, atm_nbh):
        return frozenset(_atom_keys(atm_nbh) - {atm_key})

    atm_ngb_keys_dct = _transform_items_to_values(atom_neighborhoods(xgr),
                                                  _neighbor_keys)
    return atm_ngb_keys_dct


def atom_bond_keys(xgr):
    """ bond keys, by atom
    """
    return _transform_values(atom_neighborhoods(xgr), _bond_keys)


def atom_neighborhoods(xgr):
    """ neighborhood subgraphs, by atom
    """
    bnd_keys = _bond_keys(xgr)

    def _neighborhood(atm_key):
        nbh_bnd_keys = set(filter(lambda x: atm_key in x, bnd_keys))
        return bond_induced_subgraph(xgr, nbh_bnd_keys)

    atm_keys = list(_atom_keys(xgr))
    atm_nbh_dct = dict(zip(atm_keys, map(_neighborhood, atm_keys)))
    return atm_nbh_dct


# bond properties
def bond_neighbor_keys(xgr):
    """ keys of neighboring bonds, by bond
    """
    def _neighbor_keys(bnd_key, bnd_nbh):
        return frozenset(_bond_keys(bnd_nbh) - {bnd_key})

    bnd_ngb_keys_dct = _transform_items_to_values(bond_neighborhoods(xgr),
                                                  _neighbor_keys)
    return bnd_ngb_keys_dct


def bond_neighborhoods(xgr):
    """ neighborhood subgraphs, by bond
    """
    bnd_keys = list(_bond_keys(xgr))

    def _neighborhood(bnd_key):
        nbh_bnd_keys = set(filter(lambda x: bnd_key & x, bnd_keys))
        return bond_induced_subgraph(xgr, nbh_bnd_keys)

    bnd_nbh_dct = dict(zip(bnd_keys, map(_neighborhood, bnd_keys)))
    return bnd_nbh_dct


# other properties
def branch(xgr, atm_key, bnd_key):
    """ branch extending along `bnd_key` away from `atm_key`
    """
    return bond_induced_subgraph(xgr, branch_bond_keys(xgr, atm_key, bnd_key))


def branch_bond_keys(xgr, atm_key, bnd_key):
    """ keys for branch extending along `bnd_key` away from `atm_key`
    """
    bnd_key = frozenset(bnd_key)
    assert atm_key in bnd_key
    assert bnd_key in _bond_keys(xgr)

    atm_bnd_keys_dct = atom_bond_keys(xgr)

    bnch_bnd_keys = {bnd_key}
    seen_bnd_keys = set()
    excl_bnd_keys = atm_bnd_keys_dct[atm_key] - {bnd_key}

    new_bnd_keys = {bnd_key}

    bnd_ngb_keys_dct = bond_neighbor_keys(xgr)
    while new_bnd_keys:
        new_bnd_ngb_keys = set(
            _chain(*_values_by_key(bnd_ngb_keys_dct, new_bnd_keys)))
        bnch_bnd_keys.update(new_bnd_ngb_keys - excl_bnd_keys)
        seen_bnd_keys.update(new_bnd_keys)
        new_bnd_keys = bnch_bnd_keys - seen_bnd_keys

    return frozenset(bnch_bnd_keys)


def rings(xgr):
    """ rings in the graph (minimal basis)
    """
    xgrs = [bond_induced_subgraph(xgr, bnd_keys)
            for bnd_keys in rings_bond_keys(xgr)]
    return tuple(sorted(xgrs, key=_frozen))


def rings_bond_keys(xgr):
    """ bond keys for each ring in the graph (minimal basis)
    """
    bnd_keys = _bond_keys(xgr)

    def _ring_bond_keys(rng_atm_keys):
        return frozenset(filter(lambda x: x <= rng_atm_keys, bnd_keys))

    nxg = _nxg_from_graph(xgr)
    rng_atm_keys_lst = _nxg_minimal_rings_atom_keys(nxg)
    rng_bnd_keys_lst = frozenset(map(_ring_bond_keys, rng_atm_keys_lst))
    return rng_bnd_keys_lst


def subgraph(xgr, atm_keys):
    """ the subgraph induced by a subset of the atoms
    """
    atm_keys = set(atm_keys)
    assert atm_keys <= _atom_keys(xgr)
    bnd_keys = set(filter(lambda x: x <= atm_keys, _bond_keys(xgr)))
    atm_dct = _by_key(_atoms(xgr), atm_keys)
    bnd_dct = _by_key(_bonds(xgr), bnd_keys)
    return _from_atoms_and_bonds(atm_dct, bnd_dct)


def bond_induced_subgraph(xgr, bnd_keys):
    """ the subgraph induced by a subset of the bonds
    """
    atm_keys = set(_chain(*bnd_keys))
    bnd_keys = set(bnd_keys)
    assert atm_keys <= _atom_keys(xgr)
    assert bnd_keys <= _bond_keys(xgr)
    atm_dct = _by_key(_atoms(xgr), atm_keys)
    bnd_dct = _by_key(_bonds(xgr), bnd_keys)
    return _from_atoms_and_bonds(atm_dct, bnd_dct)


# transformations
def delete_atoms(xgr, atm_keys):
    """ delete atoms from the molecular graph
    """
    all_atm_keys = _atom_keys(xgr)
    atm_keys = set(atm_keys)
    assert atm_keys <= all_atm_keys
    atm_keys_left = all_atm_keys - atm_keys
    return subgraph(xgr, atm_keys_left)
