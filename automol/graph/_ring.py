""" ring graph library
"""
import operator
import functools
from automol.graph._graph import frozen as _frozen
from automol.graph._graph import bond_keys as _bond_keys
from automol.graph._graph import (bond_induced_subgraph as
                                  _bond_induced_subgraph)
from automol.graph import _networkx


def rings(gra):
    """ rings in the graph (minimal basis)
    """
    gras = [_bond_induced_subgraph(gra, bnd_keys)
            for bnd_keys in rings_bond_keys(gra)]
    return tuple(sorted(gras, key=_frozen))


def rings_atom_keys(gra):
    """ atom keys for each ring in the graph sorted by connectivity (minimal basis)
    """
    def _sorted_ring_atom_keys(rng_bnd_keys):
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

    rng_atm_keys_lst = frozenset(
        map(_sorted_ring_atom_keys, rings_bond_keys(gra)))
    return rng_atm_keys_lst


def rings_bond_keys(gra):
    """ bond keys for each ring in the graph (minimal basis)
    """
    bnd_keys = _bond_keys(gra)

    def _ring_bond_keys(rng_atm_keys):
        return frozenset(filter(lambda x: x <= rng_atm_keys, bnd_keys))

    nxg = _networkx.from_graph(gra)
    rng_atm_keys_lst = _networkx.minimum_cycle_basis(nxg)
    rng_bnd_keys_lst = frozenset(map(_ring_bond_keys, rng_atm_keys_lst))
    return rng_bnd_keys_lst


def ring_systems(gra):
    """ polycyclic ring systems in the graph
    """
    gras = [_bond_induced_subgraph(gra, bnd_keys)
            for bnd_keys in ring_systems_bond_keys(gra)]
    return tuple(sorted(gras, key=_frozen))


def ring_systems_bond_keys(gra):
    """ polycyclic ring systems in the graph
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
    rsy_bnd_keys_lsts = _equivalence_partition(rng_bnd_keys_lst,
                                               _are_connected)
    rsy_bnd_keys_lst = [frozenset(functools.reduce(operator.or_, bnd_keys_lst))
                        for bnd_keys_lst in rsy_bnd_keys_lsts]
    return rsy_bnd_keys_lst


def _equivalence_partition(iterable, relation):
    """Partitions a set of objects into equivalence classes

    canned function taken from https://stackoverflow.com/a/38924631

    Args:
        iterable: collection of objects to be partitioned
        relation: equivalence relation. I.e. relation(o1,o2) evaluates to True
            if and only if o1 and o2 are equivalent

    Returns: classes, partitions
        classes: A sequence of sets. Each one is an equivalence class
    """
    classes = []
    for obj in iterable:  # for each object
        # find the class it is in
        found = False
        for cls in classes:
            # is it equivalent to this class?
            if relation(next(iter(cls)), obj):
                cls.add(obj)
                found = True
                break
        if not found:  # it is in a new class
            classes.append(set([obj]))
    return classes
