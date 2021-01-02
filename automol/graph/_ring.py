""" ring graph library
"""
import operator
import itertools
import functools
import more_itertools as mit
from automol.graph._graph import frozen
from automol.graph._graph import atom_count
from automol.graph._graph import atom_keys
from automol.graph._graph import bond_keys
from automol.graph._graph import atom_bond_keys
from automol.graph._graph import bond_induced_subgraph
from automol.graph import _networkx


def rings(gra):
    """ rings in the graph (minimal basis)
    """
    gras = [bond_induced_subgraph(gra, bnd_keys)
            for bnd_keys in rings_bond_keys(gra)]
    return tuple(sorted(gras, key=frozen))


def rings_atom_keys(gra):
    """ atom keys for each ring in the graph sorted by connectivity (minimal basis)
    """
    rng_atm_keys_lst = frozenset(
        map(_sorted_ring_atom_keys, rings_bond_keys(gra)))
    return rng_atm_keys_lst


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


def rings_bond_keys(gra):
    """ bond keys for each ring in the graph (minimal basis)
    """
    bnd_keys = bond_keys(gra)

    def _ring_bond_keys(rng_atm_keys):
        return frozenset(filter(lambda x: x <= rng_atm_keys, bnd_keys))

    nxg = _networkx.from_graph(gra)
    rng_atm_keys_lst = _networkx.minimum_cycle_basis(nxg)
    rng_bnd_keys_lst = frozenset(map(_ring_bond_keys, rng_atm_keys_lst))
    return rng_bnd_keys_lst


def is_ring_key_sequence(gra, keys):
    """ does this sequence of keys share a ring?
    """
    keys = set(keys)
    return any(keys <= rng_keys for rng_keys in rings_atom_keys(gra))


def cycle_ring_atom_key_to_front(keys, key):
    """ helper function to cycle ring atom keys until one is in front
    """
    assert key in keys, ("{:d} is not in {:s}".format(key, str(keys)))
    keys = tuple(itertools.islice(
        itertools.dropwhile(lambda x: x != key, itertools.cycle(keys)),
        len(keys)))
    return keys


def ring_arc_complement_atom_keys(gra, rng):
    """ non-intersecting arcs from a ring that shares segments with a graph
    """
    gra_atm_bnd_dct = atom_bond_keys(gra)
    rng_atm_bnd_dct = atom_bond_keys(rng)

    # 1. find divergence points, given by the atom at which the divergence
    # occurs and the bond followed by the ring as it diverges
    div_dct = {}

    for atm_key in atom_keys(gra) & atom_keys(rng):
        div = rng_atm_bnd_dct[atm_key] - gra_atm_bnd_dct[atm_key]
        if div:
            bnd_key, = div
            div_dct[atm_key] = bnd_key

    # 2. cycle through the ring atoms; if you meet a starting divergence, start
    # an arc; extend the arc until you meet an ending divergence; repeat until
    # all divergences are accounted for
    atm_keys = _sorted_ring_atom_keys(bond_keys(rng))

    arcs = []
    arc = []
    for atm_key, next_atm_key in mit.windowed(itertools.cycle(atm_keys), 2):
        bnd_key = frozenset({atm_key, next_atm_key})

        # if we haven't started an arc, see if we are at a starting divergence;
        # if so, start the arc now and cross the divergence from our list
        if not arc:
            if atm_key in div_dct and div_dct[atm_key] == bnd_key:
                div_dct.pop(atm_key)

                arc.append(atm_key)
        # if we've started an arc, extend it; then, check if we are at an
        # ending divergence; if so, end the arc and cross the divergence from
        # our list; add it to our list of arcs
        else:
            arc.append(atm_key)

            if next_atm_key in div_dct and div_dct[next_atm_key] == bnd_key:
                div_dct.pop(next_atm_key)

                arc.append(next_atm_key)
                arcs.append(arc)
                arc = []

        # if no divergences are left, break out of the loop
        if not div_dct:
            break

    arcs = tuple(map(tuple, arcs))
    return arcs


def ring_systems(gra):
    """ polycyclic ring systems in the graph
    """
    gras = [bond_induced_subgraph(gra, bnd_keys)
            for bnd_keys in ring_systems_bond_keys(gra)]
    return tuple(sorted(gras, key=frozen))


def ring_systems_atom_keys(gra):
    """ bond keys for polycyclic ring systems in the graph
    """
    atm_keys_lst = tuple(map(atom_keys, ring_systems(gra)))
    return atm_keys_lst


def ring_systems_bond_keys(gra):
    """ bond keys for polycyclic ring systems in the graph
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


def ring_systems_decomposed_atom_keys(gra):
    """ decomposed atom keys for polycyclic ring systems in the graph

    each ring system is decomposed into a ring and a series of arcs that can be
    used to successively construct the system
    """
    rsys = ring_systems(gra)
    decomps = tuple(map(_decompose_ring_system_atom_keys, rsys))
    return decomps


def _decompose_ring_system_atom_keys(rsy):
    """ decompose a ring system into a ring and a series of arcs
    """
    # sort from smallest to largest
    rngs_pool = sorted(
        rings(rsy), key=lambda x: atom_count(x, with_implicit=False))

    decomp = ()
    decomp_bnd_keys = set({})

    rng = rngs_pool.pop(0)
    bnd_keys = bond_keys(rng)
    atm_keys = _sorted_ring_atom_keys(bnd_keys)

    decomp += (atm_keys,)
    decomp_bnd_keys.update(bnd_keys)

    while rngs_pool:
        decomp_rsy = bond_induced_subgraph(rsy, decomp_bnd_keys)
        for idx, rng in enumerate(rngs_pool):
            arcs = ring_arc_complement_atom_keys(decomp_rsy, rng)
            if arcs:
                rngs_pool.pop(idx)
                decomp += arcs
                decomp_bnd_keys.update(bond_keys(rng))

    return decomp


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


if __name__ == '__main__':
    import automol

    ICH = automol.smiles.inchi('C1CCC2CC(CCC3C4CCC5CC4C53)CC2C1')
    GRA = automol.inchi.graph(ICH)
    print(automol.graph.rings_atom_keys(GRA))
