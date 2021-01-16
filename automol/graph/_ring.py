""" ring graph library
"""
import operator
import itertools
import functools
import more_itertools as mit
from automol import util
from automol.graph._graph_base import string
from automol.graph._graph import frozen
from automol.graph._graph import atom_count
from automol.graph._graph import atom_keys
from automol.graph._graph import bond_keys
from automol.graph._graph import remove_bonds
from automol.graph._graph import is_connected
from automol.graph._graph import atom_bond_keys
from automol.graph._graph import atom_shortest_paths
from automol.graph._graph import union_from_sequence
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
        map(sorted_ring_atom_keys_from_bond_keys, rings_bond_keys(gra)))
    return rng_atm_keys_lst


def sorted_ring_atom_keys(rng):
    """ get a ring's atom keys, sorted in order of connectivity
    """
    return sorted_ring_atom_keys_from_bond_keys(bond_keys(rng))


def sorted_ring_atom_keys_from_bond_keys(rng_bnd_keys):
    """ get a ring's atom keys, sorted in order of connectivity, from its bond
    keys
    """
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


def cycle_ring_atom_key_to_front(keys, key, end_key=None):
    """ helper function to cycle ring atom keys until one is in front

    :param keys: ring keys
    :parm key: the key to cycle to the font
    :param end_key: optionally, ensure that another key is the last key in the
        ring; note that this is only possible if key and end_key are adjacent
    """
    assert key in keys, ("{:d} is not in {:s}".format(key, str(keys)))
    keys = tuple(itertools.islice(
        itertools.dropwhile(lambda x: x != key, itertools.cycle(keys)),
        len(keys)))

    if end_key is not None and keys[-1] != end_key:
        assert keys[1] == end_key, (
            "end_key {:d} is not adjacent to {:d} in the ring"
            .format(key, end_key))
        keys = list(reversed(keys))
        keys = cycle_ring_atom_key_to_front(keys, key)

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
    atm_keys = sorted_ring_atom_keys_from_bond_keys(bond_keys(rng))

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
    rsy_bnd_keys_lsts = util.equivalence_partition(rng_bnd_keys_lst,
                                                   _are_connected)
    rsy_bnd_keys_lst = [frozenset(functools.reduce(operator.or_, bnd_keys_lst))
                        for bnd_keys_lst in rsy_bnd_keys_lsts]
    return rsy_bnd_keys_lst


def is_ring_system(gra):
    """ is this graph a ring system?
    """
    return union_from_sequence(rings(gra), check=False) == gra


def ring_system_decomposed_atom_keys(rsy, rng_keys=None, check=True):
    """ decomposed atom keys for a polycyclic ring system in a graph

    The ring system is decomposed into a ring and a series of arcs that can
    be used to successively construct the system

    :param rsy: the ring system
    :param rng_keys: keys for the first ring in the decomposition; if None, the
        smallest ring in the system will be chosen
    """
    if rng_keys is None:
        rng = sorted(rings(rsy), key=atom_count)[0]
        rng_keys = sorted_ring_atom_keys(rng)

    # check the arguments, if requested
    if check:
        # check that the graph is connected
        assert is_connected(rsy), "Ring system can't be disconnected."

        # check that the graph is actually a ring system
        assert is_ring_system(rsy), (
            "This is not a ring system graph:\n{:s}".format(string(rsy)))

        # check that rng is a subgraph of rsy
        assert set(rng_keys) <= atom_keys(rsy), (
            "{}\n^ Rings system doesn't contain ring as subgraph:\n{}"
            .format(string(rsy, one_indexed=False), str(rng_keys)))

    bnd_keys = list(mit.windowed(rng_keys + rng_keys[:1], 2))

    # Remove bonds for the ring
    rsy = remove_bonds(rsy, bnd_keys)
    keys_lst = [rng_keys]
    done_keys = set(rng_keys)

    while bond_keys(rsy):

        # Determine shortest paths for the graph with one more ring/arc deleted
        sp_dct = atom_shortest_paths(rsy)

        # The shortest path will be the next shortest arc in the system
        arc_keys = min(
            (sp_dct[i][j] for i, j in itertools.combinations(done_keys, 2)
             if j in sp_dct[i]), key=len)

        # Add this arc to the list
        keys_lst.append(arc_keys)

        # Add these keys to the list of done keys
        done_keys |= set(arc_keys)

        # Delete tbond keys for the new arc and continue to the next iteration
        bnd_keys = list(map(frozenset, mit.windowed(arc_keys, 2)))
        rsy = remove_bonds(rsy, bnd_keys)

    keys_lst = tuple(map(tuple, keys_lst))
    return keys_lst


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
    atm_keys = sorted_ring_atom_keys_from_bond_keys(bnd_keys)

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


if __name__ == '__main__':
    import automol

    ICH = automol.smiles.inchi('C1CCC2CC(CCC3C4CCC5CC4C53)CC2C1')
    GRA = automol.inchi.graph(ICH)
    print(automol.graph.rings_atom_keys(GRA))
