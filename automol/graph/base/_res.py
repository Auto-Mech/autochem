""" graph functions that depend on resonance structure

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import numpy
from automol.util import dict_
from automol.graph.base._core import subgraph
from automol.graph.base._core import implicit
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bond_keys
from automol.graph.base._core import bond_orders
from automol.graph.base._core import atom_unsaturations
from automol.graph.base._core import bond_unsaturations
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._algo import connected_components_atom_keys


def function(gra):
    """ testing
    """
    gra = implicit(gra)
    atm_unsat_dct = atom_unsaturations(gra, bond_order=False)
    pi_keys_lst = pi_system_atom_keys(gra, atm_unsat_dct=atm_unsat_dct)

    pi_keys = list(pi_keys_lst)[-1]
    bord_dcts = pi_system_resonance_bond_orders(gra, pi_keys,
                                                atm_unsat_dct=atm_unsat_dct)
    print("FINAL")
    print(len(bord_dcts))


def pi_system_resonance_bond_orders(gra, pi_keys, atm_unsat_dct=None):
    """ Determine resonances for a closed, connected pi-system

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param pi_keys: keys of an closed, connected pi system
        :type pi_keys: frozenset[int]
        :param atm_unsat_dct: atom unsaturations, by atom key
        :type atm_unsat_dct: dict[int: int]
    """
    atm_unsat_dct = (atom_unsaturations(gra, bond_order=False)
                     if atm_unsat_dct is None else atm_unsat_dct)
    # bnd_unsat_dct = bond_unsaturations(gra, bond_order=False)

    pi_sy = subgraph(gra, pi_keys)
    atm_keys = sorted(atom_keys(pi_sy))
    bnd_ord_dct = bond_orders(pi_sy)

    nrn_dct = {k: sorted(n for n in ns if k < n)
               for k, ns in atoms_neighbor_atom_keys(pi_sy).items()}

    # bord_dcts = []

    natms = len(atm_keys)
    last_idx = len(atm_keys) - 1

    paths = []

    def _recurse_paths(idx, path):
        print(f"path {path}")
        print(f"idx {idx}")
        key = atm_keys[idx]
        print(f"key {key}")
        if key not in path:
            path.append(key)
            print(f"A path {path}")

            parent_path = path
            nkeys = [n for n in nrn_dct[key] if n not in parent_path]
            print(f"A nkeys {nkeys}")
            for nkey in nkeys:
                # Don't copy unless we have to
                path = (parent_path.copy() if nkey != nkeys[-1] else
                        parent_path)
                path.append(nkey)
                print(f"A new path {path}")
                print(f"A recursion with {idx+1}")
                path = _recurse_paths(idx+1, path)

        if len(path) < natms:
            print(f"B recursion with {idx+1}")
            path = _recurse_paths(idx+1, path)
        elif path not in paths:
            print(f"Appending to paths: {path} (idx {idx})")
            paths.append(path)

        return path

    _recurse_paths(0, [])

    return paths


def old_pi_system_resonance_bond_orders(gra, pi_keys, atm_unsat_dct=None):
    """ Determine resonances for a closed, connected pi-system

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param pi_keys: keys of an closed, connected pi system
        :type pi_keys: frozenset[int]
        :param atm_unsat_dct: atom unsaturations, by atom key
        :type atm_unsat_dct: dict[int: int]
    """
    atm_unsat_dct = (atom_unsaturations(gra, bond_order=False)
                     if atm_unsat_dct is None else atm_unsat_dct)
    bnd_unsat_dct = bond_unsaturations(gra, bond_order=False)

    pi_sy = subgraph(gra, pi_keys)
    bnd_keys = bond_keys(pi_sy)
    seed_pool = {(k, 1+i) for k in bnd_keys
                 for i in range(0, bnd_unsat_dct[k]+1)}
    print(seed_pool)
    print(f"len pool: {len(seed_pool)}")

    nkeys_dct = atoms_neighbor_atom_keys(pi_sy)

    def _recurse_saturate(bord_dct, aus_dct, key, seen_keys):
        seen_keys.add(key)
        nkeys = nkeys_dct[key] - seen_keys
        for nkey in nkeys:
            bkey = frozenset({key, nkey})
            inc = min(aus_dct[key], aus_dct[nkey])
            if inc:
                bord_dct[bkey] += inc
                aus_dct[key] -= inc
                aus_dct[nkey] -= inc

                if (bkey, 1+inc) in seed_pool:
                    seed_pool.remove((bkey, 1+inc))
            bord_dct, aus_dct = _recurse_saturate(bord_dct, aus_dct, nkey,
                                                  seen_keys=seen_keys)
        return bord_dct, aus_dct

    bord_dcts = []
    min_mult = numpy.inf
    while seed_pool:
        bord_dct = {k: 1 for k in bnd_keys}

        bkey, bord = seed_pool.pop()
        print(bkey, bord)
        bord_dct[bkey] = bord

        key = sorted(bkey)[0]
        aus_dct = dict_.by_key(atm_unsat_dct, pi_keys)
        seen_keys = set()
        bord_dct, aus_dct = _recurse_saturate(bord_dct, aus_dct, key,
                                              seen_keys=seen_keys)
        mult = sum(aus_dct.values()) + 1
        print(f"mult {mult}")
        if mult < min_mult:
            print("Found lower multiplicity!")
            bord_dcts = []
            min_mult = mult

        if bord_dct not in bord_dcts and mult == min_mult:
            bord_dcts.append(bord_dct)

    return bord_dcts


def pi_system_atom_keys(gra, atm_unsat_dct=None):
    """ Extract keys for each closed, connected pi-system of a molecule

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param atm_unsat_dct: atom unsaturations, by atom key
        :type atm_unsat_dct: dict[int: int]
        :returns: keys for each closed, connected pi system
    """
    atm_unsat_dct = (atom_unsaturations(gra, bond_order=False)
                     if atm_unsat_dct is None else atm_unsat_dct)
    all_pi_keys = dict_.keys_by_value(atm_unsat_dct)
    pi_keys_lst = connected_components_atom_keys(subgraph(gra, all_pi_keys))
    return pi_keys_lst


if __name__ == '__main__':
    # # C=C[CH2]
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})

    # # C=CC=CC=C
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 2, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None)})

    # C=C-C(-[CH2])=C
    GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 0, None),
            3: ('C', 2, None), 6: ('C', 2, None)},
           {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({1, 2}): (1, None), frozenset({2, 6}): (1, None)})

    # # C1=CC=CC=C1
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({0, 5}): (1, None)})

    # # C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2 (pyrene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None), 10: ('C', 1, None), 11: ('C', 1, None),
    #         12: ('C', 0, None), 13: ('C', 0, None), 14: ('C', 1, None),
    #         15: ('C', 1, None)},
    #        {frozenset({4, 6}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({11, 12}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 5}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({14, 15}): (1, None), frozenset({10, 11}): (1, None),
    #         frozenset({3, 4}): (1, None), frozenset({3, 13}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({6, 7}): (1, None),
    #         frozenset({8, 13}): (1, None), frozenset({12, 14}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({2, 15}): (1, None),
    #         frozenset({12, 13}): (1, None)})

    # # C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7 (coronene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 0, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 0, None), 10: ('C', 0, None), 11: ('C', 0, None),
    #         12: ('C', 1, None), 13: ('C', 1, None), 14: ('C', 1, None),
    #         15: ('C', 1, None), 16: ('C', 0, None), 17: ('C', 0, None),
    #         18: ('C', 0, None), 19: ('C', 0, None), 20: ('C', 1, None),
    #         21: ('C', 1, None), 22: ('C', 1, None), 23: ('C', 1, None)},
    #        {frozenset({17, 10}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({11, 12}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 5}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({9, 4}): (1, None), frozenset({14, 15}): (1, None),
    #         frozenset({10, 11}): (1, None), frozenset({3, 4}): (1, None),
    #         frozenset({16, 17}): (1, None), frozenset({22, 23}): (1, None),
    #         frozenset({16, 15}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({18, 19}): (1, None), frozenset({18, 3}): (1, None),
    #         frozenset({11, 14}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({2, 21}): (1, None),
    #         frozenset({17, 18}): (1, None), frozenset({8, 13}): (1, None),
    #         frozenset({20, 21}): (1, None), frozenset({19, 20}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({16, 23}): (1, None),
    #         frozenset({19, 22}): (1, None), frozenset({12, 13}): (1, None)})

    function(GRA)
