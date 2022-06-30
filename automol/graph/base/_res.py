""" graph functions that depend on resonance structure

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
from automol.util import dict_
from automol.graph.base._core import subgraph
from automol.graph.base._core import implicit
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bond_orders
from automol.graph.base._core import atom_unsaturations
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

    if len(bord_dcts) <= 5:
        for bord_dct in bord_dcts:
            print(bord_dct)
            print()


def pi_system_resonance_bond_orders(gra, pi_keys, atm_unsat_dct=None):
    """ Determine resonances for a closed, connected pi-system

    Strategy:
    Sweep through the whole network from neighbor to neighbor, initially
    saturating maximally, then gradually unsaturating in all possible ways.
    Once all neighbors of an atom have been assigned, its contribution to the
    multiplicity is fixed -- add this to the multiplicity lower-bound.
    Add up the multiplicities and quit a path as soon as the multiplicity lb
    exceeds the current minimum multiplicity.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param pi_keys: keys of an closed, connected pi system
        :type pi_keys: frozenset[int]
        :param atm_unsat_dct: atom unsaturations, by atom key
        :type atm_unsat_dct: dict[int: int]
    """
    atm_unsat_dct = (atom_unsaturations(gra, bond_order=False)
                     if atm_unsat_dct is None else atm_unsat_dct)

    pi_sy = subgraph(gra, pi_keys)
    atm_keys = atom_keys(pi_sy)
    bnd_ord_dct = bond_orders(pi_sy)

    nkeys_dct = atoms_neighbor_atom_keys(pi_sy)

    bord_dcts = []
    spin_min = sum(atm_unsat_dct.values())

    call = 0

    def _recurse_expand(key, seen_keys, bord_dct, aus_dct, spin):
        nonlocal spin_min
        nonlocal bord_dcts
        nonlocal call
        call += 1
        print()
        print(f"call {call}")
        print(f"key {key}")
        print(f"seen_keys {seen_keys}")
        print(f"bord_dct {bord_dct}")
        print(f"aus_dct {aus_dct}")
        print(f"spin {spin}")

        nset = nkeys_dct[key] - seen_keys
        npool = list(itertools.chain(*(
            itertools.repeat(n, aus_dct[n]) for n in nset)))

        spin0 = spin
        aus_dct0 = aus_dct
        bord_dct0 = bord_dct

        seen_keys |= nset

        max_inc = aus_dct[key]
        if npool:
            incs = list(reversed(range(max_inc+1)))
            print(f'A{call} incs {incs}')
            for inc in incs:
                print(f'A{call} inc {inc}')
                spin = spin0 + (max_inc - inc)
                print(f'A{call} spin {spin}')
                print(f'A{call} spin_min {spin_min}')

                # Use set to avoid duplicates
                if spin <= spin_min:
                    nkeys_lst = list(set(itertools.combinations(npool, inc)))
                    print(f'A{call} nkeys_lst {nkeys_lst}')
                    for nkeys in nkeys_lst:
                        aus_dct = aus_dct0.copy()
                        bord_dct = bord_dct0.copy()
                        print(f"A{call}")
                        print(f"A{call} nkeys {nkeys}")
                        for nkey in nkeys:
                            bkey = frozenset({key, nkey})
                            bord_dct[bkey] += 1
                            aus_dct[key] -= 1
                            aus_dct[nkey] -= 1

                        print(f"A{call} aus_dct {aus_dct}")
                        print(f"A{call} bord_dct {bord_dct}")
                        print(f"A{call} seen_keys {seen_keys}")

                        if seen_keys == atm_keys:
                            spin = sum(aus_dct.values())
                            if bord_dct not in bord_dcts:
                                print(f"A{call} - New structure added!")
                                if spin == spin_min:
                                    bord_dcts.append(bord_dct)
                                elif spin < spin_min:
                                    bord_dcts = [bord_dct]
                                    spin_min = spin
                                    print(f"A{call} - Spin lowered to {spin}!")

                            continue

                        print(f"A{call} before recursion")

                        for nkey in nset:
                            _recurse_expand(nkey, seen_keys, bord_dct, aus_dct,
                                            spin)

                        print(f"A{call} after recursion")

    start_key = sorted(atom_keys(pi_sy))[0]
    seen_keys = frozenset({start_key})
    bord_dct = bnd_ord_dct
    aus_dct = atm_unsat_dct
    start_spin = 0
    _recurse_expand(start_key, seen_keys, bord_dct, aus_dct, start_spin)

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
    # C=C[CH2]
    GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 2, None)},
           {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})

    # C=C=C
    GRA = ({0: ('C', 2, None), 1: ('C', 0, None), 2: ('C', 2, None)},
           {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})

    # C=C=C=C
    GRA = ({0: ('C', 2, None), 1: ('C', 0, None), 2: ('C', 0, None),
            3: ('C', 2, None)},
           {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({1, 2}): (1, None)})

    # C=CC=CC=C
    GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 1, None),
            3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 2, None)},
           {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
            frozenset({0, 1}): (1, None)})

    # C=C-C(-[CH2])=C
    GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 0, None),
            3: ('C', 2, None), 6: ('C', 2, None)},
           {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({1, 2}): (1, None), frozenset({2, 6}): (1, None)})

    # C1=CC=CC=C1 (benzene)
    GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
            3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None)},
           {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
            frozenset({0, 1}): (1, None), frozenset({0, 5}): (1, None)})

    # # C1=CC=C2C=CC=CC2=C1 (naphthalene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({6, 7}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({8, 3}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({0, 9}): (1, None)})

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
