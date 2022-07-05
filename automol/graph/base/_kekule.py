""" graph functions associated with kekule (resonance) structures

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
from automol.util import dict_
from automol.graph.base._core import subgraph
from automol.graph.base._core import implicit
from automol.graph.base._core import atom_keys
from automol.graph.base._core import atom_hybridizations
from automol.graph.base._core import bond_keys
from automol.graph.base._core import bond_orders
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import set_bond_orders
from automol.graph.base._core import atom_unsaturations
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import dummy_atoms_neighbor_atom_key
from automol.graph.base._core import without_bond_orders
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import without_dummy_bonds
from automol.graph.base._core import without_fractional_bonds
from automol.graph.base._algo import connected_components
from automol.graph.base._algo import connected_components_atom_keys
# from automol.graph.base._algo import weighted_maximal_matching


# # core functions
def kekule(gra, max_stereo_overlap=True):
    """ One low-spin kekule graph, ignoring current bond orders

        Low-spin kekule graphs have double and triple bonds assigned to
        minimize the number of unpaired electrons.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param max_stereo_overlap: optionally, request as many stereo bonds as
            possible to have a bond order of 2
        :type max_stereo_overlap: bool
        :returns: a kekule graph
    """
    ste_bnd_keys = bond_stereo_keys(gra)

    def _count_stereo_double_bonds(gra):
        bnd_ord_dct = bond_orders(gra)
        count = sum(bnd_ord_dct[k] == 2 for k in ste_bnd_keys)
        return count

    gras = kekules(gra)
    if not max_stereo_overlap:
        gra = next(iter(gras))
    else:
        gra = max(gras, key=_count_stereo_double_bonds)

    return gra


def kekules(gra):
    """ All possible low-spin kekule graphs, ignoring current bond orders

        Low-spin kekule graphs have double and triple bonds assigned to
        minimize the number of unpaired electrons.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: all possible low-spin kekule graphs
    """
    orig_gra = gra
    bnd_ord_dcts = kekules_bond_orders(gra)
    gras = []
    for bnd_ord_dct in bnd_ord_dcts:
        gra = set_bond_orders(orig_gra, bnd_ord_dct)
        gras.append(gra)
    return tuple(gras)


def kekule_bond_orders(gra, max_stereo_overlap=True):
    """ Bond orders for one low-spin kekule graph

        Low-spin kekule graphs have double and triple bonds assigned to
        minimize the number of unpaired electrons.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param max_stereo_overlap: optionally, request as many stereo bonds as
            possible to have a bond order of 2
        :type max_stereo_overlap: bool
        :returns: a bond order dictionary
    """
    return bond_orders(kekule(gra, max_stereo_overlap=max_stereo_overlap))


def kekules_bond_orders(gra):
    """ Bond orders for all possible low-spin kekule graphs

        Low-spin kekule graphs have double and triple bonds assigned to
        minimize the number of unpaired electrons.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: bond orders for all possible low-spin kekule graphs
        :rtype: tuple[dict]
    """
    orig_bnd_ord_dct = bond_orders(gra)
    gra = implicit(without_dummy_bonds(without_fractional_bonds(gra)))
    gra = without_bond_orders(gra)

    # identify all of the independent pi systems and assign kekules to each
    pi_keys_lst = pi_system_atom_keys(gra)
    pi_bord_dcts_lst = [pi_system_kekules_bond_orders(gra, pi_keys)
                        for pi_keys in pi_keys_lst]

    bnd_ord_dcts = []
    # combine the kekules from each pi system together in all possible ways
    for bord_dcts in itertools.product(*pi_bord_dcts_lst):
        bord_dct = orig_bnd_ord_dct.copy()
        for dct in bord_dcts:
            bord_dct.update(dct)
        bnd_ord_dcts.append(bord_dct)

    if bnd_ord_dcts:
        bnd_ord_dcts = tuple(bnd_ord_dcts)
    else:
        bnd_ord_dcts = (orig_bnd_ord_dct,)

    return bnd_ord_dcts


def kekules_bond_orders_collated(gra):
    """ Bond orders for all possible low-spin kekule graphs, collated into a
        single dictionary

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: bond orders for all possible low-spin kekule graphs
        :rtype: tuple[dict]
    """
    bnd_keys = list(bond_keys(gra))
    bnd_ords_lst = list(dict_.values_by_key(d, bnd_keys) for d in
                        kekules_bond_orders(gra))
    bnd_ords_dct = dict(zip(bnd_keys, zip(*bnd_ords_lst)))
    return bnd_ords_dct


def kekules_bond_orders_averaged(gra):
    """ Bond orders for all possible low-spin kekule graphs, averaged

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: bond orders for all possible low-spin kekule graphs
        :rtype: tuple[dict]
    """
    bnd_ords_dct = kekules_bond_orders_collated(gra)
    avg_bnd_ord_dct = {k: sum(v)/len(v) for k, v in bnd_ords_dct.items()}
    return avg_bnd_ord_dct


# # derived properties
def linear_atom_keys(gra, dummy=True):
    """ atoms forming linear bonds, based on their hybridization

    :param gra: the graph
    :param dummy: whether or not to consider atoms connected to dummy atoms as
        linear, if different from what would be predicted based on their
        hybridization
    :returns: the linear atom keys
    :rtype: tuple[int]
    """
    gra = without_fractional_bonds(gra)
    atm_hyb_dct = atom_hybridizations(implicit(kekule(gra)))
    lin_atm_keys = set(dict_.keys_by_value(atm_hyb_dct, lambda x: x == 1))

    if dummy:
        dum_ngb_key_dct = dummy_atoms_neighbor_atom_key(gra)
        lin_atm_keys |= set(dum_ngb_key_dct.values())

    lin_atm_keys = tuple(sorted(lin_atm_keys))
    return lin_atm_keys


def linear_segments_atom_keys(gra, lin_keys=None):
    """ atom keys for linear segments in the graph
    """
    ngb_keys_dct = atoms_neighbor_atom_keys(without_dummy_atoms(gra))

    lin_keys = (linear_atom_keys(gra, dummy=True)
                if lin_keys is None else lin_keys)

    lin_keys = [k for k in lin_keys if len(ngb_keys_dct[k]) <= 2]

    lin_segs = connected_components(subgraph(gra, lin_keys))

    lin_keys_lst = []
    for lin_seg in lin_segs:
        lin_seg_keys = atom_keys(lin_seg)
        if len(lin_seg_keys) == 1:
            key, = lin_seg_keys
            lin_keys_lst.append([key])
        else:
            end_key1, end_key2 = sorted([
                key for key, ngb_keys in
                atoms_neighbor_atom_keys(lin_seg).items()
                if len(ngb_keys) == 1])
            ngb_keys_dct = atoms_neighbor_atom_keys(lin_seg)

            key = None
            keys = [end_key1]
            while key != end_key2:
                key, = ngb_keys_dct[keys[-1]] - set(keys)
                keys.append(key)
            lin_keys_lst.append(keys)

    lin_keys_lst = tuple(map(tuple, lin_keys_lst))
    return lin_keys_lst


# # helpers
def pi_system_atom_keys(gra):
    """ Extract keys for each closed, connected pi-system of a molecule

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: keys for each closed, connected pi system
    """
    atm_unsat_dct = atom_unsaturations(gra, bond_order=False)
    all_pi_keys = dict_.keys_by_value(atm_unsat_dct)
    pi_keys_lst = tuple(
        ks for ks in
        connected_components_atom_keys(subgraph(gra, all_pi_keys))
        if len(ks) > 1)
    return pi_keys_lst


def pi_system_kekules_bond_orders(gra, pi_keys, log=False):
    """ Determine kekules for a closed, connected pi-system

        A completely general algorithm.

        Strategy:
        Recursively visit one atom at a time and loop through all possible bond
        assignments for it with its neighbors, starting from the greatest
        number of bonds and working down to the least. Then continue the
        recursion for each neighboring atom.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param pi_keys: keys of an closed, connected pi system
        :type pi_keys: frozenset[int]
        :param log: print a debugging log of the number of recursive calls?
        :type log: bool
    """
    pi_sy = subgraph(gra, pi_keys)
    atm_keys = atom_keys(pi_sy)
    bnd_keys = bond_keys(pi_sy)

    nkeys_dct = atoms_neighbor_atom_keys(pi_sy)

    aus_dct = atom_unsaturations(gra, bond_order=False)
    spin_min = sum(aus_dct.values())

    bord_dcts = []

    calls = 0

    def _recurse_expand(key, vkeys_dct, bord_dct, aus_dct, spin):
        nonlocal bord_dcts
        nonlocal spin_min
        nonlocal calls
        calls += 1

        nset = nkeys_dct[key] - vkeys_dct[key]
        # Limit the neighbor pool to <=2 neighbors to prevent quadruple bonds
        # (only applies to [C][C])
        npool = list(itertools.chain(*(
            itertools.repeat(n, min(aus_dct[n], 2)) for n in nset)))

        spin0 = spin
        aus_dct0 = aus_dct
        bord_dct0 = bord_dct
        vkeys_dct = vkeys_dct.copy()

        max_inc = aus_dct[key]
        if nset:
            # Update the dictionary of visited keys
            vkeys_dct[key] |= nset
            for nkey in nset:
                vkeys_dct[nkey] |= {key}

            incs = list(reversed(range(max_inc+1)))
            for inc in incs:
                spin = spin0 + (max_inc - inc)

                if spin <= spin_min:
                    nkeys_lst = list(set(itertools.combinations(npool, inc)))
                    for nkeys in nkeys_lst:
                        aus_dct = aus_dct0.copy()
                        bord_dct = bord_dct0.copy()

                        for nkey in nkeys:
                            bkey = frozenset({key, nkey})
                            bord_dct[bkey] += 1
                            aus_dct[key] -= 1
                            aus_dct[nkey] -= 1

                        if vkeys_dct == nkeys_dct:
                            spin = sum(aus_dct.values())
                            if bord_dct not in bord_dcts:
                                if spin == spin_min:
                                    bord_dcts.append(bord_dct)
                                elif spin < spin_min:
                                    bord_dcts = [bord_dct]
                                    spin_min = spin

                            continue

                        for nkey in nset:
                            _recurse_expand(
                                nkey, vkeys_dct, bord_dct, aus_dct, spin)

    start_key = sorted(atom_keys(pi_sy))[0]
    vkeys_dct = {k: frozenset() for k in atm_keys}
    bord_dct = {k: 1 for k in bnd_keys}
    start_spin = 0
    _recurse_expand(start_key, vkeys_dct, bord_dct, aus_dct, start_spin)

    if log:
        print(f"calls {calls}")

    return bord_dcts


# def pi_system_kekule_bond_orders(gra, pi_keys, bnd_key=None):
#     """ Determine a kekule structure for a closed, connected pi-system
#
#         Optionally, try making a particular bond key have a double bond.
#
#         :param gra: molecular graph
#         :type gra: automol graph data structure
#         :param pi_keys: keys of an closed, connected pi system
#         :type pi_keys: frozenset[int]
#         :param bnd_key: optionally, make this bond have a double bond
#         :type bnd_key: frozenset[int]
#     """
#     pi_sy = subgraph(gra, pi_keys)
#     bnd_weight_dct = (None if bnd_key is None else
#                       {frozenset(bnd_key): 2})
#
#     bnd_keys = weighted_maximal_matching(
#         pi_sy, bnd_weight_dct=bnd_weight_dct)
#     print(bnd_keys)


if __name__ == '__main__':
    # C=C[CH2]
    GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 2, None)},
           {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})

    # C=CC=C[CH2]
    GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 1, None),
            3: ('C', 1, None), 4: ('C', 2, None)},
           {frozenset({3, 4}): (1, None), frozenset({0, 1}): (1, None),
            frozenset({2, 3}): (1, None), frozenset({1, 2}): (1, None)})

    # # C=C=C=C
    # GRA = ({0: ('C', 2, None), 1: ('C', 0, None), 2: ('C', 0, None),
    #         3: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None)})

    # # C=CC=CC=C
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 2, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None)})

    # # C=C-C(-[CH2])=C
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 2, None), 6: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({2, 6}): (1, None)})

    # # C1=CC=CC=C1 (benzene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({0, 5}): (1, None)})

    # # C12=CC=C1C=C2  _  _
    # #              ||_||_||
    # GRA = ({0: ('C', 0, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({0, 3}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({0, 5}): (1, None)})

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

    # # C1=CC=C2C=C3C=CC=CC3=CC2=C1 (anthracene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 0, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 1, None),
    #         9: ('C', 1, None), 10: ('C', 0, None), 11: ('C', 1, None),
    #         12: ('C', 0, None), 13: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({0, 13}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({11, 12}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
    #         frozenset({8, 7}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({3, 12}): (1, None), frozenset({10, 5}): (1, None),
    #         frozenset({12, 13}): (1, None), frozenset({10, 11}): (1, None)})

    # # C1=CC=C2C(=C1)C=CC3=CC=CC=C32 (phenanthrene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None), 10: ('C', 1, None), 11: ('C', 1, None),
    #         12: ('C', 1, None), 13: ('C', 0, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({4, 6}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({11, 12}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
    #         frozenset({3, 13}): (1, None), frozenset({0, 5}): (1, None),
    #         frozenset({8, 7}): (1, None), frozenset({8, 13}): (1, None),
    #         frozenset({12, 13}): (1, None), frozenset({10, 11}): (1, None)})

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

    # # [CH2]C=CCC=CC=C[CH2]
    # GRA = ({0: ('C', 2, None), 3: ('C', 1, None), 4: ('C', 1, None),
    #         5: ('C', 2, None), 6: ('C', 1, None), 7: ('C', 1, None),
    #         8: ('C', 1, None), 9: ('C', 1, None), 10: ('C', 2, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({9, 10}): (1, None),
    #         frozenset({0, 3}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
    #         frozenset({8, 7}): (1, None), frozenset({5, 6}): (1, None)})

    # [C][C]
    GRA = ({0: ('C', 0, None), 1: ('C', 0, None)},
           {frozenset({0, 1}): (1, None)})

    # C=C[CH2]
    GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 2, None)},
           {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})

    # C1=CC=CC=C1 (benzene)
    GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
            3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None)},
           {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
            frozenset({0, 1}): (1, None), frozenset({0, 5}): (1, None)})

    print(len(kekules(GRA)))
    print(kekules_bond_orders_collated(GRA))
