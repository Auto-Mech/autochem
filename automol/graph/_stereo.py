""" stereo graph library
"""
import itertools
import functools
import numpy
from automol.util import dict_
import automol.util
from automol.graph._embed_dep import rings_bond_keys
from automol.graph._graph_dep import atom_stereo_parities
from automol.graph._graph_dep import bond_stereo_parities
from automol.graph._graph_dep import set_atom_stereo_parities
from automol.graph._graph_dep import set_bond_stereo_parities
from automol.graph._graph_dep import without_bond_orders
from automol.graph._graph_dep import without_stereo_parities
from automol.graph._graph_dep import atom_bond_valences
from automol.graph._graph_dep import atoms_neighbor_atom_keys
from automol.graph._graph_dep import string
from automol.graph._graph_dep import has_stereo
from automol.graph._graph_dep import atom_stereo_keys
from automol.graph._graph_dep import bond_stereo_keys
from automol.graph._graph_dep import stereo_priority_vector
from automol.graph._graph_dep import explicit
from automol.graph._graph_dep import sp2_bond_keys
from automol.graph._graph import frozen


def stereogenic_atom_keys(gra, assigned=False):
    """ Find stereogenic atoms in this graph.

    If the `assigned` flag is set to `False`, only  unassigned stereogenic
    atoms will be detected.

    :param gra: the graph
    :param assigned: Include atoms that already have stereo assignments?
    :param assigned: bool
    :returns: the stereogenic atom keys
    :rtype: frozenset
    """
    gra = without_bond_orders(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in
    atm_keys = dict_.keys_by_value(atom_bond_valences(gra), lambda x: x == 4)
    if not assigned:
        # Remove assigned stereo keys
        atm_keys -= atom_stereo_keys(gra)

    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    def _is_stereogenic(atm_key):
        atm_ngb_keys = list(atm_ngb_keys_dct[atm_key])
        pri_vecs = [stereo_priority_vector(gra, atm_key, atm_ngb_key)
                    for atm_ngb_key in atm_ngb_keys]
        ret = not any(pv1 == pv2
                      for pv1, pv2 in itertools.combinations(pri_vecs, r=2))
        return ret

    ste_gen_atm_keys = frozenset(filter(_is_stereogenic, atm_keys))
    return ste_gen_atm_keys


def stereogenic_bond_keys(gra, assigned=False):
    """ Find stereogenic bonds in this graph.

    If the `assigned` flag is set to `False`, only  unassigned stereogenic
    bonds will be detected.

    :param gra: the graph
    :param assigned: Include bonds that already have stereo assignments?
    :param assigned: bool
    :returns: the stereogenic bond keys
    :rtype: frozenset
    """
    gra = without_bond_orders(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in

    # get candidates: planar bonds
    bnd_keys = sp2_bond_keys(gra)

    if not assigned:
        # remove bonds that already have stereo assignments
        bnd_keys -= bond_stereo_keys(gra)

    bnd_keys -= functools.reduce(  # remove double bonds in small rings
        frozenset.union,
        filter(lambda x: len(x) < 8, rings_bond_keys(gra)), frozenset())

    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    def _is_stereogenic(bnd_key):
        atm1_key, atm2_key = bnd_key

        def _is_symmetric_on_bond(atm_key, atm_ngb_key):
            atm_ngb_keys = list(atm_ngb_keys_dct[atm_key] - {atm_ngb_key})

            if not atm_ngb_keys:                # C=:O:
                ret = True
            elif len(atm_ngb_keys) == 1:        # C=N:-X
                ret = False
            else:
                assert len(atm_ngb_keys) == 2   # C=C(-X)-Y
                ret = (stereo_priority_vector(gra, atm_key, atm_ngb_keys[0]) ==
                       stereo_priority_vector(gra, atm_key, atm_ngb_keys[1]))

            return ret

        return not (_is_symmetric_on_bond(atm1_key, atm2_key) or
                    _is_symmetric_on_bond(atm2_key, atm1_key))

    ste_gen_bnd_keys = frozenset(filter(_is_stereogenic, bnd_keys))
    return ste_gen_bnd_keys


def stereomers(gra):
    """ all stereomers, ignoring this graph's assignments
    """
    bool_vals = (False, True)

    def _expand_atom_stereo(sgr):
        atm_ste_keys = stereogenic_atom_keys(sgr)
        nste_atms = len(atm_ste_keys)
        sgrs = [set_atom_stereo_parities(sgr, dict(zip(atm_ste_keys,
                                                       atm_ste_par_vals)))
                for atm_ste_par_vals
                in itertools.product(bool_vals, repeat=nste_atms)]
        return sgrs

    def _expand_bond_stereo(sgr):
        bnd_ste_keys = stereogenic_bond_keys(sgr)
        nste_bnds = len(bnd_ste_keys)
        sgrs = [set_bond_stereo_parities(sgr, dict(zip(bnd_ste_keys,
                                                       bnd_ste_par_vals)))
                for bnd_ste_par_vals
                in itertools.product(bool_vals, repeat=nste_bnds)]
        return sgrs

    last_sgrs = []
    sgrs = [without_stereo_parities(gra)]

    while sgrs != last_sgrs:
        last_sgrs = sgrs
        sgrs = list(itertools.chain(*map(_expand_atom_stereo, sgrs)))
        sgrs = list(itertools.chain(*map(_expand_bond_stereo, sgrs)))

    return tuple(sorted(sgrs, key=frozen))


def substereomers(gra):
    """ all stereomers compatible with this graph's assignments
    """
    _assigned = functools.partial(
        dict_.filter_by_value, func=lambda x: x is not None)

    known_atm_ste_par_dct = _assigned(atom_stereo_parities(gra))
    known_bnd_ste_par_dct = _assigned(bond_stereo_parities(gra))

    def _is_compatible(sgr):
        atm_ste_par_dct = _assigned(atom_stereo_parities(sgr))
        bnd_ste_par_dct = _assigned(bond_stereo_parities(sgr))
        _compat_atm_assgns = (set(known_atm_ste_par_dct.items()) <=
                              set(atm_ste_par_dct.items()))
        _compat_bnd_assgns = (set(known_bnd_ste_par_dct.items()) <=
                              set(bnd_ste_par_dct.items()))
        return _compat_atm_assgns and _compat_bnd_assgns

    sgrs = tuple(filter(_is_compatible, stereomers(gra)))
    return sgrs


def atom_stereo_sorted_neighbor_atom_keys(gra, atm_key, atm_ngb_keys):
    """ get the neighbor keys of an atom sorted by stereo priority
    """
    atm_ngb_keys = list(atm_ngb_keys)

    # explicitly create an object array because otherwise the argsort
    # interprets [()] as []
    atm_pri_vecs = numpy.empty(len(atm_ngb_keys), dtype=numpy.object_)
    atm_pri_vecs[:] = [stereo_priority_vector(gra, atm_key, atm_ngb_key)
                       for atm_ngb_key in atm_ngb_keys]

    sort_idxs = numpy.argsort(atm_pri_vecs)
    sorted_atm_ngb_keys = tuple(map(atm_ngb_keys.__getitem__, sort_idxs))
    return sorted_atm_ngb_keys


def atoms_stereo_sorted_neighbor_atom_keys(sgr):
    """ Obtain neighbor atom keys for all stereo atoms, sorted by stereo
    priority.

    Includes all stereo atoms and atoms constituting stereo bonds. For stereo
    bonds, the neighbors for each atom in the bond exclude the other atom in
    the bond.

    :param sgr: the graph
    :returns: Neighbor atom keys, sorted by stereo priority, keyed by atom.
    :rtype: dict
    """
    atm_ste_keys = atom_stereo_keys(sgr)
    bnd_ste_keys = bond_stereo_keys(sgr)
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(sgr)

    ste_atm_ngb_keys_dct = {}
    for atm_key in atm_ste_keys:
        atm_ngb_keys = atm_ngb_keys_dct[atm_key]

        ste_atm_ngb_keys_dct[atm_key] = atom_stereo_sorted_neighbor_atom_keys(
            sgr, atm_key, atm_ngb_keys)

    for bnd_key in bnd_ste_keys:
        atm1_key, atm2_key = sorted(bnd_key)

        atm1_ngb_keys = atm_ngb_keys_dct[atm1_key] - bnd_key
        atm2_ngb_keys = atm_ngb_keys_dct[atm2_key] - bnd_key

        ste_atm_ngb_keys_dct[atm1_key] = atom_stereo_sorted_neighbor_atom_keys(
            sgr, atm1_key, atm1_ngb_keys)
        ste_atm_ngb_keys_dct[atm2_key] = atom_stereo_sorted_neighbor_atom_keys(
            sgr, atm2_key, atm2_ngb_keys)

    return ste_atm_ngb_keys_dct


def to_index_based_stereo(sgr):
    """ Convert a graph to index-based stereo assignments, where parities are
    defined relative to the ordering of indices rather than the absolute stereo
    priority.

    :param sgr: a graph with absolute stereo assignments
    :returns: a graph with index-based stereo assignments
    """
    assert sgr == explicit(sgr), (
        "Not an explicit graph:\n{}".format(string(sgr, one_indexed=False)))

    abs_srt_keys_dct = atoms_stereo_sorted_neighbor_atom_keys(sgr)
    atm_ste_keys = atom_stereo_keys(sgr)
    bnd_ste_keys = bond_stereo_keys(sgr)

    abs_atm_ste_par_dct = atom_stereo_parities(sgr)
    abs_bnd_ste_par_dct = bond_stereo_parities(sgr)

    idx_atm_ste_par_dct = {}
    idx_bnd_ste_par_dct = {}

    # Determine index-based stereo assignments for atoms
    for atm_key in atm_ste_keys:
        abs_srt_keys = abs_srt_keys_dct[atm_key]
        idx_srt_keys = sorted(abs_srt_keys)

        if automol.util.is_even_permutation(idx_srt_keys, abs_srt_keys):
            idx_atm_ste_par_dct[atm_key] = abs_atm_ste_par_dct[atm_key]
        else:
            idx_atm_ste_par_dct[atm_key] = not abs_atm_ste_par_dct[atm_key]

    # Determine index-based stereo assignments for bonds
    for bnd_key in bnd_ste_keys:
        atm1_key, atm2_key = sorted(bnd_key)

        atm1_abs_srt_keys = abs_srt_keys_dct[atm1_key]
        atm2_abs_srt_keys = abs_srt_keys_dct[atm2_key]
        atm1_idx_srt_keys = sorted(atm1_abs_srt_keys)
        atm2_idx_srt_keys = sorted(atm2_abs_srt_keys)

        if not ((atm1_idx_srt_keys[0] != atm1_abs_srt_keys[0]) ^
                (atm2_idx_srt_keys[0] != atm2_abs_srt_keys[0])):
            idx_bnd_ste_par_dct[bnd_key] = abs_bnd_ste_par_dct[bnd_key]
        else:
            idx_bnd_ste_par_dct[bnd_key] = not abs_bnd_ste_par_dct[bnd_key]

    sgr = set_atom_stereo_parities(sgr, idx_atm_ste_par_dct)
    sgr = set_bond_stereo_parities(sgr, idx_bnd_ste_par_dct)
    return sgr


def from_index_based_stereo(sgr):
    """ Convert a graph from index-based stereo assignments back to absolute
    stereo assignments, where parities are independent of atom ordering.

    :param sgr: a graph with index-based stereo assignments
    :returns: a graph with absolute stereo assignments
    """
    assert sgr == explicit(sgr), (
        "Not an explicit graph:\n{}".format(string(sgr, one_indexed=False)))

    gra = without_stereo_parities(sgr)

    if has_stereo(sgr):
        atm_keys_pool = atom_stereo_keys(sgr)
        bnd_keys_pool = bond_stereo_keys(sgr)

        idx_atm_ste_par_dct = atom_stereo_parities(sgr)
        idx_bnd_ste_par_dct = bond_stereo_parities(sgr)

        atm_ngb_keys_dct = atoms_neighbor_atom_keys(sgr)

        atm_keys = set()
        bnd_keys = set()

        last_gra = None

        # Do the assignments iteratively to handle higher-order stereo
        while last_gra != gra:
            last_gra = gra

            abs_atm_ste_par_dct = {}
            abs_bnd_ste_par_dct = {}

            atm_keys.update(stereogenic_atom_keys(gra) & atm_keys_pool)
            bnd_keys.update(stereogenic_bond_keys(gra) & bnd_keys_pool)

            # Determine absolute stereo assignments for atoms
            for atm_key in atm_keys:
                abs_srt_keys = atom_stereo_sorted_neighbor_atom_keys(
                    gra, atm_key, atm_ngb_keys_dct[atm_key])
                idx_srt_keys = sorted(abs_srt_keys)

                if automol.util.is_even_permutation(idx_srt_keys,
                                                    abs_srt_keys):
                    abs_atm_ste_par_dct[atm_key] = (
                        idx_atm_ste_par_dct[atm_key])
                else:
                    abs_atm_ste_par_dct[atm_key] = (
                        not idx_atm_ste_par_dct[atm_key])

            # Determine absolute stereo assignments for bonds
            for bnd_key in bnd_keys:
                atm1_key, atm2_key = sorted(bnd_key)

                atm1_abs_srt_keys = atom_stereo_sorted_neighbor_atom_keys(
                    gra, atm1_key, atm_ngb_keys_dct[atm1_key] - bnd_key)
                atm2_abs_srt_keys = atom_stereo_sorted_neighbor_atom_keys(
                    gra, atm2_key, atm_ngb_keys_dct[atm2_key] - bnd_key)
                atm1_idx_srt_keys = sorted(atm1_abs_srt_keys)
                atm2_idx_srt_keys = sorted(atm2_abs_srt_keys)

                if not ((atm1_idx_srt_keys[0] != atm1_abs_srt_keys[0]) ^
                        (atm2_idx_srt_keys[0] != atm2_abs_srt_keys[0])):
                    abs_bnd_ste_par_dct[bnd_key] = (
                        idx_bnd_ste_par_dct[bnd_key])
                else:
                    abs_bnd_ste_par_dct[bnd_key] = (
                        not idx_bnd_ste_par_dct[bnd_key])

            gra = set_atom_stereo_parities(gra, abs_atm_ste_par_dct)
            gra = set_bond_stereo_parities(gra, abs_bnd_ste_par_dct)

    atm_ste_keys = atom_stereo_keys(gra)
    bnd_ste_keys = bond_stereo_keys(gra)
    assert atm_ste_keys == atm_keys_pool, (
        "Index-based to absolute stereo conversion failed:\n"
        "{} != {}".format(str(atm_ste_keys), str(atm_keys_pool)))
    assert bnd_ste_keys == bnd_keys_pool, (
        "Index-based to absolute stereo conversion failed:\n"
        "{} != {}".format(str(bnd_ste_keys), str(bnd_keys_pool)))

    return gra


if __name__ == '__main__':
    # # atom stereo
    # SGR1 = ({0: ('C', 0, None), 1: ('C', 0, True), 2: ('F', 0, None),
    #          3: ('O', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
    #          6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None)},
    #         {frozenset({0, 1}): (1, None), frozenset({0, 4}): (1, None),
    #          frozenset({0, 5}): (1, None), frozenset({0, 6}): (1, None),
    #          frozenset({1, 2}): (1, None), frozenset({1, 3}): (1, None),
    #          frozenset({1, 7}): (1, None), frozenset({8, 3}): (1, None)})
    # # bond stereo
    # SGR1 = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
    #          3: ('F', 0, None), 4: ('O', 0, None), 5: ('H', 0, None),
    #          6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None),
    #          9: ('H', 0, None)},
    #         {frozenset({0, 1}): (1, None), frozenset({0, 5}): (1, None),
    #          frozenset({0, 6}): (1, None), frozenset({0, 7}): (1, None),
    #          frozenset({1, 2}): (1, False), frozenset({8, 1}): (1, None),
    #          frozenset({2, 3}): (1, None), frozenset({2, 4}): (1, None),
    #          frozenset({9, 4}): (1, None)})
    # higher-order stereo
    SGR1 = ({0: ('C', 0, False), 1: ('C', 0, False), 2: ('C', 0, True),
             3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
             6: ('F', 0, None), 7: ('F', 0, None), 8: ('H', 0, None),
             9: ('H', 0, None), 10: ('H', 0, None)},
            {frozenset({1, 9}): (1, None), frozenset({2, 10}): (1, None),
             frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({2, 7}): (1, None), frozenset({0, 5}): (1, None),
             frozenset({2, 4}): (1, None), frozenset({1, 6}): (1, None),
             frozenset({1, 3}): (1, None), frozenset({0, 8}): (1, None)})

    SGR2 = to_index_based_stereo(SGR1)
    SGR1_ = from_index_based_stereo(SGR2)
    print(SGR1_)
    print(SGR1 == SGR1_)
    assert SGR1 == SGR1_
