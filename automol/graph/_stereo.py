""" stereo graph library
"""
import itertools
import functools
import numpy
from automol import dict_
from automol.graph._res import resonance_dominant_atom_hybridizations
from automol.graph._res import resonance_dominant_bond_orders
from automol.graph._ring import rings_bond_keys
from automol.graph._graph import atoms
from automol.graph._graph import bonds
from automol.graph._graph import atom_stereo_parities
from automol.graph._graph import bond_stereo_parities
from automol.graph._graph import set_atom_stereo_parities
from automol.graph._graph import set_bond_stereo_parities
from automol.graph._graph import without_bond_orders
from automol.graph._graph import without_stereo_parities
from automol.graph._graph import frozen
from automol.graph._graph import atom_bond_valences
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import explicit
from automol.graph._graph import implicit
from automol.graph._graph import backbone_keys
from automol.graph._graph import explicit_hydrogen_keys


def has_stereo(gra):
    """ does this graph have stereo of any kind?
    """
    return bool(atom_stereo_keys(gra) or bond_stereo_keys(gra))


def atom_stereo_keys(sgr):
    """ keys to atom stereo-centers
    """
    atm_ste_keys = dict_.keys_by_value(atom_stereo_parities(sgr),
                                       lambda x: x in [True, False])
    return atm_ste_keys


def bond_stereo_keys(sgr):
    """ keys to bond stereo-centers
    """
    bnd_ste_keys = dict_.keys_by_value(bond_stereo_parities(sgr),
                                       lambda x: x in [True, False])
    return bnd_ste_keys


def stereo_priority_vector(gra, atm_key, atm_ngb_key):
    """ generates a sortable one-to-one representation of the branch extending
    from `atm_key` through its bonded neighbor `atm_ngb_key`
    """
    bbn_keys = backbone_keys(gra)
    exp_hyd_keys = explicit_hydrogen_keys(gra)

    if atm_ngb_key not in bbn_keys:
        assert atm_ngb_key in exp_hyd_keys
        assert frozenset({atm_key, atm_ngb_key}) in bonds(gra)
        pri_vec = ()
    else:
        gra = implicit(gra)
        atm_dct = atoms(gra)
        bnd_dct = bonds(gra)
        assert atm_key in bbn_keys
        assert frozenset({atm_key, atm_ngb_key}) in bnd_dct

        # here, switch to an implicit graph
        atm_ngb_keys_dct = atom_neighbor_keys(gra)

        def _priority_vector(atm1_key, atm2_key, seen_keys):
            # we keep a list of seen keys to cut off cycles, avoiding infinite
            # loops

            bnd_val = bnd_dct[frozenset({atm1_key, atm2_key})]
            atm_val = atm_dct[atm2_key]

            bnd_val = _replace_nones_with_negative_infinity(bnd_val)
            atm_val = _replace_nones_with_negative_infinity(atm_val)

            if atm2_key in seen_keys:
                ret = (bnd_val,)
            else:
                seen_keys.update({atm1_key, atm2_key})
                atm3_keys = atm_ngb_keys_dct[atm2_key] - {atm1_key}
                if atm3_keys:
                    next_vals, seen_keys = zip(*[
                        _priority_vector(atm2_key, atm3_key, seen_keys)
                        for atm3_key in atm3_keys])
                    ret = (bnd_val, atm_val) + next_vals
                else:
                    ret = (bnd_val, atm_val)

            return ret, seen_keys

        pri_vec, _ = _priority_vector(atm_key, atm_ngb_key, set())

    return pri_vec


def _replace_nones_with_negative_infinity(seq):
    return [-numpy.inf if val is None else val for val in seq]


def stereogenic_atom_keys(gra):
    """ (unassigned) stereogenic atoms in this graph
    """
    gra = without_bond_orders(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in
    atm_keys = dict_.keys_by_value(atom_bond_valences(gra), lambda x: x == 4)
    atm_keys -= atom_stereo_keys(gra)

    atm_ngb_keys_dct = atom_neighbor_keys(gra)

    def _is_stereogenic(atm_key):
        atm_ngb_keys = list(atm_ngb_keys_dct[atm_key])
        pri_vecs = [stereo_priority_vector(gra, atm_key, atm_ngb_key)
                    for atm_ngb_key in atm_ngb_keys]
        return not any(pv1 == pv2
                       for pv1, pv2 in itertools.combinations(pri_vecs, r=2))

    ste_gen_atm_keys = frozenset(filter(_is_stereogenic, atm_keys))
    return ste_gen_atm_keys


def stereogenic_bond_keys(gra):
    """ (unassigned) stereogenic bonds in this graph
    """
    gra = without_bond_orders(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in

    # get candidates: planar bonds
    bnd_keys = sp2_bond_keys(gra)

    # remove bonds that already have stereo assignments
    bnd_keys -= bond_stereo_keys(gra)
    bnd_keys -= functools.reduce(  # remove double bonds in small rings
        frozenset.union,
        filter(lambda x: len(x) < 8, rings_bond_keys(gra)), frozenset())

    atm_ngb_keys_dct = atom_neighbor_keys(gra)

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


def sp2_bond_keys(gra):
    """ determine the sp2 bonds in this graph
    """
    gra = without_bond_orders(gra)
    bnd_keys = dict_.keys_by_value(
        resonance_dominant_bond_orders(gra), lambda x: 2 in x)

    # make sure both ends are sp^2 (excludes cumulenes)
    atm_hyb_dct = resonance_dominant_atom_hybridizations(gra)
    sp2_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 2)
    bnd_keys = frozenset({bnd_key for bnd_key in bnd_keys
                          if bnd_key <= sp2_atm_keys})
    return bnd_keys


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


def stereo_sorted_atom_neighbor_keys(gra, atm_key, atm_ngb_keys):
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
