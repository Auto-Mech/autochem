""" stereomer expansions
"""
from functools import partial as _partial
from itertools import chain as _chain
from itertools import product as _product
from itertools import starmap as _starmap
from itertools import combinations as _combinations
from .._math import unique as _unique
from .._dict import filter_by_value as _filter_by_value
from .._dict import keys_by_value as _keys_by_value
from .._dict import transform_values as _transform_values
from .._core import frozen as _frozen
from .._core import (atom_implicit_hydrogen_valences as
                     _atom_implicit_hydrogen_valences)
from .._core import atom_stereo_parities as _atom_stereo_parities
from .._core import bond_stereo_parities as _bond_stereo_parities
from .._core import set_atom_stereo_parities as _set_atom_stereo_parities
from .._core import set_bond_stereo_parities as _set_bond_stereo_parities
from .._core import without_bond_orders as _without_bond_orders
from .._core import without_stereo_parities as _without_stereo_parities
from .._graph import atom_bond_keys as _atom_bond_keys
from .._graph import rings_bond_keys as _rings_bond_keys
from .._graph import branch as _branch
from .._expl import backbone_isomorphic as _backbone_isomorphic
from .._expl import explicit as _explicit
from .._res import atom_bond_valences as _atom_bond_valences
from .._res import (resonance_dominant_bond_orders as
                    _resonance_dominant_bond_orders)


# properties
def is_chiral(sgr):
    """ is this stereo graph chiral?

    (ignores partial and higher-order stereo -- the latter should cover most
    cases that aren't highly symmetric)
    """
    _ans = False

    # ignore incomplete and higher-order stereo -- the first is undefined, the
    # second is hard
    if not _is_incomplete_or_higher_order(sgr):
        _ans = not _backbone_isomorphic(sgr, reflection(sgr))

    return _ans


def _is_incomplete_or_higher_order(sgr):
    return (atom_stereo_keys(sgr) !=
            stereogenic_atom_keys(_without_stereo_parities(sgr)))


def atom_stereo_keys(sgr):
    """ keys to atom stereo-centers
    """
    atm_ste_keys = _keys_by_value(_atom_stereo_parities(sgr),
                                  lambda x: x in [True, False])
    return atm_ste_keys


def bond_stereo_keys(sgr):
    """ keys to bond stereo-centers
    """
    bnd_ste_keys = _keys_by_value(_bond_stereo_parities(sgr),
                                  lambda x: x in [True, False])
    return bnd_ste_keys


def stereogenic_atom_keys(xgr):
    """ (unassigned) stereogenic atoms in this graph
    """
    xgr = _without_bond_orders(xgr)
    tet_atm_keys = (
        _keys_by_value(_atom_bond_valences(xgr), lambda x: x == 4) &
        _keys_by_value(_atom_implicit_hydrogen_valences(xgr), lambda x: x < 2)
    )
    atm_keys = tet_atm_keys - atom_stereo_keys(xgr)

    atm_bnd_keys_dct = _atom_bond_keys(xgr)

    def _is_stereogenic(atm_key):
        atm_bnd_keys = atm_bnd_keys_dct[atm_key]
        bnchs = [_branch(xgr, atm_key, bnd_key) for bnd_key in atm_bnd_keys]
        assert len(atm_bnd_keys) in (3, 4)
        _ans = not any(
            _starmap(_backbone_isomorphic, _combinations(bnchs, r=2)))
        return _ans

    ste_gen_atm_keys = frozenset(filter(_is_stereogenic, atm_keys))
    return ste_gen_atm_keys


def stereogenic_bond_keys(xgr):
    """ (unassigned) stereogenic bonds in this graph
    """
    rng_bnd_keys_lst = _rings_bond_keys(xgr)
    sm_rng_bnd_keys_lst = tuple(filter(lambda x: len(x) < 8, rng_bnd_keys_lst))

    def _is_candidate(bnd_key):
        return not any(bnd_key in rng_bnd_keys
                       for rng_bnd_keys in sm_rng_bnd_keys_lst)

    xgr = _without_bond_orders(xgr)
    dbl_bnd_keys = _keys_by_value(_resonance_dominant_bond_orders(xgr),
                                  lambda x: 2 in x)
    bnd_keys = dbl_bnd_keys - bond_stereo_keys(xgr)
    bnd_keys = set(filter(_is_candidate, bnd_keys))
    xgr = _explicit(xgr, set(_chain(*bnd_keys)))

    atm_bnd_keys_dct = _atom_bond_keys(xgr)

    def _is_symmetric_on_bond(atm_key, bnd_key):
        atm_bnd_keys = atm_bnd_keys_dct[atm_key]
        assert bnd_key in atm_bnd_keys
        atm_bnd_keys = atm_bnd_keys - {bnd_key}

        if not atm_bnd_keys:
            _ans = True
        elif len(atm_bnd_keys) == 1:
            _ans = False
        else:
            assert len(atm_bnd_keys) == 2
            bnch1, bnch2 = (
                _branch(xgr, atm_key, atm_bnd_key)
                for atm_bnd_key in atm_bnd_keys)
            _ans = _backbone_isomorphic(bnch1, bnch2)
        return _ans

    def _is_stereogenic(bnd_key):
        return not any(_is_symmetric_on_bond(atm_key, bnd_key)
                       for atm_key in bnd_key)

    ste_gen_bnd_keys = frozenset(filter(_is_stereogenic, bnd_keys))
    return ste_gen_bnd_keys


# transformations
def reflection(sgr):
    """ stereo graph reflection (inverts atom parities)
    """
    atm_pars = _filter_by_value(_atom_stereo_parities(sgr),
                                lambda val: val is not None)
    refl_atm_pars = _transform_values(atm_pars, lambda val: not val)
    return _set_atom_stereo_parities(sgr, refl_atm_pars)


def _explicit_stereo(sgr):
    """ make the hydrogens at atom and bond stereo sites explicit
    """
    atm_ste_keys = atom_stereo_keys(sgr)
    bnd_ste_keys = bond_stereo_keys(sgr)
    bnd_ste_atm_keys = set(_chain(*bnd_ste_keys))
    ste_atm_keys = atm_ste_keys | bnd_ste_atm_keys
    return _explicit(sgr, atm_keys=ste_atm_keys)


def stereomers(xgr):
    """ all stereomers, ignoring this graph's assignments
    """
    bool_vals = (False, True)

    def _expand_atom_stereo(sgr):
        atm_ste_keys = stereogenic_atom_keys(sgr)
        nste_atms = len(atm_ste_keys)
        sgrs = [_set_atom_stereo_parities(sgr, dict(zip(atm_ste_keys,
                                                        atm_ste_par_vals)))
                for atm_ste_par_vals in _product(bool_vals, repeat=nste_atms)]
        return sgrs

    def _expand_bond_stereo(sgr):
        bnd_ste_keys = stereogenic_bond_keys(sgr)
        nste_bnds = len(bnd_ste_keys)
        sgrs = [_set_bond_stereo_parities(sgr, dict(zip(bnd_ste_keys,
                                                        bnd_ste_par_vals)))
                for bnd_ste_par_vals in _product(bool_vals, repeat=nste_bnds)]
        return sgrs

    last_sgrs = []
    sgrs = [_without_stereo_parities(xgr)]

    while sgrs != last_sgrs:
        last_sgrs = sgrs
        sgrs = list(_chain(*map(_expand_atom_stereo, sgrs)))
        sgrs = list(_chain(*map(_expand_bond_stereo, sgrs)))

    return tuple(sorted(sgrs, key=_frozen))


def substereomers(xgr):
    """ all stereomers compatible with this graph's assignments
    """
    _assigned = _partial(_filter_by_value, func=lambda x: x is not None)

    known_atm_ste_par_dct = _assigned(_atom_stereo_parities(xgr))
    known_bnd_ste_par_dct = _assigned(_bond_stereo_parities(xgr))

    def _is_compatible(sgr):
        atm_ste_par_dct = _assigned(_atom_stereo_parities(sgr))
        bnd_ste_par_dct = _assigned(_bond_stereo_parities(sgr))
        _compat_atm_assgns = (set(known_atm_ste_par_dct.items()) <=
                              set(atm_ste_par_dct.items()))
        _compat_bnd_assgns = (set(known_bnd_ste_par_dct.items()) <=
                              set(bnd_ste_par_dct.items()))
        return _compat_atm_assgns and _compat_bnd_assgns

    sgrs = tuple(filter(_is_compatible, stereomers(xgr)))
    return sgrs


# comparisons
def enantiomerically_unique(xgrs):
    """ unique non-isomorphic non-enantiomeric graphs from a series
    """

    def _isomorphic_or_enantiomeric(xgr1, xgr2):
        _ans = _backbone_isomorphic(xgr1, xgr2)
        if not _ans and is_chiral(xgr1):
            _ans = _backbone_isomorphic(xgr1, reflection(xgr2))
        return _ans

    xgrs = _unique(xgrs, equiv=_isomorphic_or_enantiomeric)
    return xgrs
