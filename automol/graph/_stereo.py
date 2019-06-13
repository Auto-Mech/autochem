""" stereo graph library
"""
import itertools
import functools
import more_itertools as mit
import numpy
from qcelemental import constants as qcc
from automol import dict_
from automol import cart
from automol.graph._res import (resonance_dominant_atom_hybridizations as
                                _resonance_dominant_atom_hybridizations)
from automol.graph._res import (resonance_dominant_bond_orders as
                                _resonance_dominant_bond_orders)
from automol.graph._graph import atoms as _atoms
from automol.graph._graph import bonds as _bonds
from automol.graph._graph import atom_keys as _atom_keys
from automol.graph._graph import atom_symbols as _atom_symbols
from automol.graph._graph import atom_stereo_parities as _atom_stereo_parities
from automol.graph._graph import bond_stereo_parities as _bond_stereo_parities
from automol.graph._graph import (set_atom_stereo_parities as
                                  _set_atom_stereo_parities)
from automol.graph._graph import (set_bond_stereo_parities as
                                  _set_bond_stereo_parities)
from automol.graph._graph import without_bond_orders as _without_bond_orders
from automol.graph._graph import (without_stereo_parities as
                                  _without_stereo_parities)
from automol.graph._graph import frozen as _frozen
from automol.graph._graph import atom_bond_valences as _atom_bond_valences
from automol.graph._graph import atom_neighbor_keys as _atom_neighbor_keys
from automol.graph._graph import branch as _branch
from automol.graph._graph import explicit as _explicit
from automol.graph._graph import implicit as _implicit
from automol.graph._graph import backbone_keys as _backbone_keys
from automol.graph._graph import (explicit_hydrogen_keys as
                                  _explicit_hydrogen_keys)
from automol.graph._graph import rings_bond_keys as _rings_bond_keys
from automol.graph._graph import (rings_sorted_atom_keys as
                                  _rings_sorted_atom_keys)
from automol.graph._graph import connected_components as _connected_components


def has_stereo(xgr):
    """ does this graph have stereo of any kind?
    """
    return bool(atom_stereo_keys(xgr) or bond_stereo_keys(xgr))


def atom_stereo_keys(sgr):
    """ keys to atom stereo-centers
    """
    atm_ste_keys = dict_.keys_by_value(_atom_stereo_parities(sgr),
                                       lambda x: x in [True, False])
    return atm_ste_keys


def bond_stereo_keys(sgr):
    """ keys to bond stereo-centers
    """
    bnd_ste_keys = dict_.keys_by_value(_bond_stereo_parities(sgr),
                                       lambda x: x in [True, False])
    return bnd_ste_keys


def stereo_priority_vector(xgr, atm_key, atm_ngb_key):
    """ generates a sortable one-to-one representation of the branch extending
    from `atm_key` through its bonded neighbor `atm_ngb_key`
    """
    bbn_keys = _backbone_keys(xgr)
    exp_hyd_keys = _explicit_hydrogen_keys(xgr)

    if atm_ngb_key not in bbn_keys:
        assert atm_ngb_key in exp_hyd_keys
        assert frozenset({atm_key, atm_ngb_key}) in _bonds(xgr)
        pri_vec = ()
    else:
        xgr = _implicit(xgr)
        atm_dct = _atoms(xgr)
        bnd_dct = _bonds(xgr)
        assert atm_key in bbn_keys
        assert frozenset({atm_key, atm_ngb_key}) in bnd_dct

        # here, switch to an implicit graph
        atm_ngb_keys_dct = _atom_neighbor_keys(xgr)

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


def stereogenic_atom_keys(xgr):
    """ (unassigned) stereogenic atoms in this graph
    """
    xgr = _without_bond_orders(xgr)
    xgr = _explicit(xgr)  # for simplicity, add the explicit hydrogens back in
    atm_keys = dict_.keys_by_value(_atom_bond_valences(xgr), lambda x: x == 4)
    atm_keys -= atom_stereo_keys(xgr)

    atm_ngb_keys_dct = _atom_neighbor_keys(xgr)

    def _is_stereogenic(atm_key):
        atm_ngb_keys = list(atm_ngb_keys_dct[atm_key])
        pri_vecs = [stereo_priority_vector(xgr, atm_key, atm_ngb_key)
                    for atm_ngb_key in atm_ngb_keys]
        return not any(pv1 == pv2
                       for pv1, pv2 in itertools.combinations(pri_vecs, r=2))

    ste_gen_atm_keys = frozenset(filter(_is_stereogenic, atm_keys))
    return ste_gen_atm_keys


def stereogenic_bond_keys(xgr):
    """ (unassigned) stereogenic bonds in this graph
    """
    xgr = _without_bond_orders(xgr)
    xgr = _explicit(xgr)  # for simplicity, add the explicit hydrogens back in
    bnd_keys = dict_.keys_by_value(
        _resonance_dominant_bond_orders(xgr), lambda x: 2 in x)

    bnd_keys -= bond_stereo_keys(xgr)
    bnd_keys -= functools.reduce(  # remove double bonds in small rings
        frozenset.union,
        filter(lambda x: len(x) < 8, _rings_bond_keys(xgr)), frozenset())

    atm_ngb_keys_dct = _atom_neighbor_keys(xgr)

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
                ret = (stereo_priority_vector(xgr, atm_key, atm_ngb_keys[0]) ==
                       stereo_priority_vector(xgr, atm_key, atm_ngb_keys[1]))

            return ret

        return not (_is_symmetric_on_bond(atm1_key, atm2_key) or
                    _is_symmetric_on_bond(atm2_key, atm1_key))

    ste_gen_bnd_keys = frozenset(filter(_is_stereogenic, bnd_keys))
    return ste_gen_bnd_keys


def stereomers(xgr):
    """ all stereomers, ignoring this graph's assignments
    """
    bool_vals = (False, True)

    def _expand_atom_stereo(sgr):
        atm_ste_keys = stereogenic_atom_keys(sgr)
        nste_atms = len(atm_ste_keys)
        sgrs = [_set_atom_stereo_parities(sgr, dict(zip(atm_ste_keys,
                                                        atm_ste_par_vals)))
                for atm_ste_par_vals
                in itertools.product(bool_vals, repeat=nste_atms)]
        return sgrs

    def _expand_bond_stereo(sgr):
        bnd_ste_keys = stereogenic_bond_keys(sgr)
        nste_bnds = len(bnd_ste_keys)
        sgrs = [_set_bond_stereo_parities(sgr, dict(zip(bnd_ste_keys,
                                                        bnd_ste_par_vals)))
                for bnd_ste_par_vals
                in itertools.product(bool_vals, repeat=nste_bnds)]
        return sgrs

    last_sgrs = []
    sgrs = [_without_stereo_parities(xgr)]

    while sgrs != last_sgrs:
        last_sgrs = sgrs
        sgrs = list(itertools.chain(*map(_expand_atom_stereo, sgrs)))
        sgrs = list(itertools.chain(*map(_expand_bond_stereo, sgrs)))

    return tuple(sorted(sgrs, key=_frozen))


def substereomers(xgr):
    """ all stereomers compatible with this graph's assignments
    """
    _assigned = functools.partial(
        dict_.filter_by_value, func=lambda x: x is not None)

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


def stereo_sorted_atom_neighbor_keys(xgr, atm_key, atm_ngb_keys):
    """ get the neighbor keys of an atom sorted by stereo priority
    """
    atm_ngb_keys = list(atm_ngb_keys)

    # explicitly create an object array because otherwise the argsort
    # interprets [()] as []
    atm_pri_vecs = numpy.empty(len(atm_ngb_keys), dtype=numpy.object_)
    atm_pri_vecs[:] = [stereo_priority_vector(xgr, atm_key, atm_ngb_key)
                       for atm_ngb_key in atm_ngb_keys]

    sort_idxs = numpy.argsort(atm_pri_vecs)
    sorted_atm_ngb_keys = tuple(map(atm_ngb_keys.__getitem__, sort_idxs))
    return sorted_atm_ngb_keys


def set_stereo_from_atom_coordinates(xgr, atm_xyz_dct):
    """ set atom and bond stereo parities using a set of atomic coordinates

    (coordinate distances need not match connectivity -- what matters is the
    relative positions at stereo sites)
    """
    last_xgr = None

    # set atom and bond stereo, iterating to self-consistency
    atm_keys = set()
    bnd_keys = set()
    while last_xgr != xgr:
        last_xgr = xgr
        atm_keys.update(stereogenic_atom_keys(xgr))
        bnd_keys.update(stereogenic_bond_keys(xgr))
        xgr = _set_atom_stereo_from_coordinates(xgr, atm_keys, atm_xyz_dct)
        xgr = _set_bond_stereo_from_coordinates(xgr, bnd_keys, atm_xyz_dct)

    return xgr


def _set_atom_stereo_from_coordinates(xgr, atm_keys, atm_xyz_dct):
    assert xgr == _explicit(xgr)

    atm_pars = [
        _atom_stereo_parity_from_coordinates(xgr, atm_key, atm_xyz_dct)
        for atm_key in atm_keys]
    xgr = _set_atom_stereo_parities(xgr, dict(zip(atm_keys, atm_pars)))
    return xgr


def _set_bond_stereo_from_coordinates(xgr, bnd_keys, atm_xyz_dct):
    assert xgr == _explicit(xgr)

    bnd_pars = [
        _bond_stereo_parity_from_coordinates(xgr, bnd_key, atm_xyz_dct)
        for bnd_key in bnd_keys]
    xgr = _set_bond_stereo_parities(xgr, dict(zip(bnd_keys, bnd_pars)))
    return xgr


def _atom_stereo_parity_from_coordinates(xgr, atm_key, atm_xyz_dct):
    atm_ngb_keys_dct = _atom_neighbor_keys(xgr)
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]
    atm_ngb_keys = stereo_sorted_atom_neighbor_keys(
        xgr, atm_key, atm_ngb_keys)
    atm_ngb_xyzs = list(map(atm_xyz_dct.__getitem__, atm_ngb_keys))
    det_mat = numpy.ones((4, 4))
    det_mat[:, :3] = atm_ngb_xyzs
    det_val = numpy.linalg.det(det_mat)
    assert det_val != 0.  # for now, assume no four-atom planes
    par = det_val > 0.
    return par


def _bond_stereo_parity_from_coordinates(xgr, bnd_key, atm_xyz_dct):
    atm1_key, atm2_key = bnd_key
    atm_ngb_keys_dct = _atom_neighbor_keys(xgr)
    atm1_ngb_keys = atm_ngb_keys_dct[atm1_key] - {atm2_key}
    atm2_ngb_keys = atm_ngb_keys_dct[atm2_key] - {atm1_key}

    atm1_ngb_keys = stereo_sorted_atom_neighbor_keys(
        xgr, atm1_key, atm1_ngb_keys)
    atm2_ngb_keys = stereo_sorted_atom_neighbor_keys(
        xgr, atm2_key, atm2_ngb_keys)

    # get the top priority neighbor keys on each side
    atm1_ngb_key = atm1_ngb_keys[0]
    atm2_ngb_key = atm2_ngb_keys[0]

    # get the bond direction vectors for each
    atm1_bnd_vec = numpy.subtract(atm_xyz_dct[atm1_ngb_key],
                                  atm_xyz_dct[atm1_key])
    atm2_bnd_vec = numpy.subtract(atm_xyz_dct[atm2_ngb_key],
                                  atm_xyz_dct[atm2_key])
    dot_val = numpy.vdot(atm1_bnd_vec, atm2_bnd_vec)
    assert dot_val != 0.  # for now, assume no collinear
    par = dot_val > 0.
    return par


def atom_stereo_coordinates(sgr):
    """ stereo-specific coordinates for this molecular graph
    """
    assert sgr == _explicit(sgr)
    last_xgr = None
    xgr = _without_stereo_parities(sgr)

    full_atm_ste_par_dct = _atom_stereo_parities(sgr)
    full_bnd_ste_par_dct = _bond_stereo_parities(sgr)

    # first, get a set of non-stereo-specific coordinates for the graph
    atm_xyz_dct = _atom_coordinates(xgr)

    atm_keys = set()
    bnd_keys = set()

    while last_xgr != xgr:
        last_xgr = xgr
        atm_keys.update(stereogenic_atom_keys(xgr))
        bnd_keys.update(stereogenic_bond_keys(xgr))
        atm_ste_par_dct = {atm_key: full_atm_ste_par_dct[atm_key]
                           for atm_key in atm_keys}
        bnd_ste_par_dct = {bnd_key: full_bnd_ste_par_dct[bnd_key]
                           for bnd_key in bnd_keys}
        xgr, atm_xyz_dct = _correct_atom_stereo_coordinates(
            xgr, atm_ste_par_dct, atm_xyz_dct)
        xgr, atm_xyz_dct = _correct_bond_stereo_coordinates(
            xgr, bnd_ste_par_dct, atm_xyz_dct)

    return atm_xyz_dct


def _correct_atom_stereo_coordinates(xgr, atm_ste_par_dct, atm_xyz_dct):
    ring_atm_keys = set(itertools.chain(*_rings_sorted_atom_keys(xgr)))
    atm_ngb_keys_dct = _atom_neighbor_keys(xgr)

    atm_keys = list(atm_ste_par_dct.keys())
    for atm_key in atm_keys:
        par = atm_ste_par_dct[atm_key]
        curr_par = _atom_stereo_parity_from_coordinates(
            xgr, atm_key, atm_xyz_dct)
        atm_ngb_keys = atm_ngb_keys_dct[atm_key]

        if curr_par != par:
            # for now, we simply exclude rings from the pivot keys
            # (will not work for stereo atom at the intersection of two rings)
            atm_piv_keys = list(atm_ngb_keys - ring_atm_keys)[:2]
            assert len(atm_piv_keys) == 2
            atm3_key, atm4_key = atm_piv_keys
            atm_xyz = atm_xyz_dct[atm_key]
            atm3_xyz = atm_xyz_dct[atm3_key]
            atm4_xyz = atm_xyz_dct[atm4_key]
            rot_axis = cart.vec.unit_bisector(
                atm3_xyz, atm4_xyz, orig_xyz=atm_xyz)
            rot_ = cart.vec.rotate_(
                rot_axis, numpy.pi, orig_xyz=atm_xyz)

            rot_atm_keys = list(
                _atom_keys(_branch(xgr, atm_key, {atm_key, atm3_key})) |
                _atom_keys(_branch(xgr, atm_key, {atm_key, atm4_key})))

            rot_atm_xyzs = list(
                map(rot_, map(atm_xyz_dct.__getitem__, rot_atm_keys)))

            atm_xyz_dct.update(dict(zip(rot_atm_keys, rot_atm_xyzs)))

        new_par = _atom_stereo_parity_from_coordinates(
            xgr, atm_key, atm_xyz_dct)

        assert new_par == par
        xgr = _set_atom_stereo_parities(xgr, {atm_key: par})

    return xgr, atm_xyz_dct


def _correct_bond_stereo_coordinates(xgr, bnd_ste_par_dct, atm_xyz_dct):
    atm_ngb_keys_dct = _atom_neighbor_keys(xgr)

    bnd_keys = list(bnd_ste_par_dct.keys())
    for bnd_key in bnd_keys:
        par = bnd_ste_par_dct[bnd_key]
        curr_par = _bond_stereo_parity_from_coordinates(
            xgr, bnd_key, atm_xyz_dct)

        if curr_par != par:
            atm1_key, atm2_key = bnd_key
            atm1_xyz = atm_xyz_dct[atm1_key]
            atm2_xyz = atm_xyz_dct[atm2_key]

            atm1_ngb_keys = atm_ngb_keys_dct[atm1_key] - {atm2_key}
            atm2_ngb_keys = atm_ngb_keys_dct[atm2_key] - {atm1_key}
            atm1_ngb_keys = stereo_sorted_atom_neighbor_keys(
                xgr, atm1_key, atm1_ngb_keys)
            atm2_ngb_keys = stereo_sorted_atom_neighbor_keys(
                xgr, atm2_key, atm2_ngb_keys)

            rot_axis = numpy.subtract(atm2_xyz, atm1_xyz)
            rot_ = cart.vec.rotate_(
                rot_axis, numpy.pi, orig_xyz=atm1_xyz)

            rot_atm_keys = _atom_keys(
                _branch(xgr, atm2_key, {atm2_key, atm2_ngb_keys[0]}))

            if len(atm2_ngb_keys) > 1:
                assert len(atm2_ngb_keys) == 2
                rot_atm_keys |= _atom_keys(
                    _branch(xgr, atm2_key, {atm2_key, atm2_ngb_keys[1]}))

            rot_atm_keys = list(rot_atm_keys)

            rot_atm_xyzs = list(
                map(rot_, map(atm_xyz_dct.__getitem__, rot_atm_keys)))

            atm_xyz_dct.update(dict(zip(rot_atm_keys, rot_atm_xyzs)))

        assert _bond_stereo_parity_from_coordinates(
            xgr, bnd_key, atm_xyz_dct) == par
        xgr = _set_bond_stereo_parities(xgr, {bnd_key: par})

    return xgr, atm_xyz_dct


def _atom_coordinates(sgr):
    """ non-stereo-specific coordinates for a molecular graph
    """
    atm_xyz_dct = {}
    for idx, cnn_sgr in enumerate(_connected_components(sgr)):
        shift = 20. * idx
        cnn_atm_xyz_dct = _connected_graph_atom_coordinates(cnn_sgr)
        atm_keys = list(cnn_atm_xyz_dct.keys())
        atm_xyzs = numpy.array(list(cnn_atm_xyz_dct.values()))
        atm_xyzs += numpy.array([0., 0., shift])
        atm_xyz_dct.update(dict(zip(atm_keys, map(tuple, atm_xyzs))))
    return atm_xyz_dct


def _connected_graph_atom_coordinates(sgr):
    """ non-stereo-specific coordinates for a connected molecular graph

    (currently assumes a with at most one ring -- fix that)
    """
    assert sgr == _explicit(sgr)

    atm_keys = _atom_keys(sgr)
    rng_atm_keys_lst = _rings_sorted_atom_keys(sgr)

    if len(atm_keys) == 1:
        atm1_key, = atm_keys
        atm_xyz_dct = {}
        atm_xyz_dct[atm1_key] = (0., 0., 0.)
    elif len(atm_keys) == 2:
        atm1_key, atm2_key = atm_keys

        atm1_xyz = (0., 0., 0.)

        dist = _bond_distance(sgr, atm2_key, atm1_key)
        atm2_xyz = cart.vec.from_internals(dist=dist, xyz1=atm1_xyz)
        atm_xyz_dct = {}
        atm_xyz_dct[atm1_key] = tuple(atm1_xyz)
        atm_xyz_dct[atm2_key] = tuple(atm2_xyz)
    elif not rng_atm_keys_lst:
        atm_ngb_keys_dct = _atom_neighbor_keys(sgr)

        # grab the first three coordinates along the longest chain
        max_chain = longest_chain(sgr)

        atm1_key, atm2_key, atm3_key = max_chain[:3]

        atm1_xyz = (0., 0., 0.)

        dist = _bond_distance(sgr, atm2_key, atm1_key)
        atm2_xyz = cart.vec.from_internals(dist=dist, xyz1=atm1_xyz)

        dist = _bond_distance(sgr, atm3_key, atm2_key)
        ang = _bond_angle(sgr, atm3_key, atm2_key, atm1_key)
        atm3_xyz = cart.vec.from_internals(dist=dist, xyz1=atm2_xyz,
                                           ang=ang, xyz2=atm1_xyz)

        atm_xyz_dct = {}
        atm_xyz_dct[atm1_key] = tuple(atm1_xyz)
        atm_xyz_dct[atm2_key] = tuple(atm2_xyz)
        atm_xyz_dct[atm3_key] = tuple(atm3_xyz)

        atm_xyz_dct = _extend_atom_coordinates(
            sgr, atm1_key, atm2_key, atm3_key, atm_xyz_dct)

        atm_ngb_keys_dct = _atom_neighbor_keys(sgr)
        atm4_key = sorted(atm_ngb_keys_dct[atm3_key] - {atm2_key})[0]
        atm_xyz_dct = _extend_atom_coordinates(
            sgr, atm4_key, atm3_key, atm2_key, atm_xyz_dct)
    elif len(rng_atm_keys_lst) == 1:
        # for now, we'll assume only one ring
        rng_atm_keys, = rng_atm_keys_lst
        rng_atm_xyzs = _polygon_coordinates(num=len(rng_atm_keys))

        atm_xyz_dct = dict(zip(rng_atm_keys, rng_atm_xyzs))

        assert len(rng_atm_keys) >= 3
        trip_iter = mit.windowed(rng_atm_keys[-2:] + rng_atm_keys, 3)
        for atm1_key, atm2_key, atm3_key in trip_iter:
            atm_xyz_dct = _extend_atom_coordinates(
                sgr, atm1_key, atm2_key, atm3_key, atm_xyz_dct)
    else:
        raise NotImplementedError("This algorithm is currently not implemented"
                                  "for more than one ring.")

    return atm_xyz_dct


def _polygon_coordinates(num):
    """ main formula: side / 2 = radius * sin(2 * pi / 2 / 2)
    """
    side = 1.5 * qcc.conversion_factor('angstrom', 'bohr')
    rad = side / 2. / numpy.sin(numpy.pi / num)
    angs = [2. * numpy.pi * idx / num for idx in range(num)]
    xyzs = [(rad * numpy.cos(ang), rad * numpy.sin(ang), 0.) for ang in angs]
    return tuple(xyzs)


def _extend_atom_coordinates(sgr, atm1_key, atm2_key, atm3_key, atm_xyz_dct):
    assert {atm1_key, atm2_key, atm3_key} <= set(atm_xyz_dct)
    atm1_xyz = atm_xyz_dct[atm1_key]
    atm2_xyz = atm_xyz_dct[atm2_key]
    atm3_xyz = atm_xyz_dct[atm3_key]

    atm_ngb_keys_dct = _atom_neighbor_keys(sgr)

    fix_atm_keys = set(atm_xyz_dct)

    all_atm4_keys = atm_ngb_keys_dct[atm3_key] - {atm2_key}
    fix_atm4_keys = list(all_atm4_keys & fix_atm_keys)
    atm4_keys = list(all_atm4_keys - fix_atm_keys)

    # set the initial dihedral angle value
    dih_incr = _dihedral_increment(sgr, atm1_key, atm2_key, atm3_key)
    if fix_atm4_keys:
        fix_atm4_dih_vals = list(
            cart.vec.dihedral_angle(
                atm1_xyz, atm2_xyz, atm3_xyz, atm_xyz_dct[fix_atm4_key])
            for fix_atm4_key in fix_atm4_keys)
        dih = max(fix_atm4_dih_vals)
    else:
        dih = numpy.pi - dih_incr

    for atm4_key in atm4_keys:
        dist = _bond_distance(sgr, atm3_key, atm4_key)
        ang = _bond_angle(sgr, atm2_key, atm3_key, atm4_key)
        dih += dih_incr

        atm4_xyz = cart.vec.from_internals(dist=dist, xyz1=atm3_xyz,
                                           ang=ang, xyz2=atm2_xyz,
                                           dih=dih, xyz3=atm1_xyz)
        atm_xyz_dct[atm4_key] = atm4_xyz

        atm_xyz_dct = _extend_atom_coordinates(
            sgr, atm2_key, atm3_key, atm4_key, atm_xyz_dct)

    return atm_xyz_dct


def _bond_distance(xgr, atm1_key, atm2_key):
    """ predicted bond distance

    (currently crude, but could easily be made more sophisticated
    """
    atm_sym_dct = _atom_symbols(xgr)
    assert atm2_key in _atom_neighbor_keys(xgr)[atm1_key]

    atm1_sym = atm_sym_dct[atm1_key]
    atm2_sym = atm_sym_dct[atm2_key]

    if atm1_sym == 'H' or atm2_sym == 'H':
        dist = 1.1 * qcc.conversion_factor('angstrom', 'bohr')
    else:
        dist = 1.5 * qcc.conversion_factor('angstrom', 'bohr')

    return dist


def _bond_angle(sgr, atm1_key, atm2_key, atm3_key):
    """ predict the bond angles an atom makes with its neighbors
    """
    atm_ngb_keys_dct = _atom_neighbor_keys(sgr)
    atm2_ngb_keys = atm_ngb_keys_dct[atm2_key]
    assert {atm1_key, atm3_key} <= atm2_ngb_keys

    atm_hyb_dct = _resonance_dominant_atom_hybridizations(sgr)
    atm2_hyb = atm_hyb_dct[atm2_key]

    if atm2_hyb == 3:
        ang = 109.5 * qcc.conversion_factor('degree', 'radian')
    elif atm2_hyb == 2:
        ang = 120.0 * qcc.conversion_factor('degree', 'radian')
    else:
        assert atm2_hyb == 1
        ang = 180.0 * qcc.conversion_factor('degree', 'radian')

    return ang


def _dihedral_increment(sgr, atm1_key, atm2_key, atm3_key):
    """ predict dihedral increment for atoms attached to `atm3_key`
    """
    atm_ngb_keys_dct = _atom_neighbor_keys(sgr)
    atm2_ngb_keys = atm_ngb_keys_dct[atm2_key]
    atm3_ngb_keys = atm_ngb_keys_dct[atm3_key]
    assert atm2_key in atm3_ngb_keys
    assert atm1_key in atm2_ngb_keys - {atm3_key}

    atm_hyb_dct = _resonance_dominant_atom_hybridizations(sgr)
    atm3_hyb = atm_hyb_dct[atm3_key]
    dih_incr = 2 * numpy.pi / atm3_hyb if not atm3_hyb == 0 else 0.
    return dih_incr


def longest_chain(xgr):
    """ longest chain in the graph
    """
    atm_keys = _atom_keys(xgr)

    max_chain = max((_longest_chain(xgr, atm_key) for atm_key in atm_keys),
                    key=len)
    return max_chain


def _longest_chain(xgr, atm_key):
    atm_ngb_keys_dct = _atom_neighbor_keys(xgr)
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    chains_lst = []
    next_chains_lst = [[atm_key, atm_ngb_key] for atm_ngb_key in atm_ngb_keys]

    while True:
        chains_lst = next_chains_lst
        next_chains_lst = []
        for chain in chains_lst:
            atm_ngb_keys = atm_ngb_keys_dct[chain[-1]]
            next_atm_keys = sorted(atm_ngb_keys - set(chain))
            for next_atm_key in next_atm_keys:
                next_chains_lst.append(chain + [next_atm_key])

        if not next_chains_lst:
            break

    max_chain = tuple(chains_lst[0])
    return max_chain
