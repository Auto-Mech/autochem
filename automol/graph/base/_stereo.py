""" graph functions that depend on stereo assignments

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import functools
import itertools
import numpy
from automol import util
import automol.geom.base    # !!!!
from automol.util import dict_
from automol.graph.base._core import atoms
from automol.graph.base._core import bonds
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bond_keys
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import set_atom_stereo_parities
from automol.graph.base._core import set_bond_stereo_parities
from automol.graph.base._core import string
from automol.graph.base._core import frozen
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import has_stereo
from automol.graph.base._core import atom_bond_valences
from automol.graph.base._core import backbone_keys
from automol.graph.base._core import explicit_hydrogen_keys
from automol.graph.base._core import without_bond_orders
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._core import implicit
from automol.graph.base._core import explicit
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._algo import rings_atom_keys
from automol.graph.base._algo import rings_bond_keys
from automol.graph.base._algo import branch
from automol.graph.base._resonance import sp2_bond_keys


# # core functions
def stereo_priority_vector(gra, atm_key, atm_ngb_key):
    """ generates a sortable one-to-one representation of the branch extending
    from `atm_key` through its bonded neighbor `atm_ngb_key`
    """
    bbn_keys = backbone_keys(gra)
    exp_hyd_keys = explicit_hydrogen_keys(gra)

    if atm_ngb_key not in bbn_keys:
        assert atm_ngb_key in exp_hyd_keys
        assert frozenset({atm_key, atm_ngb_key}) in bond_keys(gra)
        pri_vec = ()
    else:
        gra = implicit(gra)
        atm_dct = atoms(gra)
        bnd_dct = bonds(gra)
        assert atm_key in bbn_keys
        assert frozenset({atm_key, atm_ngb_key}) in bnd_dct

        # here, switch to an implicit graph
        atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

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


# # index-based stereo conversions
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

        if util.is_even_permutation(idx_srt_keys, abs_srt_keys):
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

                if util.is_even_permutation(idx_srt_keys, abs_srt_keys):
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


# # derived properties
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


# # stereo setting code
def set_stereo_from_geometry(gra, geo, geo_idx_dct=None):
    """ set graph stereo from a geometry

    (coordinate distances need not match connectivity -- what matters is the
    relative positions at stereo sites)
    """
    gra = without_stereo_parities(gra)
    last_gra = None

    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {atm_key: idx for idx, atm_key in enumerate(atm_keys)})

    # set atom and bond stereo, iterating to self-consistency
    atm_keys = set()
    bnd_keys = set()
    while last_gra != gra:
        last_gra = gra
        atm_keys.update(stereogenic_atom_keys(gra))
        bnd_keys.update(stereogenic_bond_keys(gra))
        gra = _set_atom_stereo_from_geometry(gra, atm_keys, geo, geo_idx_dct)
        gra = _set_bond_stereo_from_geometry(gra, bnd_keys, geo, geo_idx_dct)

    return gra


# # stereo parity evaluation code
def atom_stereo_parity_from_geometry(gra, atm_key, geo, geo_idx_dct):
    """ get the current stereo parity of an atom from its geometry
    """
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    # sort the neighbor keys by stereo priority
    atm_ngb_keys = atom_stereo_sorted_neighbor_atom_keys(
        gra, atm_key, atm_ngb_keys)

    # determine the parity based on the coordinates
    xyzs = automol.geom.base.coordinates(geo)
    atm_ngb_idxs = dict_.values_by_key(geo_idx_dct, atm_ngb_keys)
    atm_ngb_xyzs = [xyzs[idx] for idx in atm_ngb_idxs]
    det_mat = numpy.ones((4, 4))
    det_mat[:, :3] = atm_ngb_xyzs
    det_val = numpy.linalg.det(det_mat)
    assert det_val != 0.  # for now, assume no four-atom planes
    par = det_val > 0.
    return par


def bond_stereo_parity_from_geometry(gra, bnd_key, geo, geo_idx_dct):
    """ get the current stereo parity of a bond from its geometry
    """
    atm1_key, atm2_key = bnd_key
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    atm1_ngb_keys = atm_ngb_keys_dct[atm1_key] - {atm2_key}
    atm2_ngb_keys = atm_ngb_keys_dct[atm2_key] - {atm1_key}

    atm1_ngb_keys = atom_stereo_sorted_neighbor_atom_keys(
        gra, atm1_key, atm1_ngb_keys)
    atm2_ngb_keys = atom_stereo_sorted_neighbor_atom_keys(
        gra, atm2_key, atm2_ngb_keys)

    # get the top priority neighbor keys on each side
    atm1_ngb_key = atm1_ngb_keys[0]
    atm2_ngb_key = atm2_ngb_keys[0]

    # determine the parity based on the coordinates
    xyzs = automol.geom.base.coordinates(geo)
    atm1_xyz = xyzs[geo_idx_dct[atm1_key]]
    atm2_xyz = xyzs[geo_idx_dct[atm2_key]]
    atm1_ngb_xyz = xyzs[geo_idx_dct[atm1_ngb_key]]
    atm2_ngb_xyz = xyzs[geo_idx_dct[atm2_ngb_key]]
    atm1_bnd_vec = numpy.subtract(atm1_ngb_xyz, atm1_xyz)
    atm2_bnd_vec = numpy.subtract(atm2_ngb_xyz, atm2_xyz)
    dot_val = numpy.vdot(atm1_bnd_vec, atm2_bnd_vec)
    assert dot_val != 0.  # for now, assume no collinear
    par = dot_val > 0.
    return par


# helpers
def _set_atom_stereo_from_geometry(gra, atm_keys, geo, geo_idx_dct):
    assert gra == explicit(gra)

    atm_pars = [
        atom_stereo_parity_from_geometry(gra, atm_key, geo, geo_idx_dct)
        for atm_key in atm_keys]
    gra = set_atom_stereo_parities(gra, dict(zip(atm_keys, atm_pars)))
    return gra


def _set_bond_stereo_from_geometry(gra, bnd_keys, geo, geo_idx_dct):
    assert gra == explicit(gra)

    bnd_pars = [
        bond_stereo_parity_from_geometry(gra, bnd_key, geo, geo_idx_dct)
        for bnd_key in bnd_keys]
    gra = set_bond_stereo_parities(gra, dict(zip(bnd_keys, bnd_pars)))
    return gra


def _stereo_corrected_geometry(sgr, geo, geo_idx_dct):
    """ correct the stereo parities of a geometry

    (works iterately to handle cases of higher-order stereo)
    """
    assert sgr == explicit(sgr)
    gra = without_stereo_parities(sgr)

    if has_stereo(sgr):
        full_atm_ste_par_dct = atom_stereo_parities(sgr)
        full_bnd_ste_par_dct = bond_stereo_parities(sgr)

        atm_keys = set()
        bnd_keys = set()

        last_gra = None

        while last_gra != gra:
            last_gra = gra

            atm_keys.update(stereogenic_atom_keys(gra))
            bnd_keys.update(stereogenic_bond_keys(gra))

            atm_ste_par_dct = {atm_key: full_atm_ste_par_dct[atm_key]
                               for atm_key in atm_keys}
            bnd_ste_par_dct = {bnd_key: full_bnd_ste_par_dct[bnd_key]
                               for bnd_key in bnd_keys}
            geo, gra = _atom_stereo_corrected_geometry(
                gra, atm_ste_par_dct, geo, geo_idx_dct)
            geo, gra = _bond_stereo_corrected_geometry(
                gra, bnd_ste_par_dct, geo, geo_idx_dct)

    return geo


def _atom_stereo_corrected_geometry(gra, atm_ste_par_dct, geo, geo_idx_dct):
    """ correct the atom stereo parities of a geometry, for a subset of atoms
    """
    ring_atm_keys = set(itertools.chain(*rings_atom_keys(gra)))
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    atm_keys = list(atm_ste_par_dct.keys())
    for atm_key in atm_keys:
        par = atm_ste_par_dct[atm_key]
        curr_par = atom_stereo_parity_from_geometry(
            gra, atm_key, geo, geo_idx_dct)

        if curr_par != par:
            atm_ngb_keys = atm_ngb_keys_dct[atm_key]
            # for now, we simply exclude rings from the pivot keys
            # (will not work for stereo atom at the intersection of two rings)
            atm_piv_keys = list(atm_ngb_keys - ring_atm_keys)[:2]
            assert len(atm_piv_keys) == 2
            atm3_key, atm4_key = atm_piv_keys

            # get coordinates
            xyzs = automol.geom.base.coordinates(geo)
            atm_xyz = xyzs[geo_idx_dct[atm_key]]
            atm3_xyz = xyzs[geo_idx_dct[atm3_key]]
            atm4_xyz = xyzs[geo_idx_dct[atm4_key]]

            # do the rotation
            rot_axis = util.vec.unit_bisector(
                atm3_xyz, atm4_xyz, orig_xyz=atm_xyz)

            rot_atm_keys = (
                atom_keys(branch(gra, atm_key, {atm_key, atm3_key})) |
                atom_keys(branch(gra, atm_key, {atm_key, atm4_key})))
            rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

            geo = automol.geom.rotate(
                geo, rot_axis, numpy.pi, orig_xyz=atm_xyz, idxs=rot_idxs)

        assert atom_stereo_parity_from_geometry(
            gra, atm_key, geo, geo_idx_dct) == par
        gra = set_atom_stereo_parities(gra, {atm_key: par})

    return geo, gra


def _bond_stereo_corrected_geometry(gra, bnd_ste_par_dct, geo, geo_idx_dct):
    """ correct the bond stereo parities of a geometry, for a subset of bonds
    """
    bnd_keys = list(bnd_ste_par_dct.keys())
    for bnd_key in bnd_keys:
        par = bnd_ste_par_dct[bnd_key]
        curr_par = bond_stereo_parity_from_geometry(
            gra, bnd_key, geo, geo_idx_dct)

        if curr_par != par:
            xyzs = automol.geom.base.coordinates(geo)

            atm1_key, atm2_key = bnd_key
            atm1_xyz = xyzs[geo_idx_dct[atm1_key]]
            atm2_xyz = xyzs[geo_idx_dct[atm2_key]]

            rot_axis = numpy.subtract(atm2_xyz, atm1_xyz)

            rot_atm_keys = atom_keys(
                branch(gra, atm1_key, {atm1_key, atm2_key}))

            rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

            geo = automol.geom.rotate(
                geo, rot_axis, numpy.pi, orig_xyz=atm1_xyz, idxs=rot_idxs)

        assert bond_stereo_parity_from_geometry(
            gra, bnd_key, geo, geo_idx_dct) == par
        gra = set_bond_stereo_parities(gra, {bnd_key: par})

    return geo, gra


def _replace_nones_with_negative_infinity(seq):
    return [-numpy.inf if val is None else val for val in seq]
