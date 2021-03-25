""" code for setting graph stereo from a geometry
"""
import itertools
import numpy
from automol.util import dict_
import automol.geom
import automol.create.geom
from automol import util
from automol.graph._ring import rings_atom_keys
from automol.graph._stereo import has_stereo
from automol.graph._stereo import stereogenic_atom_keys
from automol.graph._stereo import stereogenic_bond_keys
from automol.graph._stereo import atom_stereo_sorted_neighbor_atom_keys
from automol.graph._graph import atom_keys
from automol.graph._graph import atoms_neighbor_atom_keys
from automol.graph._graph import explicit
from automol.graph._graph import branch
from automol.graph._graph import atom_stereo_parities
from automol.graph._graph import bond_stereo_parities
from automol.graph._graph import without_stereo_parities
from automol.graph._graph import set_atom_stereo_parities
from automol.graph._graph import set_bond_stereo_parities


# stereo setting code
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


def _set_atom_stereo_from_geometry(gra, atm_keys, geo, geo_idx_dct):
    assert gra == explicit(gra)

    atm_pars = [
        _atom_stereo_parity_from_geometry(gra, atm_key, geo, geo_idx_dct)
        for atm_key in atm_keys]
    gra = set_atom_stereo_parities(gra, dict(zip(atm_keys, atm_pars)))
    return gra


def _set_bond_stereo_from_geometry(gra, bnd_keys, geo, geo_idx_dct):
    assert gra == explicit(gra)

    bnd_pars = [
        _bond_stereo_parity_from_geometry(gra, bnd_key, geo, geo_idx_dct)
        for bnd_key in bnd_keys]
    gra = set_bond_stereo_parities(gra, dict(zip(bnd_keys, bnd_pars)))
    return gra


# stereo parity evaluation code
def _atom_stereo_parity_from_geometry(gra, atm_key, geo, geo_idx_dct):
    """ get the current stereo parity of an atom from its geometry
    """
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    # sort the neighbor keys by stereo priority
    atm_ngb_keys = atom_stereo_sorted_neighbor_atom_keys(
        gra, atm_key, atm_ngb_keys)

    # determine the parity based on the coordinates
    xyzs = automol.geom.coordinates(geo)
    atm_ngb_idxs = dict_.values_by_key(geo_idx_dct, atm_ngb_keys)
    atm_ngb_xyzs = [xyzs[idx] for idx in atm_ngb_idxs]
    det_mat = numpy.ones((4, 4))
    det_mat[:, :3] = atm_ngb_xyzs
    det_val = numpy.linalg.det(det_mat)
    assert det_val != 0.  # for now, assume no four-atom planes
    par = det_val > 0.
    return par


def _bond_stereo_parity_from_geometry(gra, bnd_key, geo, geo_idx_dct):
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
    xyzs = automol.geom.coordinates(geo)
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


# stereo correction code
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
        curr_par = _atom_stereo_parity_from_geometry(
            gra, atm_key, geo, geo_idx_dct)

        if curr_par != par:
            atm_ngb_keys = atm_ngb_keys_dct[atm_key]
            # for now, we simply exclude rings from the pivot keys
            # (will not work for stereo atom at the intersection of two rings)
            atm_piv_keys = list(atm_ngb_keys - ring_atm_keys)[:2]
            assert len(atm_piv_keys) == 2
            atm3_key, atm4_key = atm_piv_keys

            # get coordinates
            xyzs = automol.geom.coordinates(geo)
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

        assert _atom_stereo_parity_from_geometry(
            gra, atm_key, geo, geo_idx_dct) == par
        gra = set_atom_stereo_parities(gra, {atm_key: par})

    return geo, gra


def _bond_stereo_corrected_geometry(gra, bnd_ste_par_dct, geo, geo_idx_dct):
    """ correct the bond stereo parities of a geometry, for a subset of bonds
    """
    bnd_keys = list(bnd_ste_par_dct.keys())
    for bnd_key in bnd_keys:
        par = bnd_ste_par_dct[bnd_key]
        curr_par = _bond_stereo_parity_from_geometry(
            gra, bnd_key, geo, geo_idx_dct)

        if curr_par != par:
            xyzs = automol.geom.coordinates(geo)

            atm1_key, atm2_key = bnd_key
            atm1_xyz = xyzs[geo_idx_dct[atm1_key]]
            atm2_xyz = xyzs[geo_idx_dct[atm2_key]]

            rot_axis = numpy.subtract(atm2_xyz, atm1_xyz)

            rot_atm_keys = atom_keys(
                branch(gra, atm1_key, {atm1_key, atm2_key}))

            rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

            geo = automol.geom.rotate(
                geo, rot_axis, numpy.pi, orig_xyz=atm1_xyz, idxs=rot_idxs)

        assert _bond_stereo_parity_from_geometry(
            gra, bnd_key, geo, geo_idx_dct) == par
        gra = set_bond_stereo_parities(gra, {bnd_key: par})

    return geo, gra
