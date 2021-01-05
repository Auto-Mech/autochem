""" stereo graph library
"""
import itertools
import more_itertools as mit
import numpy
from qcelemental import constants as qcc
from automol import dict_
import automol.geom
import automol.zmatrix
import automol.create.geom
from automol import cart
from automol import formula
from automol.graph._res import resonance_dominant_atom_hybridizations
from automol.graph._ring import rings_atom_keys
from automol.graph._stereo import has_stereo
from automol.graph._stereo import stereogenic_atom_keys
from automol.graph._stereo import stereogenic_bond_keys
from automol.graph._stereo import stereo_sorted_atom_neighbor_keys
from automol.graph._graph import atom_keys
from automol.graph._graph import atom_symbols
from automol.graph._graph import add_bonded_atom
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import explicit
from automol.graph._graph import branch
from automol.graph._graph import longest_chain
from automol.graph._graph import atom_stereo_parities
from automol.graph._graph import bond_stereo_parities
from automol.graph._graph import without_stereo_parities
from automol.graph._graph import set_atom_stereo_parities
from automol.graph._graph import set_bond_stereo_parities
from automol.graph._graph import connected_components


def heuristic_geometry(gra):
    """ stereo-specific coordinates for a molecular geometry

    (need not be connected)
    """
    assert gra == explicit(gra)
    gra_iter = iter(connected_components(gra))
    gra_ = next(gra_iter)
    geo_, geo_idx_dct_ = _connected_heuristic_geometry(gra_)

    geo = geo_
    geo_idx_dct = geo_idx_dct_
    for gra_ in gra_iter:
        geo_, geo_idx_dct_ = _connected_heuristic_geometry(gra_)

        natms = automol.geom.count(geo)
        geo_idx_dct_ = dict_.transform_values(geo_idx_dct_, (natms).__add__)
        geo_idx_dct.update(geo_idx_dct_)

        geo = automol.geom.join(geo, geo_)

    return geo, geo_idx_dct


def _connected_heuristic_geometry(gra):
    """ stereo-specific coordinates for a connected molecular geometry
    """
    assert gra == explicit(gra)

    atm_keys = sorted(atom_keys(gra))

    zma, zma_key_dct = connected_heuristic_zmatrix(gra)

    geo = automol.zmatrix.geometry(zma)
    idxs = dict_.values_by_key(zma_key_dct, atm_keys)
    geo = automol.geom.from_subset(geo, idxs)
    geo_idx_dct = {atm_key: idx for idx, atm_key in enumerate(atm_keys)}

    return geo, geo_idx_dct


def connected_heuristic_zmatrix(gra):
    """ stereo-specific coordinates for a connected molecular graph

    (currently unable to handle rings -- fix that)
    """
    assert gra == explicit(gra)

    atm_sym_dct = atom_symbols(gra)

    # this will contain triplets of adjacent atoms from which to continue
    # filling out the z-matrix, after it has been started
    triplets = []

    # 1. start the z-matrix and set the lists of triplets
    rng_atm_keys_lst = rings_atom_keys(gra)
    if not rng_atm_keys_lst:
        # find the first heavy atom in the longest chain (if there isn't one,
        # we are dealing with atomic or molecular hydrogen, which will be
        # captured by the last two cases)
        chain = longest_chain(gra)

        if atm_sym_dct[chain[0]] != 'H':
            chain = list(reversed(chain))

        if len(chain) > 1:
            atm_key = chain[1]
        else:
            atm_key = chain[0]

        # determine the z-matrix of the starting atom and its neighbors
        zma, zma_key_dct, dummy_atm_key, gra = _start_zmatrix_from_atom(
            gra, atm_key)

        # since this is the first heavy atom in the longest chain, we only need
        # to follow one branch from this atom to complete the z-matrix; this
        # will be the branch extending toward the next heavy atom in the chai
        if len(chain) > 3:
            atm1_key, atm2_key, atm3_key = chain[:3]

            # if we inserted a dummy atom on the starting geometry, we should
            # use that as atom 1 in the triplet, rather than
            if dummy_atm_key is not None:
                atm1_key = dummy_atm_key

            triplets = [(atm1_key, atm2_key, atm3_key)]
    elif len(rng_atm_keys_lst) == 1:
        rng_atm_keys, = rng_atm_keys_lst

        zma, zma_key_dct = _start_zmatrix_from_ring(gra, rng_atm_keys)

        triplets += list(mit.windowed(rng_atm_keys[-2:] + rng_atm_keys, 3))
    else:
        # currently, multiple rings are not implemented
        raise NotImplementedError

    # 2. complete the z-matrix by looping over triplets
    for atm1_key, atm2_key, atm3_key in triplets:
        zma, zma_key_dct, gra = _complete_zmatrix_for_branch(
            gra, atm1_key, atm2_key, atm3_key, zma, zma_key_dct)

    # 3. convert to Cartesian geometry for stereo correction
    geo = automol.zmatrix.geometry(zma)
    geo_idx_dct = zma_key_dct
    geo = _stereo_corrected_geometry(gra, geo, geo_idx_dct)

    # 4. convert back to z-matrix, keeping the original z-matrix structure
    vma = automol.zmatrix.var_(zma)
    zma = automol.zmatrix.from_geometry(vma, geo)

    return zma, zma_key_dct


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
    atm_ngb_keys_dct = atom_neighbor_keys(gra)
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    # sort the neighbor keys by stereo priority
    atm_ngb_keys = stereo_sorted_atom_neighbor_keys(
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
    atm_ngb_keys_dct = atom_neighbor_keys(gra)
    atm1_ngb_keys = atm_ngb_keys_dct[atm1_key] - {atm2_key}
    atm2_ngb_keys = atm_ngb_keys_dct[atm2_key] - {atm1_key}

    atm1_ngb_keys = stereo_sorted_atom_neighbor_keys(
        gra, atm1_key, atm1_ngb_keys)
    atm2_ngb_keys = stereo_sorted_atom_neighbor_keys(
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
    atm_ngb_keys_dct = atom_neighbor_keys(gra)

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
            rot_axis = cart.vec.unit_bisector(
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


# connected graph heuristic z-matrix code
TET_ANG = numpy.arccos(-1/3)
TRI_ANG = 2 * numpy.pi / 3
LIN_ANG = numpy.pi
RIT_ANG = numpy.pi / 2


def _start_zmatrix_from_ring(gra, rng_atm_keys):
    """ generates a z-matrix for a ring
    """
    # the key dictionary can be constructed immediately
    zma_key_dct = {
        atm_key: zma_key for zma_key, atm_key in enumerate(rng_atm_keys)}

    # now, build the z-matrix
    natms = len(rng_atm_keys)
    atm_sym_dct = atom_symbols(gra)

    dist_val = 1.5 * qcc.conversion_factor('angstrom', 'bohr')
    ang_val = (natms - 2.) * numpy.pi / natms
    dih_val = 0.

    # 1. construct the z-matrix for a 3-atom system
    key_mat = [[None, None, None],
               [0, None, None],
               [1, 0, None]]

    name_mat = [[None, None, None],
                ['R1', None, None],
                ['R2', 'A2', None]]

    syms = dict_.values_by_key(atm_sym_dct, rng_atm_keys[:3])

    val_dct = {'R1': dist_val, 'R2': dist_val, 'A2': ang_val}

    zma = automol.zmatrix.from_data(syms, key_mat, name_mat, val_dct)

    # 2. append z-matrix rows for the remaining atoms
    for row, rng_atm_key in enumerate(rng_atm_keys):
        if row > 2:
            sym = atm_sym_dct[rng_atm_key]
            dist_name = automol.zmatrix.new_distance_name(zma)
            ang_name = automol.zmatrix.new_central_angle_name(zma)
            dih_name = automol.zmatrix.new_dihedral_angle_name(zma)

            zma = automol.zmatrix.append(
                zma, sym,
                [row-1, row-2, row-3],
                [dist_name, ang_name, dih_name],
                [dist_val, ang_val, dih_val])

    return zma, zma_key_dct


def _start_zmatrix_from_atom(gra, atm_key):
    """ generates a z-matrix for a single atom and its neighbors

    returns a z-matrix and a dictionary mapping atom keys to rows of the
    z-matrix
    """
    atm_ngb_keys_dct = atom_neighbor_keys(gra)

    dummy_atm_key = None

    atm_sym_dct = atom_symbols(gra)

    # sort hydrogens to be first
    atm_ngb_keys = sorted(atm_ngb_keys_dct[atm_key])
    atm_ngb_syms = list(map(atm_sym_dct.__getitem__, atm_ngb_keys))
    srt = formula.argsort_symbols(atm_ngb_syms, syms_first=('H', 'C'))
    atm_ngb_keys = list(map(atm_ngb_keys.__getitem__, srt))

    atm_keys = [atm_key] + atm_ngb_keys

    zma_key_dct = {
        atm_key: zma_key for zma_key, atm_key in enumerate(atm_keys)}

    syms = list(map(atm_sym_dct.__getitem__, atm_keys))

    natms = len(atm_keys)

    key_mat = [[None, None, None],
               [0, None, None],
               [0, 1, None],
               [0, 1, 2],
               [0, 1, 2]][:natms]

    name_mat = [[None, None, None],
                ['R1', None, None],
                ['R2', 'A2', None],
                ['R3', 'A3', 'D3'],
                ['R4', 'A4', 'D4']][:natms]

    atm_hyb = resonance_dominant_atom_hybridizations(gra)[atm_key]

    # z-matrix coordinate values
    val_dct = {}

    # determine bond distances
    for dist_name, atm_ngb_key in zip(['R1', 'R2', 'R3', 'R4'], atm_ngb_keys):
        dist_val = _bond_distance(gra, atm_key, atm_ngb_key)
        val_dct[dist_name] = dist_val

    # determine bond angles and dihedral angles
    if atm_hyb == 3:
        # for sp3 atoms, use z-matrix for a tetrahedral structure
        val_dct.update({'A2': TET_ANG, 'A3': TET_ANG, 'A4': TET_ANG,
                        'D3': +TRI_ANG, 'D4': -TRI_ANG})
    elif atm_hyb == 2:
        # for sp2 atoms, use z-matrix for a trigonal planar structure
        val_dct.update({'A2': TRI_ANG, 'A3': TRI_ANG, 'D3': LIN_ANG})
    elif atm_hyb == 1 and natms > 2:
        # for sp1 atoms, uze z-matrix for a linear structure

        # if there's only one neighbor, we don't need to set any angles at all
        assert natms == 3

        # insert a dummy atom
        syms.insert(2, 'X')
        key_mat.append([0, 2, 1])
        name_mat.append(['R3', 'A3', 'D3'])

        # note: we need to insert the dummy atom bond distance at R2 and shift
        # over the current value
        val_dct.update({'R2': 1.0, 'R3': val_dct['R2'],
                        'A2': RIT_ANG, 'A3': RIT_ANG, 'D3': LIN_ANG})

        gra, dummy_atm_key = add_bonded_atom(gra, 'X', atm_key, bnd_ord=0)

        zma_key_dct[dummy_atm_key] = 2
        zma_key_dct[atm_keys[2]] = 3

    vma = automol.vmatrix.from_data(syms, key_mat, name_mat)

    names = automol.vmatrix.names(vma)
    val_dct = {name: val_dct[name] for name in names}

    zma = automol.zmatrix.from_data(syms, key_mat, name_mat, val_dct)

    gra_with_dummies = gra
    return zma, zma_key_dct, dummy_atm_key, gra_with_dummies


def _complete_zmatrix_for_branch(gra, atm1_key, atm2_key, atm3_key, zma,
                                 zma_key_dct):
    """ core function for generating geometries; starting from three
    neighboring atoms at the end of the z-matrix, fills out the z-matrix for
    the rest of the branch extending out from the third atom

    returns a z-matrix and a dictionary mapping atom keys to rows of the
    z-matrix
    """
    atm_sym_dct = atom_symbols(gra)
    atm_ngb_keys_dct = atom_neighbor_keys(gra)

    atm1_zma_key = zma_key_dct[atm1_key]
    atm2_zma_key = zma_key_dct[atm2_key]
    atm3_zma_key = zma_key_dct[atm3_key]

    atm4_keys = atm_ngb_keys_dct[atm3_key] - {atm2_key}

    # first handle the linear bond case, inserting a dummy atom
    atm_hyb = resonance_dominant_atom_hybridizations(gra)[atm3_key]
    if atm_hyb == 1 and atm4_keys:
        assert len(atm4_keys) == 1
        atm4_key, = atm4_keys

        # first, insert the dummy atom
        dist4_name = automol.zmatrix.new_distance_name(zma)
        dist4_val = 1.

        ang4_name = automol.zmatrix.new_central_angle_name(zma)
        ang4_val = RIT_ANG

        dih4_name = automol.zmatrix.new_dihedral_angle_name(zma)
        dih4_val = 0.

        gra, dummy_atm_key = add_bonded_atom(gra, 'X', atm3_key, bnd_ord=0)
        dummy_atm_zma_key = automol.zmatrix.count(zma)

        zma_key_dct[dummy_atm_key] = dummy_atm_zma_key

        zma = automol.zmatrix.append(
            zma, 'X',
            [atm3_zma_key, atm2_zma_key, atm1_zma_key],
            [dist4_name, ang4_name, dih4_name],
            [dist4_val, ang4_val, dih4_val])

        # shift the keys to include the dummy atom
        atm1_key, atm2_key = atm2_key, dummy_atm_key
        atm1_zma_key, atm2_zma_key = atm2_zma_key, dummy_atm_zma_key

        atm4_sym = atm_sym_dct[atm4_key]

        dist4_name = automol.zmatrix.new_distance_name(zma)
        dist4_val = _bond_distance(gra, atm3_key, atm4_key)

        ang4_name = automol.zmatrix.new_central_angle_name(zma)
        ang4_val = RIT_ANG

        dih4_name = automol.zmatrix.new_dihedral_angle_name(zma)
        dih4_val = LIN_ANG

        zma_key_dct[atm4_key] = automol.zmatrix.count(zma)

        zma = automol.zmatrix.append(
            zma, atm4_sym,
            [atm3_zma_key, atm2_zma_key, atm1_zma_key],
            [dist4_name, ang4_name, dih4_name],
            [dist4_val, ang4_val, dih4_val])

        geo = automol.zmatrix.geometry(zma)

        # recursion 1
        zma, zma_key_dct, gra = _complete_zmatrix_for_branch(
            gra, atm2_key, atm3_key, atm4_key, zma, zma_key_dct)

    # from here on, we can assume a non-linear bond
    else:
        dih_incr = _dihedral_increment(gra, atm1_key, atm2_key, atm3_key)

        fixed_atm4_keys = set(atm4_keys) & set(zma_key_dct)
        if fixed_atm4_keys:
            geo = automol.zmatrix.geometry(zma)
            atm4_zma_keys = list(map(zma_key_dct.__getitem__, fixed_atm4_keys))
            dih4_vals = [
                automol.geom.dihedral_angle(geo, atm1_zma_key, atm2_zma_key,
                                            atm3_zma_key, atm4_zma_key)
                for atm4_zma_key in atm4_zma_keys]
            dih4_val = max(dih4_vals)
        else:
            dih4_val = numpy.pi - dih_incr

        # subtract off the fixed atom 4 keys in case we're dealing with a ring
        for atm4_key in atm4_keys - fixed_atm4_keys:
            atm4_sym = atm_sym_dct[atm4_key]

            dist4_name = automol.zmatrix.new_distance_name(zma)
            dist4_val = _bond_distance(gra, atm3_key, atm4_key)

            ang4_name = automol.zmatrix.new_central_angle_name(zma)
            ang4_val = _bond_angle(gra, atm2_key, atm3_key, atm4_key)

            dih4_name = automol.zmatrix.new_dihedral_angle_name(zma)

            dih4_val += dih_incr
            dih4_val = numpy.mod(dih4_val, 2*numpy.pi)

            zma_key_dct[atm4_key] = automol.zmatrix.count(zma)

            zma = automol.zmatrix.append(
                zma, atm4_sym,
                [atm3_zma_key, atm2_zma_key, atm1_zma_key],
                [dist4_name, ang4_name, dih4_name],
                [dist4_val, ang4_val, dih4_val])

            # recursion 2
            zma, zma_key_dct, gra = _complete_zmatrix_for_branch(
                gra, atm2_key, atm3_key, atm4_key, zma, zma_key_dct)

    gra_with_dummies = gra

    return zma, zma_key_dct, gra_with_dummies


def _bond_distance(gra, atm1_key, atm2_key, check=True):
    """ predicted bond distance

    (currently crude, but could easily be made more sophisticated
    """
    atm_sym_dct = atom_symbols(gra)

    if check:
        assert atm2_key in atom_neighbor_keys(gra)[atm1_key]

    atm1_sym = atm_sym_dct[atm1_key]
    atm2_sym = atm_sym_dct[atm2_key]

    if atm1_sym == 'H' or atm2_sym == 'H':
        dist = 1.1 * qcc.conversion_factor('angstrom', 'bohr')
    else:
        dist = 1.5 * qcc.conversion_factor('angstrom', 'bohr')

    return dist


def _bond_angle(sgr, atm1_key, atm2_key, atm3_key, check=True):
    """ predict the bond angles an atom makes with its neighbors
    """
    if check:
        atm_ngb_keys_dct = atom_neighbor_keys(sgr)
        atm2_ngb_keys = atm_ngb_keys_dct[atm2_key]
        assert {atm1_key, atm3_key} <= atm2_ngb_keys

    atm_hyb_dct = resonance_dominant_atom_hybridizations(sgr)
    atm2_hyb = atm_hyb_dct[atm2_key]

    if atm2_hyb == 3:
        ang_deg = 109.5
    elif atm2_hyb == 2:
        ang_deg = 120.0
    else:
        assert atm2_hyb == 1
        ang_deg = 180.0

    ang = ang_deg * qcc.conversion_factor('degree', 'radian')
    return ang


def _dihedral_increment(sgr, atm1_key, atm2_key, atm3_key, check=False):
    """ predict dihedral increment for atoms attached to `atm3_key`
    """
    if check:
        atm_ngb_keys_dct = atom_neighbor_keys(sgr)
        atm2_ngb_keys = atm_ngb_keys_dct[atm2_key]
        atm3_ngb_keys = atm_ngb_keys_dct[atm3_key]
        assert atm2_key in atm3_ngb_keys
        assert atm1_key in atm2_ngb_keys - {atm3_key}

    atm_hyb_dct = resonance_dominant_atom_hybridizations(sgr)
    atm3_hyb = atm_hyb_dct[atm3_key]
    dih_incr = 2 * numpy.pi / atm3_hyb if not atm3_hyb == 0 else 0.
    return dih_incr


if __name__ == '__main__':
    import automol.graph

    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 1, None)},
    #        {frozenset({4, 5}): (1, False), frozenset({0, 2}): (1, None),
    #         frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
    #         frozenset({1, 3}): (1, None)})
    # GRA2 = ({0: ('C', 1, None), 1: ('C', 3, None), 2: ('C', 0, None),
    #          3: ('C', 1, False), 4: ('O', 1, None)},
    #         {frozenset({3, 4}): (1, None), frozenset({0, 2}): (1, None),
    #          frozenset({1, 3}): (1, None), frozenset({2, 3}): (1, None)})
    # GRA2 = automol.graph.transform_keys(GRA2, lambda x: x+6)
    # GRA = automol.graph.union(GRA, GRA2)

    # GRA = ({0: ('Ne', 0, None)}, {})
    # GRA = ({0: ('C', 1, None), 1: ('C', 0, None), 2: ('O', 1, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})
    # GRA = ({0: ('C', 1, False), 1: ('C', 1, False), 2: ('F', 0, None),
    #         3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
    #         frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
    #         frozenset({1, 5}): (1, None)})
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
    #        {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
    #         frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
    #         frozenset({2, 5}): (1, True)})
    # GRA = ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 1, False),
    #         3: ('O', 1, None), 4: ('O', 0, None), 5: ('O', 0, None),
    #         6: ('O', 0, None)},
    #        {frozenset({4, 6}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({2, 6}): (1, None), frozenset({0, 2}): (1, None),
    #         frozenset({1, 5}): (1, None), frozenset({3, 5}): (1, None)})
    # GRA = ({0: ('C', 2, None), 1: ('C', 2, None), 2: ('C', 1, True),
    #         3: ('O', 1, None), 4: ('O', 0, None), 5: ('O', 0, None)},
    #        {frozenset({1, 2}): (1, None), frozenset({0, 2}): (1, None),
    #         frozenset({0, 4}): (1, None), frozenset({2, 4}): (1, None),
    #         frozenset({1, 5}): (1, None), frozenset({3, 5}): (1, None)})
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None)},
    #        {frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({0, 2}): (1, None), frozenset({2, 4}): (1, None),
    #         frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})
    # GRA = ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 1, None),
    #         3: ('C', 1, True), 4: ('O', 0, None), 5: ('O', 1, None),
    #         6: ('O', 0, None)},
    #        {frozenset({2, 3}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({3, 6}): (1, None), frozenset({2, 4}): (1, None),
    #         frozenset({5, 6}): (1, None), frozenset({1, 3}): (1, None)})
    # GRA = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
    #         3: ('C', 1, True), 4: ('C', 1, False), 5: ('O', 0, None),
    #         6: ('O', 1, None), 7: ('O', 0, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({1, 4}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({0, 3}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({2, 5}): (1, None),
    #         frozenset({4, 7}): (1, None)})

    # GRA = ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
    #         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, True)},
    #        {frozenset({1, 4}): (1, None), frozenset({0, 3}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({5, 6}): (1, False),
    #         frozenset({3, 5}): (1, False), frozenset({2, 7}): (1, None),
    #         frozenset({4, 7}): (1, None)})

    GRA = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, False),
            3: ('C', 1, True), 4: ('C', 1, True), 5: ('O', 1, None),
            6: ('O', 0, None), 7: ('O', 0, None)},
           {frozenset({3, 4}): (1, None), frozenset({2, 6}): (1, None),
            frozenset({0, 2}): (1, None), frozenset({3, 6}): (1, None),
            frozenset({2, 4}): (1, None), frozenset({1, 3}): (1, None),
            frozenset({5, 7}): (1, None), frozenset({4, 7}): (1, None)})

    SGR = explicit(GRA)
    GEO, GEO_IDX_DCT = heuristic_geometry(SGR)

    # GEO = (
    #     ('C', (-4.3870588134, -1.233231672517, 0.143749726309016)),
    #     ('C', (3.430304171771, -2.162836645393, 0.06129774977456508)),
    #     ('C', (-2.228277354885, 0.1940343942502, -1.0788747507898575)),
    #     ('C', (1.642516616594, -0.2666792157215, -1.1217150146524835)),
    #     ('C', (-0.128929182575, 1.167922277159, 0.6891116967734786)),
    #     ('O', (2.44416862025, 4.3088944278, 2.2811617948010614)),
    #     ('O', (-0.4839163664245, -1.530973627393, -2.314392676447687)),
    #     ('O', (0.1546996690877, 3.88389710099, 0.8186943317590577)),
    #     ('H', (-3.69835752847, -2.73865502820, 1.385388988006255)),
    #     ('H', (-5.55611791656, 0.04579783881378, 1.2732325623382277)),
    #     ('H', (-5.59710374244, -2.09637223693, -1.2950602194746856)),
    #     ('H', (2.43356758828, -3.46358524781, 1.3249230383597128)),
    #     ('H', (4.89255603965, -1.19599087753, 1.1592100844834574)),
    #     ('H', (4.37095524985, -3.286477042474, -1.3988006492215568)),
    #     ('H', (-2.909935258917, 1.611421581693, -2.430728150238205)),
    #     ('H', (2.608452006545, 0.954321649398, -2.49176226302396)),
    #     ('H', (-0.2313534706629, 0.353123800972, 2.5932645521380446)),
    #     ('H', (3.243829672371, 5.45538852286, 1.0816218339892072)))

    XGR = set_stereo_from_geometry(SGR, GEO, geo_idx_dct=GEO_IDX_DCT)

    print(automol.geom.string(GEO))
    REF_ATM_STE_PAR_DCT = dict_.filter_by_value(
        atom_stereo_parities(GRA), lambda x: x is not None)
    ATM_STE_PAR_DCT = dict_.filter_by_value(
        atom_stereo_parities(XGR), lambda x: x is not None)
    print(REF_ATM_STE_PAR_DCT)
    print(ATM_STE_PAR_DCT)
