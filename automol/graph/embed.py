""" geometry embedding using the distance geometry algorithm

This algorithm generates approximate geometries by a randomized guess at the
distance matrix, within heuristic distance bounds, which is then converted to
an approximate geometry that satifies these distances.

Blaney, J. M.; Dixon, J. S. “Distance Geometry in Molecular Modeling”.
Reviews in Computational Chemistry; VCH: New York, 1994.

The steps in the algorithm are as follows:

    1. Generate distance bounds matrix B. (Here, B is replaced with L and U).
    2. Narrow the bounds in B by triangle smoothing.
    3. Generate a distance matrix D by uniform sampling within the bounds.
    4. Generate the metric matrix G (matrix of position vector dot products).
    5. Diagonalize G and determine the principal components.
    6. The three largest eigenvectors and eigenvalues of G can be used to
    generate x, y, z coordinates for the molecule which approximately
    correspond to the distance matrix D.
    7. Do error refinement to clean up the structure and enforce correct
    chirality.

Step 1 is performed by the function distance_bounds_matrices() below.

The function sample_raw_distance_geometry() below returns the result of step
1-6. The actual work for these steps is handled in a separate module,
automol.embed.
"""

import numbers
import itertools
import more_itertools as mit
import numpy
from phydat import phycon
from phydat import ptab
from automol import error
from automol import embed
from automol.util import dict_
import automol.geom.base
from automol.graph.base import atom_keys
from automol.graph.base import bond_keys
from automol.graph.base import atom_stereo_keys
from automol.graph.base import bond_stereo_keys
from automol.graph.base import atom_neighbor_atom_key
from automol.graph.base import atom_shortest_paths
from automol.graph.base import atom_symbols
from automol.graph.base import atom_stereo_parities
from automol.graph.base import bond_stereo_parities
from automol.graph.base import string
from automol.graph.base import explicit
from automol.graph.base import atom_neighborhoods
from automol.graph.base import bond_neighborhoods
from automol.graph.base import atoms_neighbor_atom_keys
from automol.graph.base import subgraph
from automol.graph.base import atom_stereo_sorted_neighbor_atom_keys
from automol.graph.base import resonance_dominant_atom_hybridizations
from automol.graph.base import atom_van_der_waals_radius
from automol.graph.base import sp2_bond_keys
from automol.graph.base import rings_atom_keys
from automol.graph.base import atom_stereo_parity_from_geometry
from automol.graph.base import bond_stereo_parity_from_geometry


# bond distances
XY_DIST = 1.5       # angstroms
XH_DIST = 1.1       # angstroms

# bond angles
TET_ANG = 109.4712  # degrees
TRI_ANG = 120.      # degrees
LIN_ANG = 180.      # degrees


# # geometry embedding functions
def geometry(gra, keys=None, ntries=5, max_dist_err=0.2):
    """ sample a qualitatively-correct stereo geometry

    :param gra: the graph, which may or may not have stereo
    :param keys: graph keys, in the order in which they should appear in the
        geometry
    :param ntries: number of tries for finding a valid geometry
    :param max_dist_err: maximum distance error convergence threshold

    Qualitatively-correct means it has the right connectivity and the right
    stero parities, but its bond lengths and bond angles may not be
    quantitatively realistic
    """
    assert gra == explicit(gra), (
        "Graph => geometry conversion requires explicit hydrogens!\n"
        "Use automol.graph.explicit() to convert to an explicit graph.")

    # 0. Get keys and symbols
    symb_dct = atom_symbols(gra)

    keys = sorted(atom_keys(gra)) if keys is None else keys
    symbs = tuple(map(symb_dct.__getitem__, keys))

    # 1. Generate bounds matrices
    lmat, umat = distance_bounds_matrices(gra, keys)
    chi_dct = chirality_constraint_bounds(gra, keys)
    pla_dct = planarity_constraint_bounds(gra, keys)
    conv1_ = qualitative_convergence_checker_(gra, keys)
    conv2_ = embed.distance_convergence_checker_(lmat, umat, max_dist_err)

    def conv_(xmat, err, grad):
        return conv1_(xmat, err, grad) & conv2_(xmat, err, grad)

    # 2. Generate coordinates with correct stereo, trying a few times
    for _ in range(ntries):
        xmat = embed.sample_raw_distance_coordinates(lmat, umat, dim4=True)
        xmat, conv = embed.cleaned_up_coordinates(
            xmat, lmat, umat, pla_dct=pla_dct, chi_dct=chi_dct, conv_=conv_)
        if conv:
            break

    if not conv:
        raise error.FailedGeometryGenerationError(f'Bad gra {string(gra)}')

    # 3. Generate a geometry data structure from the coordinates
    xyzs = xmat[:, :3]
    geo = automol.geom.base.from_data(symbs, xyzs, angstrom=True)

    return geo


def fake_stereo_geometry(gra, ntries=5, max_dist_err=0.5):
    """ generate a fake stereo geometry
    """
    # determine stereo "groups" with geometrically interdependent chirality
    atm_ngbs_dct = atom_neighborhoods(gra)
    bnd_ngbs_dct = bond_neighborhoods(gra)
    atm_ste_keys = atom_stereo_keys(gra)
    bnd_ste_keys = bond_stereo_keys(gra)
    atm_ste_groups = list(
        map(atom_keys, map(atm_ngbs_dct.__getitem__, atm_ste_keys)))
    bnd_ste_groups = list(
        map(atom_keys, map(bnd_ngbs_dct.__getitem__, bnd_ste_keys)))

    ste_groups = _aggregate_connected_groups(atm_ste_groups + bnd_ste_groups)

    ste_groups = list(map(sorted, ste_groups))

    natms = 0
    geo_idx_dct = {}
    geo = ()
    for group in ste_groups:
        group_geo = geometry(
            gra, keys=group, ntries=ntries,
            max_dist_err=max_dist_err)
        group_natms = len(group)

        idxs = list(range(natms, natms+group_natms))
        geo_idx_dct.update(dict(zip(group, idxs)))

        natms += group_natms
        geo = automol.geom.base.join(geo, group_geo)

    return geo, geo_idx_dct


def _aggregate_connected_groups(keys_lst):
    """ from a list of groups, join the ones with common atoms
    """
    groups = []
    for keys in keys_lst:
        idxs = list(idx for idx, group in enumerate(groups)
                    if group & keys)
        if not idxs:
            groups.append(keys)
        else:
            groups[idxs[0]] |= keys
            for idx in idxs[1:]:
                groups[idxs[0]] |= groups[idx]
                groups.pop(idx)

    groups = tuple(groups)

    return groups


# # convergence checking
def qualitative_convergence_checker_(gra, keys, rqq_bond_max=1.8,
                                     rqh_bond_max=1.3, rhh_bond_max=1.1,
                                     bond_nobond_diff=0.3):
    """ a convergence checker for error minimization, checking that the
    geometry is qualitatively correct (correct connectivity and stereo)
    """
    symb_dct = atom_symbols(gra)
    pairs = set(map(frozenset, itertools.combinations(keys, 2)))

    bnd_keys = pairs & bond_keys(gra)
    nob_keys = pairs - bond_keys(gra)

    nob_symbs = tuple(tuple(map(symb_dct.__getitem__, nob_key))
                      for nob_key in nob_keys)
    bnd_symbs = tuple(tuple(map(symb_dct.__getitem__, bnd_key))
                      for bnd_key in bnd_keys)
    nob_idxs = tuple(tuple(map(keys.index, nob_key)) for nob_key in nob_keys)
    bnd_idxs = tuple(tuple(map(keys.index, bnd_key)) for bnd_key in bnd_keys)

    bnd_udists = tuple((rqq_bond_max if 'H' not in symb else
                        rhh_bond_max if set(symb) == {'H'} else
                        rqh_bond_max) for symb in bnd_symbs)

    diff = bond_nobond_diff
    nob_ldists = tuple((rqq_bond_max+diff if 'H' not in symb else
                        rhh_bond_max+diff if set(symb) == {'H'} else
                        rqh_bond_max+diff) for symb in nob_symbs)

    bnd_idxs += tuple(map(tuple, map(reversed, bnd_idxs)))
    bnd_idx_vecs = tuple(map(list, zip(*bnd_idxs)))
    bnd_udists *= 2

    nob_idxs += tuple(map(tuple, map(reversed, nob_idxs)))
    nob_idx_vecs = tuple(map(list, zip(*nob_idxs)))
    nob_ldists *= 2

    symbs = tuple(map(symb_dct.__getitem__, keys))
    geo_idx_dct = dict(map(reversed, enumerate(keys)))
    atm_ste_keys = atom_stereo_keys(gra) & set(keys)
    bnd_ste_keys = bond_stereo_keys(gra) & bnd_keys
    atm_ste_par_dct = atom_stereo_parities(gra)
    bnd_ste_par_dct = bond_stereo_parities(gra)

    def _is_converged(xmat, err, grad):
        assert err and numpy.any(grad)
        xyzs = xmat[:, :3]
        dmat = embed.distance_matrix_from_coordinates(xyzs)

        # check for correct connectivity
        connectivity_check = (
            (numpy.all(dmat[bnd_idx_vecs] < bnd_udists)
             if bnd_udists else True) and
            (numpy.all(dmat[nob_idx_vecs] > nob_ldists)
             if nob_ldists else True))

        # check for correct stereo parities
        geo = automol.geom.base.from_data(symbs, xyzs, angstrom=True)
        atom_stereo_check = all(
            (atom_stereo_parity_from_geometry(gra, atm_key, geo, geo_idx_dct)
             == atm_ste_par_dct[atm_key])
            for atm_key in atm_ste_keys)

        bond_stereo_check = all(
            (bond_stereo_parity_from_geometry(gra, bnd_key, geo, geo_idx_dct)
             == bnd_ste_par_dct[bnd_key])
            for bnd_key in bnd_ste_keys)

        return connectivity_check and atom_stereo_check and bond_stereo_check

    return _is_converged


# # bounds matrices
def distance_bounds_matrices(gra, keys, sp_dct=None):
    """ initial distance bounds matrices

    :param gra: molecular graph
    :param keys: atom keys specifying the order of indices in the matrix
    :param sp_dct: a 2d dictionary giving the shortest path between any pair of
        atoms in the graph
    """
    assert set(keys) <= set(atom_keys(gra))

    sub_gra = subgraph(gra, keys, stereo=True)
    sp_dct = atom_shortest_paths(sub_gra) if sp_dct is None else sp_dct

    bounds_ = path_distance_bounds_(gra)

    natms = len(keys)
    umat = numpy.zeros((natms, natms))
    lmat = numpy.zeros((natms, natms))
    for (idx1, key1), (idx2, key2) in itertools.combinations(
            enumerate(keys), 2):
        if key2 in sp_dct[key1]:
            path = sp_dct[key1][key2]
            ldist, udist = bounds_(path)
            lmat[idx1, idx2] = lmat[idx2, idx1] = ldist
            umat[idx1, idx2] = umat[idx2, idx1] = udist
        else:
            # they are disconnected
            lmat[idx1, idx2] = lmat[idx2, idx1] = closest_approach(
                gra, key1, key2)
            umat[idx1, idx2] = umat[idx2, idx1] = 999

        assert lmat[idx1, idx2] <= umat[idx1, idx2], (
            "Lower bound exceeds upper bound. This is a bug!\n"
            f"{string(gra, one_indexed=False)}\npath: {str(path)}\n")

    return lmat, umat


def join_distance_bounds_matrices(gra, keys, dist_range_dct, geos=None,
                                  relax_angles=False, relax_torsions=False,
                                  sp_dct=None, angstrom=True):
    """ distance bounds matrices for joining multiple geometries

    :param gra: molecular graph:
    :param keys: atom keys specifying the order of indices in the matrix
    :param dist_range_dct: distance ranges for specific atoms in the graph
    :param geos: (optional) geometries which will be used to fix the bond
        angles, bond distances, and chiralities of all connected atoms in the
        graph; if `relax_torsions` is False, the 4-atom dihedral angles will be
        allowed to vary as well
    :param relax_torsions: whether or not to allow torsions to change from
        their value in the reactant geometries
    :param sp_dct: a 2d dictionary giving the shortest path between any pair of
        atoms in the graph
    """
    sp_dct = atom_shortest_paths(gra) if sp_dct is None else sp_dct

    natms = len(keys)

    lmat, umat = distance_bounds_matrices(gra, keys)

    # save the current values so that we can overwrite the fixed torsions below
    lmat_old = numpy.copy(lmat)
    umat_old = numpy.copy(umat)

    # 1. set known geometric parameters
    if geos:
        xmats = [automol.geom.base.coordinates(geo, angstrom=angstrom)
                 for geo in geos]
        dmats = list(map(embed.distance_matrix_from_coordinates, xmats))

        start = 0
        for dmat in dmats:
            dim, _ = numpy.shape(dmat)
            end = start + dim

            lmat[start:end, start:end] = dmat
            umat[start:end, start:end] = dmat

            start = end

    # 2. reopen bounds on the torsions from the reactant
    if relax_torsions:
        tors_ijs = [[i, j] for i, j in itertools.combinations(range(natms), 2)
                    if j in sp_dct[i] and len(sp_dct[i][j]) >= 4]
        tors_ijs += list(map(list, map(reversed, tors_ijs)))

        tors_idxs = tuple(map(list, zip(*tors_ijs)))

        lmat[tors_idxs] = lmat_old[tors_idxs]
        umat[tors_idxs] = umat_old[tors_idxs]

    # 3. reopen bounds on the angles from the reactant
    if relax_angles:
        ang_ijs = [[i, j] for i, j in itertools.combinations(range(natms), 2)
                   if j in sp_dct[i] and len(sp_dct[i][j]) >= 3]
        ang_ijs += list(map(list, map(reversed, ang_ijs)))

        ang_idxs = tuple(map(list, zip(*ang_ijs)))

        lmat[ang_idxs] = lmat_old[ang_idxs]
        umat[ang_idxs] = umat_old[ang_idxs]

    # 4. set distance bounds for the forming bonds
    for bnd, (ldist, udist) in dist_range_dct.items():
        idx1 = tuple(bnd)
        idx2 = tuple(reversed(idx1))
        lmat[idx1] = lmat[idx2] = ldist
        umat[idx1] = umat[idx2] = udist

    return lmat, umat


def ts_distance_bounds_matrices(gra, keys, frm_bnds_dct, rct_geos=None,
                                relax_torsions=False, sp_dct=None):
    """ distance bounds matrices for a transition state

    :param gra: molecular graph
    :param keys: atom keys specifying the order of indices in the matrix
    :param frm_bnds_dct: bounds for bonds formed in the reaction; the keys are
        bond keys and the values are lower and upper distance bounds for those
        forming bonds
    :param rct_geos: reactant geometries; must follow the same order as the
        atoms in `keys
    :param relax_torsions: whether or not to allow torsions to change from
        their value in the reactant geometries
    :param sp_dct: a 2d dictionary giving the shortest path between any pair of
        atoms in the graph
    """
    sp_dct = atom_shortest_paths(gra) if sp_dct is None else sp_dct

    natms = len(keys)

    lmat, umat = distance_bounds_matrices(gra, keys)

    # save the current values so that we can overwrite the fixed torsions below
    lmat_old = numpy.copy(lmat)
    umat_old = numpy.copy(umat)

    # 2. set known geometric parameters
    if rct_geos:
        xmats = [automol.geom.base.coordinates(geo, angstrom=True)
                 for geo in rct_geos]
        dmats = list(map(embed.distance_matrix_from_coordinates, xmats))

        start = 0
        for dmat in dmats:
            dim, _ = numpy.shape(dmat)
            end = start + dim

            lmat[start:end, start:end] = dmat
            umat[start:end, start:end] = dmat

            start = end

    # 3. reopen bounds on the torsions from the reactant
    if relax_torsions:
        tors_ijs = [[i, j] for i, j in itertools.combinations(range(natms), 2)
                    if j in sp_dct[i] and len(sp_dct[i][j]) >= 4]
        tors_ijs += list(map(list, map(reversed, tors_ijs)))

        tors_idxs = tuple(map(list, zip(*tors_ijs)))

        lmat[tors_idxs] = lmat_old[tors_idxs]
        umat[tors_idxs] = umat_old[tors_idxs]

    # 1. set distance bounds for the forming bonds
    for bnd, (ldist, udist) in frm_bnds_dct.items():
        idx1 = tuple(bnd)
        idx2 = tuple(reversed(idx1))
        lmat[idx1] = lmat[idx2] = ldist
        umat[idx1] = umat[idx2] = udist

    return lmat, umat


# # constraint dictionaries
def chirality_constraint_bounds(gra, keys):
    """ bounds for enforcing chirality restrictions
    """
    ste_keys = set(atom_stereo_keys(gra)) & set(keys)
    par_dct = atom_stereo_parities(gra)
    ngb_key_dct = atoms_neighbor_atom_keys(gra)

    def _chirality_constraint(key):
        ngb_keys = ngb_key_dct[key]
        ngb_keys = atom_stereo_sorted_neighbor_atom_keys(gra, key, ngb_keys)
        idxs = tuple(map(keys.index, ngb_keys))
        vol_range = (-999., -7.) if par_dct[key] else (+7., +999.)
        return idxs, vol_range

    chi_dct = dict(map(_chirality_constraint, ste_keys))
    return chi_dct


def planarity_constraint_bounds(gra, keys):
    """ bounds for enforcing planarity restrictions
    """
    ngb_key_dct = atoms_neighbor_atom_keys(gra)
    ngb_dct = bond_neighborhoods(gra)
    bnd_keys = [bnd_key for bnd_key in sp2_bond_keys(gra)
                if atom_keys(ngb_dct[bnd_key]) <= set(keys)]

    def _planarity_constraints(bnd_key):
        key1, key2 = sorted(bnd_key)
        key1ab = sorted(ngb_key_dct[key1] - {key2})
        key2ab = sorted(ngb_key_dct[key2] - {key1})

        lst = []

        # I don't think the order of the keys matters, but I tried to be
        # roughly consistent with Figure 8 in the Blaney Dixon paper
        if len(key1ab) == 2 and len(key2ab) == 2:
            lst.append(tuple(map(keys.index, key1ab + key2ab)))
        if len(key1ab) == 2:
            lst.append(tuple(map(keys.index, [key1, key2] + key1ab)))
        if len(key2ab) == 2:
            lst.append(tuple(map(keys.index, [key1, key2] + key2ab)))
        if (len(key1ab) == 2 and len(key2ab) == 1) or (
                len(key1ab) == 1 and len(key2ab) == 2):
            lst.append(tuple(map(keys.index, [key1] + key1ab + key2ab)))
            lst.append(tuple(map(keys.index, [key2] + key1ab + key2ab)))
        if len(key1ab) == 1 and len(key2ab) == 1:
            lst.append(tuple(map(keys.index, [key1, key2] + key1ab + key2ab)))

        return tuple(lst)

    const_dct = {
        idxs: (-0.5, +0.5) for idxs in
        itertools.chain(*map(_planarity_constraints, bnd_keys))}

    return const_dct


# # distance range dictionaries
def distance_ranges_from_coordinates(gra, dist_dct, ang_dct=None, dih_dct=None,
                                     angstrom=True, degree=True,
                                     rings_keys=(), keys=None, check=False):
    """ generate a set of distance ranges from coordinate values

    :param gra: molecular graph
    atom keys specifying the order of indices in the matrix
    :param dist_dct: a dictionary of desired distances for certain atoms; the
        keys are pairs of atoms, the values are distances in angstroms
    :type dist_dct: dict[(int, int): float]
    :param ang_dct: a dictionary of desired angles for certain atoms; the keys
        are triples of atoms; if the first or last element in a triple is None,
        an appopriate neighboring atom will be found
    :param dih_dct: a dictionary of desired angles for certain atoms; the keys
        are quadruples of atoms; if the first or last element in a triple is
        None, an appopriate neighboring atom will be found
    :type dist_dct: dict[(int, int, int): float]
    :param rings_keys: keys for rings in the graph; angle ranges will
        automatically be set to allow ring formation
    :param keys: set of keys that can be used to fill in the angle keys; if
        None, all graph keys will be considered available for use
    :param check: check the angle keys to make sure they can all be filled in?
    """
    keys = atom_keys(gra) if keys is None else keys

    ang_dct = {} if ang_dct is None else ang_dct
    dih_dct = {} if dih_dct is None else dih_dct

    # Fill in angle keys
    ang_key_filler_ = angle_key_filler_(gra, keys, check=check)
    ang_dct = dict_.transform_keys(ang_dct, ang_key_filler_)
    if None in ang_dct:
        ang_dct.pop(None)

    # Fill in dihedral keys
    dih_dct = dict_.transform_keys(dih_dct, ang_key_filler_)
    if None in dih_dct:
        dih_dct.pop(None)

    # Convert angles into distances
    dist_dct = dict_.transform_keys(dist_dct, frozenset)
    for (key1, key2, key3), a123 in ang_dct.items():
        a123 *= phycon.DEG2RAD if degree else 1.

        k12 = frozenset({key1, key2})
        k23 = frozenset({key2, key3})
        k13 = frozenset({key1, key3})

        d12 = (dist_dct[k12] if k12 in dist_dct else
               heuristic_bond_distance(gra, key1, key2, angstrom=angstrom))
        d23 = (dist_dct[k23] if k23 in dist_dct else
               heuristic_bond_distance(gra, key2, key3, angstrom=angstrom))

        d13 = numpy.sqrt(d12**2 + d23**2 - 2*d12*d23*numpy.cos(a123))

        dist_dct[k13] = d13

    # Convert convert fixed distances into ranges
    dist_range_dct = {k: (d, d) for k, d in dist_dct.items()}

    # Convert dihedrals into distances
    for (key1, key2, key3, key4), val in dih_dct.items():
        # Allow user to leave dihedrals open-ended, as a lower or upper bound
        if isinstance(val, numbers.Number):
            d1234 = val
        else:
            assert hasattr(val, '__len__') and len(val) == 2
            d1234 = next(v for v in val if v is not None)

        d1234 *= phycon.DEG2RAD if degree else 1.

        k12 = frozenset({key1, key2})
        k23 = frozenset({key2, key3})
        k34 = frozenset({key3, key4})
        k13 = frozenset({key1, key3})
        k24 = frozenset({key2, key4})
        k14 = frozenset({key1, key4})

        d12 = (dist_dct[k12] if k12 in dist_dct else
               heuristic_bond_distance(gra, key1, key2, angstrom=angstrom))
        d23 = (dist_dct[k23] if k23 in dist_dct else
               heuristic_bond_distance(gra, key2, key3, angstrom=angstrom))
        d34 = (dist_dct[k34] if k34 in dist_dct else
               heuristic_bond_distance(gra, key3, key4, angstrom=angstrom))
        d13 = (dist_dct[k13] if k13 in dist_dct else
               heuristic_bond_distance(gra, key1, key3, angstrom=angstrom))
        d24 = (dist_dct[k24] if k24 in dist_dct else
               heuristic_bond_distance(gra, key2, key4, angstrom=angstrom))

        term1 = (d12**2 + d23**2 - d13**2)*(d23**2 + d34**2 - d24**2)
        term2 = 2*d23**2 * (d13**2 + d24**2 - d23**2)
        denom = numpy.sqrt(
            (4*d12**2 * d23**2 - (d12**2 + d23**2 - d13**2)**2) *
            (4*d23**2 * d34**2 - (d23**2 + d34**2 - d24**2)**2))

        d14 = numpy.sqrt((term1 + term2 - numpy.cos(d1234) * denom) /
                         (2 * d23**2))

        if isinstance(val, numbers.Number) or val[0] == val[1]:
            dist_range_dct[k14] = (d14, d14)
        elif val[0] is None:
            ld14 = closest_approach(gra, key1, key4)
            dist_range_dct[k14] = (ld14, d14)
        elif val[1] is None:
            ud14 = 999.
            dist_range_dct[k14] = (d14, ud14)
        else:
            raise ValueError(f"Invalid dih_dict: {str(dih_dct)}")

    for rng_keys in rings_keys:
        assert hasattr(keys, '__iter__'), (
            "Please pass in rings keys as a list of lists")

        rsz = len(rng_keys)
        a123 = (rsz - 2.) * 180. / rsz
        la123 = (a123 - 10.) * phycon.DEG2RAD
        ua123 = (a123 + 10.) * phycon.DEG2RAD

        for key1, key2, key3 in mit.windowed(rng_keys + rng_keys[:2], 3):
            k12 = frozenset({key1, key2})
            k23 = frozenset({key2, key3})
            k13 = frozenset({key1, key3})

            d12 = (dist_dct[k12] if k12 in dist_dct else
                   heuristic_bond_distance(gra, key1, key2, angstrom=angstrom))
            d23 = (dist_dct[k23] if k23 in dist_dct else
                   heuristic_bond_distance(gra, key2, key3, angstrom=angstrom))

            ld13 = numpy.sqrt(d12**2 + d23**2 - 2*d12*d23*numpy.cos(la123))
            ud13 = numpy.sqrt(d12**2 + d23**2 - 2*d12*d23*numpy.cos(ua123))
            dist_range_dct[k13] = (ld13, ud13)

    return dist_range_dct


# # heuristic coordinate values
def heuristic_bond_distance(gra, key1, key2, angstrom=True, check=False):
    """ heuristic bond distance (in angstroms)
    """
    if check:
        assert key1 in atoms_neighbor_atom_keys(gra)[key2]

    symb_dct = atom_symbols(gra)
    symb1 = symb_dct[key1]
    symb2 = symb_dct[key2]

    if ptab.to_number(symb1) == 1 or ptab.to_number(symb2) == 1:
        dist = XH_DIST
    else:
        dist = XY_DIST

    dist *= 1 if angstrom else phycon.ANG2BOHR

    return dist


def heuristic_bond_angle(gra, key1, key2, key3, degree=False, check=False,
                         hyb_dct=None):
    """ heuristic bond angle

    If being reused multiple times, you can speed this up by passing in the
    hybridizations, so they don't need to be recalculated
    """
    if check:
        assert {key1, key3} <= set(atoms_neighbor_atom_keys(gra)[key2])

    if hyb_dct is None:
        hyb_dct = resonance_dominant_atom_hybridizations(gra)

    hyb2 = hyb_dct[key2]
    if hyb2 == 3:
        ang = TET_ANG
    elif hyb2 == 2:
        ang = TRI_ANG
    else:
        assert hyb2 == 1
        ang = LIN_ANG

    ang *= 1 if degree else phycon.DEG2RAD

    return ang


def heuristic_bond_angle_distance(gra, key1, key2, key3,
                                  d12=None, d23=None, angstrom=True,
                                  a123=None, degree=True, hyb_dct=None):
    """ heuristic distance between atoms at two ends of a bond angle

    :param angstrom: whether or not to return the distance in angstroms
    :type angstrom: bool
    :param a123: (optional) specify the value of the angle
    :type a123: float
    :param degree: units for the angle, if specified
    :type degree: bool

    uses the law of cosines:

        d13 = sqrt(d12^2 + d23^2 - 2*d12*d23*cos(a123))
    """
    if a123 is None:
        a123 = heuristic_bond_angle(
            gra, key1, key2, key3, degree=False, hyb_dct=hyb_dct)
    else:
        a123 = a123 * phycon.DEG2RAD if degree else a123

    d12 = heuristic_bond_distance(gra, key1, key2, angstrom=angstrom)
    d23 = heuristic_bond_distance(gra, key2, key3, angstrom=angstrom)
    d13 = numpy.sqrt(d12**2 + d23**2 - 2*d12*d23*numpy.cos(a123))
    return d13


def heuristic_torsion_angle_distance(gra, key1, key2, key3, key4,
                                     angstrom=True,
                                     a123=None, a234=None, d1234=None,
                                     degree=True, hyb_dct=None):
    """ heuristic max distance between atoms at two ends of a torsion angle

    Formula rearranged from eq. 2.15 in the following paper. Note that in the
    denominator of this formula there is a typo: d02 and d12 should be switched

    Havel, T. F.; Kunz, I.D.; Crippen, G. M.; "The Theory and Practice of
    Distance Geometry"; Bulletin of Mathematical Biology; Vol. 45, No. 5.,
    pp.  665-720, 1983.
    """
    if d1234 is None:
        d1234 = numpy.pi
    else:
        d1234 *= phycon.DEG2RAD if degree else 1.

    d12 = heuristic_bond_distance(gra, key1, key2, angstrom=angstrom)
    d23 = heuristic_bond_distance(gra, key2, key3, angstrom=angstrom)
    d34 = heuristic_bond_distance(gra, key3, key4, angstrom=angstrom)
    d13 = heuristic_bond_angle_distance(
        gra, key1, key2, key3, angstrom=angstrom, a123=a123, hyb_dct=hyb_dct)
    d24 = heuristic_bond_angle_distance(
        gra, key2, key3, key4, angstrom=angstrom, a123=a234, hyb_dct=hyb_dct)

    term1 = (d12**2 + d23**2 - d13**2)*(d23**2 + d34**2 - d24**2)
    term2 = 2*d23**2 * (d13**2 + d24**2 - d23**2)
    denom = numpy.sqrt((4*d12**2 * d23**2 - (d12**2 + d23**2 - d13**2)**2) *
                       (4*d23**2 * d34**2 - (d23**2 + d34**2 - d24**2)**2))

    d14 = numpy.sqrt((term1 + term2 - numpy.cos(d1234) * denom) /
                     (2 * d23**2))
    return d14


def closest_approach(gra, key1, key2):
    """ closest approach between atoms, based on their van der Waals radii

    Warning: The scaling factor on the van der waals radii was arbitrarily
    chosen based on limited tests and may need to be lowered
    """
    vdw_scaling_factor = 0.75
    dist = (atom_van_der_waals_radius(gra, key1) +
            atom_van_der_waals_radius(gra, key2)) * vdw_scaling_factor
    return dist


def path_distance_bounds_(gra):
    """ upper distance bound between two ends of a path

    :param gra: molecular graph
    :param path: the shortest path between two atoms
    :type path: list or tuple
    """

    hyb_dct = resonance_dominant_atom_hybridizations(gra)
    rng_keys_lst = rings_atom_keys(gra)
    atm_ngb_keys = atoms_neighbor_atom_keys(gra)

    ste_bnd_keys = bond_stereo_keys(gra)
    bnd_par_dct = bond_stereo_parities(gra)

    def _distance_bounds(path):
        # if the path is 0, the atoms are disconnected and could be arbitrarily
        # far apart
        rsz = shared_ring_size(path, rng_keys_lst)
        if len(path) == 1:
            ldist = udist = 0
        elif len(path) == 2:
            ldist = udist = heuristic_bond_distance(gra, *path)
        elif len(path) == 3:
            if rsz == 0:
                ldist = udist = heuristic_bond_angle_distance(
                    gra, *path, hyb_dct=hyb_dct)
            else:
                a123 = (rsz - 2.) * 180. / rsz
                la123 = a123 - 10.
                ua123 = a123 + 10.
                lrdist = heuristic_bond_angle_distance(gra, *path, a123=la123)
                urdist = heuristic_bond_angle_distance(gra, *path, a123=ua123)
                odist = heuristic_bond_angle_distance(
                    gra, *path, hyb_dct=hyb_dct)
                ldist = min(lrdist, odist)
                udist = max(urdist, odist)
        elif len(path) == 4:
            if rsz == 0:
                key2, key3 = path[1:3]
                bnd_key23 = frozenset({key2, key3})
                # handle bond stereo here
                if bnd_key23 in ste_bnd_keys:
                    key2_ngbs = atom_stereo_sorted_neighbor_atom_keys(
                        gra, key2, atm_ngb_keys[key2]-{key3})
                    key3_ngbs = atom_stereo_sorted_neighbor_atom_keys(
                        gra, key3, atm_ngb_keys[key3]-{key2})
                    pos2 = key2_ngbs.index(path[0])
                    pos3 = key3_ngbs.index(path[-1])
                    cis = bnd_par_dct[bnd_key23] != (pos2 != pos3)
                    dih = 0. if cis else 180.
                    ldist = udist = heuristic_torsion_angle_distance(
                        gra, *path, d1234=dih, degree=True, hyb_dct=hyb_dct)
                else:
                    ldist = heuristic_torsion_angle_distance(
                        gra, *path, d1234=0., degree=True, hyb_dct=hyb_dct)
                    udist = heuristic_torsion_angle_distance(
                        gra, *path, d1234=180., degree=True, hyb_dct=hyb_dct)
            else:
                ang = (rsz - 2.) * 180. / rsz
                rdist = heuristic_torsion_angle_distance(
                    gra, *path, d1234=0., degree=True, a123=ang, a234=ang)
                cdist = heuristic_torsion_angle_distance(
                    gra, *path, d1234=0., degree=True, hyb_dct=hyb_dct)
                tdist = heuristic_torsion_angle_distance(
                    gra, *path, d1234=180., degree=True, hyb_dct=hyb_dct)
                ldist = min(rdist, cdist)
                udist = max(rdist, tdist)
        # otherwise, just do the sum of the distances between atoms along the
        # path
        else:
            # we can't handle disconnected points, because in that case the
            # path is [] and there is no way to recover the keys
            assert len(path) > 2
            ldist = closest_approach(gra, path[0], path[-1])
            udist = 999.

        return ldist, udist

    return _distance_bounds


# # helpers
def angle_key_filler_(gra, keys=None, check=True):
    """ returns a function that fills in the first or last element of an angle
    key in a dictionary with a neighboring atom
    (works for central or dihedral angles)
    """
    keys = atom_keys(gra) if keys is None else keys

    def _fill_in_angle_key(ang_key):
        end1_key = ang_key[0]
        end2_key = ang_key[-1]
        mid_keys = list(ang_key[1:-1])
        assert not any(k is None for k in mid_keys)

        if end1_key is None:
            end1_key = atom_neighbor_atom_key(
                gra, mid_keys[0], excl_atm_keys=[end2_key]+mid_keys,
                incl_atm_keys=keys)
        if end2_key is None:
            end2_key = atom_neighbor_atom_key(
                gra, mid_keys[-1], excl_atm_keys=[end1_key]+mid_keys,
                incl_atm_keys=keys)

        ang_key = [end1_key] + mid_keys + [end2_key]

        if any(k is None for k in ang_key):
            if check:
                raise ValueError(
                    f"Angle key {str(ang_key)} couldn't be filled in")
            ang_key = None
        else:
            ang_key = tuple(ang_key)

        return ang_key

    return _fill_in_angle_key


def shared_ring_size(keys, rng_keys_lst):
    """ determine whether these keys share a ring and, if so, determine the
    size
    """
    rng_keys = next((rng_keys for rng_keys in rng_keys_lst
                     if set(keys) <= set(rng_keys)), ())
    natms = len(rng_keys)
    return natms
