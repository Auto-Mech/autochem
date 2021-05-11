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
from automol import embed
from automol.util import dict_
from automol.graph.geom import coordinates as geom_coordinates
from automol.graph._graph_dep import atom_keys
from automol.graph._graph_dep import atom_neighbor_atom_key
from automol.graph._embed_dep import atom_shortest_paths
from automol.graph._embed_dep import closest_approach
from automol.graph._embed_dep import (
    fake_stereo_geometry as _fake_stereo_geometry)
from automol.graph._embed_dep import geometry as _geometry
from automol.graph._embed_dep import(
    distance_bounds_matrices as _distance_bounds_matrices)
from automol.graph._embed_dep import (
    chirality_constraint_bounds as _chirality_constraint_bounds)
from automol.graph._embed_dep import (
    planarity_constraint_bounds as _planarity_constraint_bounds)
from automol.graph._embed_dep import heuristic_bond_distance


# bond distances
XY_DIST = 1.5       # angstroms
XH_DIST = 1.1       # angstroms

# bond angles
TET_ANG = 109.4712  # degrees
TRI_ANG = 120.      # degrees
LIN_ANG = 180.      # degrees


def fake_stereo_geometry(gra, ntries=5, max_dist_err=0.5):
    return _fake_stereo_geometry(
        gra, ntries=ntries, max_dist_err=max_dist_err)


def geometry(gra, keys=None, ntries=5, max_dist_err=0.2):
    return _geometry(
        gra, keys=keys, ntries=ntries, max_dist_err=max_dist_err)


def distance_bounds_matrices(gra, keys, sp_dct=None):
    return _distance_bounds_matrices(gra, keys, sp_dct=sp_dct)


def chirality_constraint_bounds(gra, keys):
    return _chirality_constraint_bounds(gra, keys)


def planarity_constraint_bounds(gra, keys):
    return _planarity_constraint_bounds(gra, keys)


def join_distance_bounds_matrices(gra, keys, dist_range_dct, geos=None,
                                  relax_angles=False, relax_torsions=False,
                                  sp_dct=None, angstrom=True):
    """ distance bounds matrices for joining multiple geometries

    :param gra: molecular graph:wq
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
        xmats = [geom_coordinates(geo, angstrom=angstrom)
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
        xmats = [geom_coordinates(geo, angstrom=True)
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
            raise ValueError("Invalid dih_dict: {}".format(str(dih_dct)))

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
                raise ValueError("Angle key {} couldn't be filled in"
                                 .format(str(ang_key)))
            ang_key = None
        else:
            ang_key = tuple(ang_key)

        return ang_key

    return _fill_in_angle_key


#
# if __name__ == '__main__':
#     import automol
#     ICH = 'InChI=1S/C5H10O3/c1-4-2-5(8-4)3-7-6/h4-6H,2-3H2,1H3/t4-,5-/m1/s1'
#     GEO = automol.inchi.geometry(ICH)
#     GRA = automol.geom.graph(GEO)
#     GROUP = [0, 1, 2, 3, 4, 7, 15, 16]
#     SUBGEO = automol.geom.from_subset(GEO, GROUP)
#     print(automol.geom.string(SUBGEO))
