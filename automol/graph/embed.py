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
import itertools
import numpy
import qcelemental as qce
from qcelemental import constants as qcc
from qcelemental import periodictable as pt
import automol.create.geom
import automol.geom
from automol import embed
from automol.graph._graph_base import string
from automol.graph._graph import atom_symbols
from automol.graph._graph import atom_keys
from automol.graph._graph import atom_shortest_paths
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import bond_stereo_parities
from automol.graph._ring import rings_atom_keys
from automol.graph._res import resonance_dominant_atom_hybridizations
from automol.graph._stereo import bond_stereo_keys
from automol.graph._stereo import sp2_bond_keys
from automol.graph._stereo import stereo_sorted_atom_neighbor_keys


ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
DEG2RAD = qcc.conversion_factor('degree', 'radian')
RAD2DEG = qcc.conversion_factor('radian', 'degree')

# bond distances
XY_DIST = 1.5       # angstroms
XH_DIST = 1.1       # angstroms

# bond angles
TET_ANG = 109.4712  # degrees
TRI_ANG = 120.      # degrees
LIN_ANG = 180.      # degrees


def heuristic_bond_distance(gra, key1, key2, angstrom=True, check=False):
    """ heuristic bond distance (in angstroms)
    """
    if check:
        assert key1 in atom_neighbor_keys(gra)[key2]

    sym_dct = atom_symbols(gra)
    sym1 = sym_dct[key1]
    sym2 = sym_dct[key2]

    if pt.to_Z(sym1) == 1 or pt.to_Z(sym2) == 1:
        dist = XH_DIST
    else:
        dist = XY_DIST

    dist *= 1 if angstrom else ANG2BOHR

    return dist


def heuristic_bond_angle(gra, key1, key2, key3, degree=False, check=False,
                         hyb_dct=None):
    """ heuristic bond angle

    If being reused multiple times, you can speed this up by passing in the
    hybridizations, so they don't need to be recalculated
    """
    if check:
        assert {key1, key3} <= set(atom_neighbor_keys(gra)[key2])

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

    ang *= 1 if degree else DEG2RAD

    return ang


def heuristic_bond_angle_distance(gra, key1, key2, key3, angstrom=True,
                                  ang=None, degree=True, hyb_dct=None):
    """ heuristic distance between atoms at two ends of a bond angle

    :param angstrom: whether or not to return the distance in angstroms
    :type angstrom: bool
    :param ang: (optional) specify the value of the angle
    :type ang: float
    :param degree: units for the angle, if specified
    :type degree: bool

    uses the law of cosines:

        d13 = sqrt(d12^2 + d23^2 - 2*d12*d23*cos(a123))
    """
    if ang is None:
        a123 = heuristic_bond_angle(
            gra, key1, key2, key3, degree=False, hyb_dct=hyb_dct)
    else:
        a123 = ang * DEG2RAD if degree else ang

    d12 = heuristic_bond_distance(gra, key1, key2, angstrom=angstrom)
    d23 = heuristic_bond_distance(gra, key2, key3, angstrom=angstrom)
    d13 = numpy.sqrt(d12**2 + d23**2 - 2*d12*d23*numpy.cos(a123))
    return d13


def heuristic_torsion_angle_distance(gra, key1, key2, key3, key4, cis=False,
                                     angstrom=True,
                                     ang1=None, ang2=None, degree=True,
                                     hyb_dct=None):
    """ heuristic max distance between atoms at two ends of a torsion angle

    (Formulas & implementation have been verified)

    max implies that the torsion angle is 180 degrees

    uses law of cosines:

        d14 = sqrt(d13^2 + d34^2 - 2*d13*d34*cos(a134))

    in this case, d13 is determined from heuristic_bond_angle_distance and a134
    is given by

        a134 = a234 + a132

    if the atoms are trans and

        a134 = a234 - a132

    if they are cis

    we can find a132 from the law of cosines:
        a132 = arccos((d13^2 + d23^2 - d12^2)/(2*d13*d23))
    """
    if ang2 is None:
        a234 = heuristic_bond_angle(
            gra, key2, key3, key4, degree=False, hyb_dct=hyb_dct)
    else:
        a234 = ang2 * DEG2RAD if degree else ang2

    d12 = heuristic_bond_distance(gra, key1, key2, angstrom=angstrom)
    d23 = heuristic_bond_distance(gra, key2, key3, angstrom=angstrom)
    d34 = heuristic_bond_distance(gra, key3, key4, angstrom=angstrom)
    d13 = heuristic_bond_angle_distance(
        gra, key1, key2, key3, angstrom=angstrom, ang=ang1, hyb_dct=hyb_dct)

    a132 = numpy.arccos((d13**2 + d23**2 - d12**2)/(2*d13*d23))
    a134 = a234 - a132 if cis else a234 + a132

    d14 = numpy.sqrt(d13**2 + d34**2 - 2*d13*d34*numpy.cos(a134))
    return d14


def van_der_waals_radius(gra, key):
    """ van der waals radius for an atom in the graph (in angstroms)
    """
    sym_dct = atom_symbols(gra)
    sym = sym_dct[key]
    rad = qce.vdwradii.get(sym, units='angstrom')
    return rad


def shared_ring_size(keys, rng_keys_lst):
    """ determine whether these keys share a ring and, if so, determine the
    size
    """
    rng_keys = next((rng_keys for rng_keys in rng_keys_lst
                     if set(keys) <= set(rng_keys)), ())
    natms = len(rng_keys)
    return natms


def closest_approach(gra, key1, key2):
    """ closest approach between atoms, based on their van der Waals radii

    Warning: The scaling factor on the van der waals radii was arbitrarily
    chosen based on limited tests and may need to be lowered
    """
    vdw_scaling_factor = 0.8
    dist = (van_der_waals_radius(gra, key1) +
            van_der_waals_radius(gra, key2)) * vdw_scaling_factor
    return dist


def path_distance_bounds_(gra):
    """ upper distance bound between two ends of a path

    :param gra: molecular graph
    :param path: the shortest path between two atoms
    :type path: list or tuple
    """

    hyb_dct = resonance_dominant_atom_hybridizations(gra)
    rng_keys_lst = rings_atom_keys(gra)
    atm_ngb_keys = atom_neighbor_keys(gra)

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
                ang = (rsz - 2.) * 180. / rsz
                rdist = heuristic_bond_angle_distance(gra, *path, ang=ang)
                odist = heuristic_bond_angle_distance(
                    gra, *path, hyb_dct=hyb_dct)
                ldist = min(rdist, odist)
                udist = max(rdist, odist)
        elif len(path) == 4:
            if rsz == 0:
                key2, key3 = path[1:3]
                bnd_key23 = frozenset({key2, key3})
                # handle bond stereo here
                if bnd_key23 in ste_bnd_keys:
                    key2_ngbs = stereo_sorted_atom_neighbor_keys(
                        gra, key2, atm_ngb_keys[key2]-{key3})
                    key3_ngbs = stereo_sorted_atom_neighbor_keys(
                        gra, key3, atm_ngb_keys[key3]-{key2})
                    pos2 = key2_ngbs.index(path[0])
                    pos3 = key3_ngbs.index(path[-1])
                    cis = bnd_par_dct[bnd_key23] != (pos2 != pos3)
                    ldist = udist = heuristic_torsion_angle_distance(
                        gra, *path, cis=cis, hyb_dct=hyb_dct)
                else:
                    ldist = heuristic_torsion_angle_distance(
                        gra, *path, cis=True, hyb_dct=hyb_dct)
                    udist = heuristic_torsion_angle_distance(
                        gra, *path, cis=False, hyb_dct=hyb_dct)
            else:
                ang = (rsz - 2.) * 180. / rsz
                rdist = heuristic_torsion_angle_distance(
                    gra, *path, cis=True, ang1=ang, ang2=ang)
                cdist = heuristic_torsion_angle_distance(
                    gra, *path, cis=True, hyb_dct=hyb_dct)
                tdist = heuristic_torsion_angle_distance(
                    gra, *path, cis=False, hyb_dct=hyb_dct)
                ldist = min(rdist, cdist)
                udist = max(rdist, tdist)
        # otherwise, just do the sum of the distances between atoms along the
        # path
        else:
            # we can't handle disconnected points, because in that case the
            # path is [] and there is no way to recover the keys
            assert len(path) > 2
            ldist = closest_approach(gra, path[0], path[-1])
            udist = 999

        return ldist, udist

    return _distance_bounds


def distance_bounds_matrices(gra, keys, sp_dct=None):
    """ initial distance bounds matrices

    :param gra: molecular graph
    :param keys: atom keys specifying the order of indices in the matrix
    :param sp_dct: a 2d dictionary giving the shortest path between any pair of
        atoms in the graph
    """
    assert set(keys) == set(atom_keys(gra))

    sp_dct = atom_shortest_paths(gra) if sp_dct is None else sp_dct

    bounds_ = path_distance_bounds_(gra)

    natms = len(keys)
    umat = numpy.zeros((natms, natms))
    lmat = numpy.zeros((natms, natms))
    for i, j in itertools.combinations(range(natms), 2):
        if j in sp_dct[i]:
            path = sp_dct[i][j]
            ldist, udist = bounds_(path)
            lmat[i, j] = lmat[j, i] = ldist
            umat[i, j] = umat[j, i] = udist
        else:
            # they are disconnected
            lmat[i, j] = lmat[j, i] = closest_approach(gra, keys[i], keys[j])
            umat[i, j] = umat[j, i] = 999

        assert lmat[i, j] <= umat[i, j], (
            "Lower bound exceeds upper bound. This is a bug!\n"
            "{}\npath: {}\n"
            .format(string(gra, one_indexed=False), str(path)))

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
        xmats = [automol.geom.coordinates(geo, angstrom=True)
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


# def sample_raw_distance_coordinates(gra, keys, dim4=False):
#     """ sample raw (uncorrected) distance coordinates
#     """
#     # 1. Generate distance bounds matrices, L and U
#     lmat, umat = distance_bounds_matrices(gra, keys)
#
#     # 2-6. Triangle-smooth the bounds matrices
#     xmat = embed.sample_raw_distance_coordinates(lmat, umat, dim4=dim4)
#
#     return xmat


def chirality_constraint_bounds(gra, keys):
    """ bounds for enforcing chirality restrictions
    """
    print(gra)
    print(keys)


def planarity_constraint_bounds(gra, keys):
    """ bounds for enforcing planarity restrictions
    """
    ngb_key_dct = atom_neighbor_keys(gra)

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
        idxs: (-0.2, 0.2) for idxs in
        itertools.chain(*map(_planarity_constraints, sp2_bond_keys(gra)))}

    return const_dct


def sample_raw_distance_geometry(gra, keys):
    """ sample a raw (uncorrected) distance geometry
    """
    # 1. Generate distance bounds matrices, L and U
    lmat, umat = distance_bounds_matrices(gra, keys)

    # 2-6. Generate coordinates from bounds matrices
    xmat = embed.sample_raw_distance_coordinates(lmat, umat, dim4=True)

    sym_dct = atom_symbols(gra)
    syms = list(map(sym_dct.__getitem__, keys))
    geo = automol.create.geom.from_data(syms, xmat, angstrom=True)

    return geo


if __name__ == '__main__':
    import automol
    # ICH = ('InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)'
    #        '4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/'
    #        't10-,11+,13-,16-,17-/m0/s1')
    # GEO = automol.inchi.geometry(ICH)
    # GRA = automol.geom.graph(GEO)
    GRA = automol.graph.from_string("""
        atoms:
          1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          3: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          4: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          5: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          6: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          7: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          8: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          9: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          10: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: true}
          11: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: true}
          12: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          13: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: true}
          14: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          15: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          16: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: false}
          17: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: false}
          18: {symbol: N, implicit_hydrogen_valence: 0, stereo_parity: null}
          19: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
          20: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
          21: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
          22: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          23: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          24: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          25: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          26: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          27: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          28: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          29: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          30: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          31: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          32: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          33: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          34: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          35: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          36: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          37: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          38: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          39: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
          40: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
        bonds:
          1-18: {order: 1, stereo_parity: null}
          1-22: {order: 1, stereo_parity: null}
          1-23: {order: 1, stereo_parity: null}
          1-24: {order: 1, stereo_parity: null}
          2-4: {order: 1, stereo_parity: null}
          2-9: {order: 1, stereo_parity: null}
          2-25: {order: 1, stereo_parity: null}
          3-5: {order: 1, stereo_parity: null}
          3-10: {order: 1, stereo_parity: null}
          3-26: {order: 1, stereo_parity: null}
          4-12: {order: 1, stereo_parity: null}
          4-27: {order: 1, stereo_parity: null}
          5-13: {order: 1, stereo_parity: null}
          5-28: {order: 1, stereo_parity: null}
          6-7: {order: 1, stereo_parity: null}
          6-17: {order: 1, stereo_parity: null}
          6-29: {order: 1, stereo_parity: null}
          6-30: {order: 1, stereo_parity: null}
          7-18: {order: 1, stereo_parity: null}
          7-31: {order: 1, stereo_parity: null}
          7-32: {order: 1, stereo_parity: null}
          8-9: {order: 1, stereo_parity: null}
          8-11: {order: 1, stereo_parity: null}
          8-33: {order: 1, stereo_parity: null}
          8-34: {order: 1, stereo_parity: null}
          9-14: {order: 1, stereo_parity: null}
          10-11: {order: 1, stereo_parity: null}
          10-17: {order: 1, stereo_parity: null}
          10-35: {order: 1, stereo_parity: null}
          11-18: {order: 1, stereo_parity: null}
          11-36: {order: 1, stereo_parity: null}
          12-15: {order: 1, stereo_parity: null}
          12-19: {order: 1, stereo_parity: null}
          13-16: {order: 1, stereo_parity: null}
          13-20: {order: 1, stereo_parity: null}
          13-37: {order: 1, stereo_parity: null}
          14-15: {order: 1, stereo_parity: null}
          14-17: {order: 1, stereo_parity: null}
          15-21: {order: 1, stereo_parity: null}
          16-17: {order: 1, stereo_parity: null}
          16-21: {order: 1, stereo_parity: null}
          16-38: {order: 1, stereo_parity: null}
          19-39: {order: 1, stereo_parity: null}
          20-40: {order: 1, stereo_parity: null}
    """)
    KEYS = sorted(automol.graph.atom_keys(GRA))
    P_DCT = planarity_constraint_bounds(GRA, KEYS)
    print(P_DCT)
