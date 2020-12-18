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
2-6. The actual work for these steps is handled in a separate module,
automol.embed.
"""
import itertools
import numpy
import qcelemental as qce
from qcelemental import constants as qcc
from qcelemental import periodictable as pt
import automol.create.geom
from automol import embed
from automol.graph._graph_base import string
from automol.graph._graph import atom_symbols
from automol.graph._graph import atom_keys
from automol.graph._graph import atom_shortest_paths
from automol.graph._graph import atom_neighbor_keys
from automol.graph._ring import rings_atom_keys
from automol.graph._res import resonance_dominant_atom_hybridizations

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
DEG2RAD = qcc.conversion_factor('degree', 'radian')
# RAD2DEG = qcc.conversion_factor('radian', 'degree')

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

    :param path: the shortest path between two atoms
    :type path: list or tuple
    """

    rng_keys_lst = rings_atom_keys(gra)
    hyb_dct = resonance_dominant_atom_hybridizations(gra)

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


def distance_bounds_matrices(gra, keys):
    """ initial distance bounds matrices
    """
    assert set(keys) == set(atom_keys(gra))

    sp_dct = atom_shortest_paths(gra)

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


def sample_raw_distance_coordinates(gra, keys, dim4=False):
    """ sample raw (uncorrected) distance coordinates
    """
    # 1. Generate distance bounds matrices, L and U
    lmat, umat = distance_bounds_matrices(gra, keys)

    # 2. Triangle-smooth the bounds matrices
    lmat, umat = embed.triangle_smooth_bounds_matrices(lmat, umat)

    # 3. Generate a distance matrix D by sampling within the bounds
    dmat = embed.sample_distance_matrix(lmat, umat)

    # 4. Generate the metric matrix G
    gmat = embed.metric_matrix(dmat)

    # 5-6. Generate coordinates from the metric matrix
    xmat = embed.coordinates_from_metric_matrix(gmat, dim4=dim4)

    return xmat


def sample_raw_distance_geometry(gra, keys):
    """ sample a raw (uncorrected) distance geometry
    """
    xmat = sample_raw_distance_coordinates(gra, keys, dim4=False)

    sym_dct = atom_symbols(gra)
    syms = list(map(sym_dct.__getitem__, keys))
    geo = automol.create.geom.from_data(syms, xmat, angstrom=True)

    return geo


if __name__ == '__main__':
    import automol
    ICH = automol.smiles.inchi('CO')

    GRA = automol.inchi.graph(ICH)
    GRA = automol.graph.explicit(GRA)
    KEYS = sorted(automol.graph.atom_keys(GRA))

    GEO = sample_raw_distance_geometry(GRA, KEYS)
    print(automol.geom.string(GEO))
