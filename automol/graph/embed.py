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
    7. (Not yet implemented) Do error refinement to clean up the structure and
    enforce correct chirality.
"""
import itertools
import more_itertools as mit
import numpy
import qcelemental as qce
from qcelemental import constants as qcc
from qcelemental import periodictable as pt
from automol.graph._graph_base import string
from automol.graph._graph import atom_symbols
from automol.graph._graph import atom_keys
from automol.graph._graph import atom_shortest_paths
from automol.graph._graph import atom_neighbor_keys
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


def heuristic_bond_distance(gra, key1, key2, angstrom=True, check=True):
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


def heuristic_bond_angle(gra, key1, key2, key3, degree=False, check=True):
    """ heuristic bond angle
    """
    if check:
        assert {key1, key3} <= set(atom_neighbor_keys(gra)[key2])

    hyb2 = resonance_dominant_atom_hybridizations(gra)[key2]
    if hyb2 == 3:
        ang = TET_ANG
    elif hyb2 == 2:
        ang = TRI_ANG
    else:
        assert hyb2 == 1
        ang = LIN_ANG

    ang *= 1 if degree else DEG2RAD

    return ang


def heuristic_bond_angle_distance(gra, key1, key2, key3, angstrom=True):
    """ heuristic distance between atoms at two ends of a bond angle

    uses the law of cosines:

        d13 = sqrt(d12^2 + d23^2 - 2*d12*d23*cos(a123))
    """
    d12 = heuristic_bond_distance(gra, key1, key2, angstrom=angstrom)
    d23 = heuristic_bond_distance(gra, key2, key3, angstrom=angstrom)
    a123 = heuristic_bond_angle(gra, key1, key2, key3, degree=False)
    d13 = numpy.sqrt(d12**2 + d23**2 - 2*d12*d23*numpy.cos(a123))
    return d13


def heuristic_torsion_angle_distance(gra, key1, key2, key3, key4, cis=False,
                                     angstrom=True):
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
    d12 = heuristic_bond_distance(gra, key1, key2, angstrom=angstrom)
    d23 = heuristic_bond_distance(gra, key2, key3, angstrom=angstrom)
    d34 = heuristic_bond_distance(gra, key3, key4, angstrom=angstrom)
    d13 = heuristic_bond_angle_distance(gra, key1, key2, key3,
                                        angstrom=angstrom)

    a132 = numpy.arccos((d13**2 + d23**2 - d12**2)/(2*d13*d23))
    a234 = heuristic_bond_angle(gra, key2, key3, key4, degree=False)
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


def closest_approach(gra, key1, key2):
    """ closest approach between atoms, based on their van der Waals radii

    Warning: The scaling factor on the van der waals radii was arbitrarily
    chosen based on limited tests and may need to be lowered
    """
    vdw_scaling_factor = 0.8
    dist = (van_der_waals_radius(gra, key1) +
            van_der_waals_radius(gra, key2)) * vdw_scaling_factor
    return dist


def upper_distance_bound(gra, path):
    """ upper distance bound between two ends of a path

    :param path: the shortest path between two atoms
    :type path: list or tuple
    """
    # if the path is 0, the atoms are disconnected and could be arbitrarily far
    # apart
    if len(path) == 1:
        dist = 0
    elif len(path) == 2:
        dist = heuristic_bond_distance(gra, *path)
    elif len(path) == 3:
        dist = heuristic_bond_angle_distance(gra, *path)
    elif len(path) == 4:
        dist = heuristic_torsion_angle_distance(gra, *path, cis=False)
    # otherwise, just do the sum of the distances between atoms along the path
    else:
        assert len(path) > 2
        dist = sum((heuristic_bond_distance(gra, *bond)
                    for bond in mit.windowed(path, 2)))

    return dist


def lower_distance_bound(gra, path):
    """ lower distance bound between two ends of a path

    :param path: the shortest path between two atoms
    :type path: list or tuple
    """
    if len(path) == 1:
        dist = 0
    elif len(path) == 2:
        dist = heuristic_bond_distance(gra, *path)
    elif len(path) == 3:
        dist = heuristic_bond_angle_distance(gra, *path)
    elif len(path) == 4:
        dist = heuristic_torsion_angle_distance(gra, *path, cis=True)
    # otherwise, just use the VDW radii for the two atoms
    else:
        assert len(path) > 2
        dist = closest_approach(gra, path[0], path[-1])

    return dist


def bounds_matrices(gra, keys):
    """ initial distance bounds matrices
    """
    assert set(keys) == set(atom_keys(gra))

    sp_dct = atom_shortest_paths(gra)

    natms = len(keys)
    umat = numpy.zeros((natms, natms))
    lmat = numpy.zeros((natms, natms))
    for i, j in itertools.combinations(range(natms), 2):
        if j in sp_dct[i]:
            path = sp_dct[i][j]
            umat[i, j] = umat[j, i] = upper_distance_bound(gra, path)
            lmat[i, j] = lmat[j, i] = lower_distance_bound(gra, path)
        else:
            # they are disconnected
            umat[i, j] = umat[j, i] = 999
            lmat[i, j] = lmat[j, i] = closest_approach(gra, keys[i], keys[j])

        assert lmat[i, j] <= umat[i, j], (
            "Lower bound exceeds upper bound. This is a bug!\n"
            "{}\npath: {}\n"
            .format(string(gra, one_indexed=False), str(path)))

    return lmat, umat


def triangle_smooth_bounds_matrices(lmat, umat):
    """ smoothing of the bounds matrix by triangle inequality

    Dress, A. W. M.; Havel, T. F. "Shortest-Path Problems and Molecular
    Conformation"; Discrete Applied Mathematics (1988) 19 p. 129-144.

    This algorithm is directly from p. 8 in the paper.
    """
    natms = len(umat)

    for k in range(natms):
        for i, j in itertools.combinations(range(natms), 2):
            if umat[i, j] > umat[i, k] + umat[k, j]:
                umat[i, j] = umat[i, k] + umat[k, j]
            if lmat[i, j] < lmat[i, k] - umat[k, j]:
                lmat[i, j] = lmat[i, k] - umat[k, j]
            if lmat[i, j] < lmat[j, k] - umat[k, i]:
                lmat[i, j] = lmat[j, k] - umat[k, i]

            assert lmat[i, j] <= umat[i, j], (
                "Lower bound exceeds upper bound. Something is wrong!")

    return lmat, umat


def random_distance_matrix(lmat, umat):
    """ determine a random distance matrix based on the bounds matrices

    That is, a random guess at d_ij = |r_i - r_j|
    """
    dmat = numpy.random.uniform(lmat, umat)
    # The sampling will not come out symmetric, so replace the lower triangle
    # with upper triangle values
    tril = numpy.tril_indices_from(dmat)
    dmat[tril] = (dmat.T)[tril]
    return dmat


def central_distance_vector(dmat):
    """ get the vector of distances from the center (average position)

    The "center" in this case is the average of the position vectors. The
    elements of this vector are therefore dc_i = |r_c - r_i|.

    See the following paper for a derivation of this formula.

    Crippen, G. M.; Havel, T. F. "Stable Calculation of Coordinates from
    Distance Information"; Acta Cryst. (1978) A34 p. 282-284

    I verified this against the alternative formula:
        dc_i^2 = 1/(2n^2) sum_j sum_k (d_ij^2 + d_ik^2 - d_jk^2)
    from page 284 of the paper.
    """
    natms = len(dmat)

    dcvec = numpy.zeros((natms,))

    for i in range(natms):
        sum_dij2 = sum(dmat[i, j]**2 for j in range(natms))
        sum_djk2 = sum(dmat[j, k]**2 for j, k in
                       itertools.combinations(range(natms), 2))
        dci2 = sum_dij2/natms - sum_djk2/(natms**2)

        assert dci2 > 0

        dcvec[i] = numpy.sqrt(dci2)

    return dcvec


def metric_matrix(dmat):
    """ the matrix of of position vector dot products, with a central origin

    "Central" in this case mean the average of the position vectors. So these
    elements are g_ij = (r_i - r_c).(r_j - r_c) where r_c is the average
    position (or center of mass, if all atoms have the same mass).

    See the following paper for a derivation of this formula.

    Crippen, G. M.; Havel, T. F. "Stable Calculation of Coordinates from
    Distance Information"; Acta Cryst. (1978) A34 p. 282-284
    """
    natms = len(dmat)

    dcvec = central_distance_vector(dmat)

    gmat = numpy.eye(natms)

    for i, j in itertools.product(range(natms), range(natms)):
        gmat[i, j] = (dcvec[i]**2 + dcvec[j]**2 - dmat[i, j]**2)/2.

    return gmat


def coordinates_from_metric_matrix(gmat, dim4=False):
    """ determine molecule coordinates from the metric matrix
    """
    dim = 3 if not dim4 else 4

    vals, vecs = numpy.linalg.eigh(gmat)
    vals = vals[::-1]
    vecs = vecs[:, ::-1]
    vals = vals[:dim]
    vecs = vecs[:, :dim]
    lvec = numpy.sqrt(numpy.abs(vals))

    xmat = vecs @ numpy.diag(lvec)

    return xmat


def metric_matrix_from_coordinates(xmat):
    """ determine the metric matrix from coordinates

    (for testing purposes only!)
    """
    return xmat @ xmat.T


def distance_matrix_from_coordinates(xmat):
    """ determine the distance matrix from coordinates

    (for testing purposes only!)
    """
    natms = len(xmat)

    dmat = numpy.zeros((natms, natms))

    for i, j in itertools.product(range(natms), range(natms)):
        dmat[i, j] = numpy.linalg.norm(xmat[i] - xmat[j])

    return dmat


if __name__ == '__main__':
    import automol
    ICH = automol.smiles.inchi('OO')
    GRA = automol.inchi.graph(ICH)
    GRA = automol.graph.explicit(GRA)
    KEYS = sorted(automol.graph.atom_keys(GRA))
    LMAT, UMAT = bounds_matrices(GRA, KEYS)
    LMAT, UMAT = triangle_smooth_bounds_matrices(LMAT, UMAT)
    print('lower:')
    print(numpy.round(LMAT, 2))
    print('upper:')
    print(numpy.round(UMAT, 2))
    DMAT = random_distance_matrix(LMAT, UMAT)
    print('distance:')
    print(numpy.round(DMAT, 2))
    GMAT = metric_matrix(DMAT)
    print('metric:')
    print(numpy.round(GMAT, 2))
    XMAT = coordinates_from_metric_matrix(GMAT, dim4=True)
    print(numpy.round(XMAT, 3))
    GMAT1 = metric_matrix_from_coordinates(XMAT)
    DMAT1 = distance_matrix_from_coordinates(XMAT)
    print('recalculated distance:')
    print(numpy.round(DMAT1, 2))
    print('recalculated metric:')
    print(numpy.round(GMAT1, 2))
    print(numpy.amax(GMAT1 - GMAT))
    print(numpy.amax(DMAT1 - DMAT))
    SYM_DCT = automol.graph.atom_symbols(GRA)
    SYMS = list(map(SYM_DCT.__getitem__, KEYS))
    GEO = automol.geom.from_data(SYMS, XMAT[:, :3], angstrom=True)
    print(automol.geom.string(GEO))
