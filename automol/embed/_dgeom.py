"""implements the basic distance geometry algoirthm.

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

This module handles the numerical work for steps 2-6.

The function sample_raw_distance_geometry() below returns the result of step
2-6. The actual work for these steps is handled in a separate module,
automol.embed.

Step 1 is performed based on the molecular graph and is a submodule of
automol.graph.

Step 7 is a whole separate algorithm which is undertaken in a separate file in
this module.
"""
import itertools

import numpy
from _collections_abc import Sequence

Sequence2D = Sequence[Sequence[float]]

MatrixLike = numpy.ndarray | Sequence2D


def sample_raw_distance_coordinates(
    lmat: MatrixLike, umat: MatrixLike, dim4=True
) -> MatrixLike:
    """Sample raw (uncorrected) distance coordinates.
    :param lmat: Lower-bound distance matrix
    :param umat: Upper-bound distance matrix.
    :return: Matrix with raw distance coordinates.
    """
    # 2. Triangle-smooth the bounds matrices
    lmat, umat = triangle_smooth_bounds_matrices(lmat, umat)

    # 3. Generate a distance matrix D by sampling within the bounds
    dmat = sample_distance_matrix(lmat, umat)

    # 4. Generate the metric matrix G
    gmat = metric_matrix(dmat)

    # 5-6. Generate coordinates from the metric matrix
    xmat = coordinates_from_metric_matrix(gmat, dim4=dim4)

    return xmat


def triangle_smooth_bounds_matrices(lmat: MatrixLike, umat: MatrixLike) -> MatrixLike:
    """Smoothing of the bounds matrix by triangle inequality.

    Dress, A. W. M.; Havel, T. F. "Shortest-Path Problems and Molecular
    Conformation"; Discrete Applied Mathematics (1988) 19 p. 129-144.

    This algorithm is directly from p. 8 in the paper.
    :param lmat: Lower-bound distance matrix
    :param umat: Upper-bound distance matrix
    :return: Smooth matrix
    """
    lmat, umat = map(numpy.array, (lmat, umat))

    natms = len(umat)

    for k in range(natms):
        for i, j in itertools.combinations(range(natms), 2):
            if umat[i, j] > umat[i, k] + umat[k, j]:
                umat[i, j] = umat[i, k] + umat[k, j]
            if lmat[i, j] < lmat[i, k] - umat[k, j]:
                lmat[i, j] = lmat[i, k] - umat[k, j]
            if lmat[i, j] < lmat[j, k] - umat[k, i]:
                lmat[i, j] = lmat[j, k] - umat[k, i]

            if numpy.allclose(lmat[i, j], umat[i, j]):
                lmat[i, j] = umat[i, j] = numpy.round(lmat[i, j], 6)

            assert lmat[i, j] <= umat[i, j], (
                "Lower bound exceeds upper bound. Something is wrong!",
                lmat[i, j],
                ">",
                umat[i, j],
            )

    return lmat, umat


def sample_distance_matrix(lmat: MatrixLike, umat: MatrixLike) -> MatrixLike:
    """Determine a random distance matrix based on the bounds matrices.

    That is, a random guess at d_ij = |r_i - r_j|
    :param lmat: Lower-bound distance matrix
    :param umat: Upper-bound distance matrix
    :return: Sample matrix
    """
    lmat, umat = map(numpy.array, (lmat, umat))

    dmat = numpy.random.uniform(lmat, umat)
    # The sampling will not come out symmetric, so replace the lower triangle
    # with upper triangle values
    tril = numpy.tril_indices_from(dmat)
    dmat[tril] = (dmat.T)[tril]
    return dmat


def distances_from_center(dmat: MatrixLike) -> MatrixLike:
    """Get the vector of distances from the center (average position).

    The "center" in this case is the average of the position vectors. The
    elements of this vector are therefore dc_i = |r_c - r_i|.

    See the following paper for a derivation of this formula.

    Crippen, G. M.; Havel, T. F. "Stable Calculation of Coordinates from
    Distance Information"; Acta Cryst. (1978) A34 p. 282-284

    I verified this against the alternative formula:
        dc_i^2 = 1/(2n^2) sum_j sum_k (d_ij^2 + d_ik^2 - d_jk^2)
    from page 284 of the paper.
    :dmat: Distance Matrix
    :return: Distances from center
    """
    dmat = numpy.array(dmat)

    natms = len(dmat)

    dcvec = numpy.zeros((natms,))

    for i in range(natms):
        sum_dij2 = sum(dmat[i, j] ** 2 for j in range(natms))
        sum_djk2 = sum(
            dmat[j, k] ** 2 for j, k in itertools.combinations(range(natms), 2)
        )

        dci2 = numpy.abs(sum_dij2 / natms - sum_djk2 / (natms**2))

        dcvec[i] = numpy.sqrt(dci2)

    return dcvec


def metric_matrix(dmat: MatrixLike) -> MatrixLike:
    """Compute the matrix of position vector dot products, with a central origin.

    "Central" in this case mean the average of the position vectors. So these
    elements are g_ij = (r_i - r_c).(r_j - r_c) where r_c is the average
    position (or center of mass, if all atoms have the same mass).

    See the following paper for a derivation of this formula.

    Crippen, G. M.; Havel, T. F. "Stable Calculation of Coordinates from
    Distance Information"; Acta Cryst. (1978) A34 p. 282-284
    """
    dmat = numpy.array(dmat)

    natms = len(dmat)

    dcvec = distances_from_center(dmat)

    gmat = numpy.eye(natms)

    for i, j in itertools.product(range(natms), range(natms)):
        gmat[i, j] = (dcvec[i] ** 2 + dcvec[j] ** 2 - dmat[i, j] ** 2) / 2.0

    return gmat


def coordinates_from_metric_matrix(gmat: MatrixLike, dim4: bool = False) -> MatrixLike:
    """Determine molecule coordinates from the metric matrix."""
    gmat = numpy.array(gmat)

    dim = 3 if not dim4 else 4

    vals, vecs = numpy.linalg.eigh(gmat)
    vals = vals[::-1]
    vecs = vecs[:, ::-1]
    vals = vals[:dim]
    vecs = vecs[:, :dim]
    lvec = numpy.sqrt(numpy.abs(vals))

    xmat = vecs @ numpy.diag(lvec)

    return xmat


def metric_matrix_from_coordinates(xmat: MatrixLike) -> MatrixLike:
    """Determine the metric matrix from coordinates.

    (for testing purposes only!)
    """
    xmat = numpy.array(xmat)
    return xmat @ xmat.T


def distance_matrix_from_coordinates(xmat: MatrixLike, dim4: bool = True) -> MatrixLike:
    """Determine the distance matrix from coordinates.

    (for testing purposes only!)
    """
    xmat = numpy.array(xmat)
    if not dim4:
        xmat = xmat[:, :3]

    natms = len(xmat)

    dmat = numpy.zeros((natms, natms))

    for i, j in itertools.product(range(natms), range(natms)):
        dmat[i, j] = numpy.linalg.norm(xmat[i] - xmat[j])

    return dmat


def greatest_distance_errors(
    dmat: MatrixLike, lmat: MatrixLike, umat: MatrixLike, count: int = 10
):
    """Get the indices of the maximum distance errors.
    :dmat: Distance matrix
    :param lmat: Lower-bound distance matrix
    :param umat: Upper-bound distance matrix
    return: Dictionary of greatest distance errores.
    """
    lerrs = (lmat - dmat) * (lmat >= dmat)
    uerrs = (dmat - umat) * (dmat >= umat)
    errs = numpy.maximum(lerrs, uerrs)
    idx_vecs = numpy.unravel_index(numpy.argsort(-errs, axis=None), errs.shape)
    vals = errs[idx_vecs][:count]
    idxs = tuple(map(tuple, zip(*idx_vecs, strict=True)))[:count]
    err_dct = dict(zip(idxs, vals, strict=True))
    return err_dct
