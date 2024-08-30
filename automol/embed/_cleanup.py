"""Implements geometry purification for the distance geometry algorithm.

This is used to clean up the structure of the molecule and enforce correct
chirality, by minimizing an error function.

The error function is minimized using the conjugate gradients method:

    1. Calculate the steepest direction: dn = - grad(Err(xn))
    2. Calculate coeff Bn as Bn = min(0, Bn^PR)

        Bn^PR = dn.(dn - dn-1) / dn-1.dn-1

    3. Update the conjugate direction: sn = dn + Bn*sn-1
    4. Perform a line search to optimize An = argmin Err(xn + An*sn)
    5. Update the position: xn+1 = xn + An*sn
    6. Recalculate the gradient. If the gradient is small, quit. Otherwise,
    return to step 1.

Style convention: Functions returning callables (i.e. functions returning
functions) and callable variables are marked by trailing underscores.

See Havel, T. F.; "Distance Geometry: Theory, Algorithms, and Chemical
Applications"; Encyclopedia of Computational Chemistry (2002) for details.

The error function formulas are found on page 11 of this paper. The gradient
formulas I derived myself. On pages 12-13, the use of four-dimensional
coordinates to improve convergence is described.
"""
import logging

import numpy
import scipy.optimize
from _collections_abc import Callable, Sequence

from ._dgeom import (
    distance_matrix_from_coordinates,
    greatest_distance_errors,
)
from ._findif import central_difference

SignedVolumeContraints = dict[tuple[int, int, int, int] : tuple[float, float]]


Sequence2D = Sequence[Sequence[float]]
NDArrayLike2D = numpy.ndarray | Sequence2D

# uncomment this line if you want the logging statements while debugging:
# logging.basicConfig(format='%(message)s', level=logging.DEBUG)

X = numpy.newaxis


def volume(xmat: NDArrayLike2D, idxs: list) -> float:
    """Calculate signed tetrahedral volume for a tetrad of atoms.

    for a tetrad of four atoms (1, 2, 3, 4) around a central atom, the signed
    volume formula of this tetrahedral pyramid is given by

        d12 . (d13 x d14)

    where dij = rj - ri, . is the dot product, and x is the cross product

    :param xmat: The matrix of XYZ coordinates
    :param idxs: A tetrad of four atom indices defining the volume
    :return: The volume
    """
    xmat = numpy.array(xmat)
    idxs = list(idxs)
    xyzs = xmat[:, :3][idxs]
    d12 = xyzs[1] - xyzs[0]
    d13 = xyzs[2] - xyzs[0]
    d14 = xyzs[3] - xyzs[0]
    # Negate the sign to match the way we calculate this elsewhere
    vol = numpy.negative(numpy.dot(d12, numpy.cross(d13, d14)))
    return vol


def volume_gradient(xmat: NDArrayLike2D, idxs: list) -> numpy.ndarray:
    """Calculate the tetrahedral volume gradient for a tetrad of atoms.
    :param xmat: The matrix of XYZ coordinates
    :param idxs: A tetrad of four atom indices defining the volume
    :return: The volume gradient.
    """
    xmat = numpy.array(xmat)
    idxs = list(idxs)
    xyzs = xmat[:, :3][idxs]

    grad = numpy.zeros_like(xmat)
    grad[idxs[0], :3] = numpy.cross(xyzs[1], xyzs[3] - xyzs[2]) - numpy.cross(
        xyzs[2], xyzs[3]
    )
    grad[idxs[1], :3] = +numpy.cross(xyzs[2] - xyzs[0], xyzs[3] - xyzs[0])
    grad[idxs[2], :3] = -numpy.cross(xyzs[1] - xyzs[0], xyzs[3] - xyzs[0])
    grad[idxs[3], :3] = +numpy.cross(xyzs[1] - xyzs[0], xyzs[2] - xyzs[0])

    # Negate the sign to match the way we calculate this elsewhere
    return numpy.negative(grad)


def error_function_(
    lmat: NDArrayLike2D,
    umat: NDArrayLike2D,
    chi_dct: SignedVolumeContraints = None,
    pla_dct: SignedVolumeContraints = None,
    wdist: float = 1.0,
    wchip: float = 1.0,
    wdim4: float = 1.0,
    leps: float = 0.1,
    ueps: float = 0.1,
    log: bool = False,
) -> Callable[[NDArrayLike2D], float]:
    """Compute the embedding error function.

    :param lmat: Lower-bound distance matrix
    :param umat: Upper-bound distance matrix
    :param chi_dct: Chirality constraints; the keys are tuples of four atoms,
        the values are lower and upper bounds on the four-point signed volume
        of these atoms
    :param pla_dct: Planarity constraints; the keys are tuples of four atoms,
        the values are lower and upper bounds on the four-point signed volume
        of these atoms
    :param wdist: Weight on the distance constraint
    :param wchip: Weight on the chirality/planarity constraint
    :param wdim4: Weight on the fourth dimension constraint
    :param leps: Denominator epsilon for lower bound distances
    :param ueps: Denominator epsilon for upper bound distances
    """
    triu = numpy.triu_indices_from(lmat)
    chi_dct = {} if chi_dct is None else chi_dct
    pla_dct = {} if pla_dct is None else pla_dct
    chip_dct = {**chi_dct, **pla_dct}

    def _function(xmat: NDArrayLike2D):
        dmat = distance_matrix_from_coordinates(xmat)

        # distance error (equation 61 in the paper referenced above)
        ltf = ((lmat**2 - dmat**2) / (leps**2 + dmat**2))[triu]
        utf = ((dmat**2 - umat**2) / (ueps**2 + umat**2))[triu]
        ltf *= ltf > 0.0
        utf *= utf > 0.0
        dist_err = wdist * (numpy.vdot(utf, utf) + numpy.vdot(ltf, ltf))

        if log:
            print("Error report:")
            print("\tDistance error:", dist_err)
            lerrs = (lmat - dmat) * (lmat > dmat)
            uerrs = (dmat - umat) * (dmat > umat)
            lsrt_idxs = numpy.unravel_index(
                numpy.argsort(numpy.ravel(-lerrs)), lerrs.shape
            )
            usrt_idxs = numpy.unravel_index(
                numpy.argsort(numpy.ravel(-uerrs)), uerrs.shape
            )
            print("\tGreatest lower-bound errors:")
            for idxs in list(zip(*lsrt_idxs, strict=True))[:5]:
                print("\t\t", idxs, lerrs[idxs])
            print("\tGreatest upper-bound errors:")
            for idxs in list(zip(*usrt_idxs, strict=True))[:5]:
                print("\t\t", idxs, uerrs[idxs])

        # chirality/planarity error (equation 62 in the paper referenced above)
        if chip_dct:
            vols = numpy.array([volume(xmat, idxs) for idxs in chip_dct.keys()])
            lvols, uvols = map(numpy.array, zip(*chip_dct.values(), strict=True))
            ltv = (lvols - vols) * (vols < lvols)
            utv = (vols - uvols) * (vols > uvols)
            chip_err = wchip * (numpy.vdot(ltv, ltv) + numpy.vdot(utv, utv))
        else:
            chip_err = 0.0

        if log:
            print("\tChirality/planarity error:", chip_err)

        # fourth-dimension error
        if numpy.shape(xmat)[1] == 4:
            dim4_err = wdim4 * numpy.vdot(xmat[:, 3], xmat[:, 3])
        else:
            dim4_err = 0.0

        if log:
            print("\tFourth dimension error:", dim4_err)
            print()

        return dist_err + chip_err + dim4_err

    return _function


def error_function_gradient_(
    lmat: NDArrayLike2D,
    umat: NDArrayLike2D,
    chi_dct: SignedVolumeContraints = None,
    pla_dct: SignedVolumeContraints = None,
    wdist: float = 1.0,
    wchip: float = 1.0,
    wdim4: float = 1.0,
    leps: float = 0.1,
    ueps: float = 0.1,
) -> Callable[[NDArrayLike2D], numpy.ndarray]:
    """Check the embedding error function gradient.

    :param lmat: Lower-bound distance matrix
    :param umat: Upper-bound distance matrix
    :param chi_dct: Chirality constraints; the keys are tuples of four atoms,
        the values are lower and upper bounds on the four-point signed volume
        of these atoms
    :param pla_dct: Planarity constraints; the keys are tuples of four atoms,
        the values are lower and upper bounds on the four-point signed volume
        of these atoms
    :param wdist: Weight on the distance constraint
    :param wchip: Weight on the chirality/planarity constraint
    :param wdim4: Weight on the fourth dimension constraint
    :param leps: Denominator epsilon for lower bound distances
    :param ueps: Denominator epsilon for upper bound distances
    :return: Error function gradient
    """
    chi_dct = {} if chi_dct is None else chi_dct
    pla_dct = {} if pla_dct is None else pla_dct
    chip_dct = {**chi_dct, **pla_dct}

    def _gradient(xmat):
        dmat = distance_matrix_from_coordinates(xmat)

        # distance error gradient
        utf = (dmat**2 - umat**2) / (ueps**2 + umat**2)
        ltf = (lmat**2 - dmat**2) / (leps**2 + dmat**2)
        utg = (+4.0 * utf / (ueps**2 + umat**2)) * (utf > 0.0)
        ltg = (
            (-4.0 * ltf / (leps**2 + dmat**2) ** 2)
            * (leps**2 + lmat**2)
            * (ltf > 0.0)
        )
        utg = utg[:, :, X]
        ltg = ltg[:, :, X]
        xmx = xmat[:, X, :] - xmat[X, :, :]
        dist_grad = numpy.sum(xmx * (ltg + utg), axis=1)
        dist_grad *= wdist

        # chirality/planarity error gradient
        if chip_dct:
            vols = numpy.array([volume(xmat, idxs) for idxs in chip_dct.keys()])
            vol_grads = numpy.array(
                [volume_gradient(xmat, idxs) for idxs in chip_dct.keys()]
            )
            lvols, uvols = map(numpy.array, zip(*chip_dct.values(), strict=True))
            ltv = (lvols - vols) * (vols < lvols)
            utv = (vols - uvols) * (vols > uvols)
            ltg = -2.0 * ltv[:, X, X] * vol_grads
            utg = +2.0 * utv[:, X, X] * vol_grads
            chip_grad = numpy.sum(ltg + utg, axis=0)
            chip_grad *= wchip
        else:
            chip_grad = numpy.zeros_like(xmat)

        # fourth-dimension error gradient
        if numpy.shape(xmat)[1] == 4:
            dim3_zeros = numpy.zeros_like(xmat[:, :3])
            dim4_grad = numpy.hstack([dim3_zeros, 2 * xmat[:, 3:]])
            dim4_grad *= wdim4
        else:
            dim4_grad = numpy.zeros_like(xmat)

        return dist_grad + chip_grad + dim4_grad

    return _gradient


def error_function_numerical_gradient_(
    lmat: NDArrayLike2D,
    umat: NDArrayLike2D,
    chi_dct: SignedVolumeContraints = None,
    pla_dct: SignedVolumeContraints = None,
    wdist: float = 1.0,
    wchip: float = 1.0,
    wdim4: float = 1.0,
    leps: float = 0.1,
    ueps: float = 0.1,
) -> Callable[[NDArrayLike2D], numpy.ndarray]:
    """Check the gradient of the distance error function.

    (For testing purposes only; Used to check the analytic gradient formula.)
    :param lmat: Lower-bound distance matrix
    :param umat: Upper-bound distance matrix
    :param chi_dct: Chirality constraints; the keys are tuples of four atoms,
        the values are lower and upper bounds on the four-point signed volume
        of these atoms
    :param pla_dct: Planarity constraints; the keys are tuples of four atoms,
        the values are lower and upper bounds on the four-point signed volume
        of these atoms
    :param wdist: Weight on the distance constraint
    :param wchip: Weight on the chirality/planarity constraint
    :param wdim4: Weight on the fourth dimension constraint
    :param leps: Denominator epsilon for lower bound distances
    :param ueps: Denominator epsilon for upper bound distances
    :return:Gradient of the distance error function
    """
    erf_ = error_function_(
        lmat,
        umat,
        chi_dct=chi_dct,
        pla_dct=pla_dct,
        wdist=wdist,
        wchip=wchip,
        wdim4=wdim4,
        leps=leps,
        ueps=ueps,
    )

    def _gradient(xmat):
        grad = central_difference(erf_, xmat, npts=11)
        return grad

    return _gradient


def polak_ribiere_beta(sd1: numpy.ndarray, sd0: numpy.ndarray) -> float:
    """Determine the conjugate gradient alpha coefficient
    from an error-minimizing line search.

    (See https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method)

    :param sd1: The steepest-descent (negative gradient) direction of the current step
    :return: The alpha  coefficient
    """
    return numpy.vdot(sd1, sd1 - sd0) / numpy.vdot(sd0, sd0)


def line_search_alpha(
    err_: Callable[[NDArrayLike2D], float], sd1: numpy.ndarray, cd1: numpy.ndarray
):
    """Perform a line search to determine the alpha coefficient."""

    # define the objective function
    def _function_of_alpha(alpha):
        return err_(sd1 + alpha * cd1)

    # do the line search and make sure it worked
    res = scipy.optimize.minimize_scalar(_function_of_alpha)
    if not res.success:
        res = scipy.optimize.minimize_scalar(_function_of_alpha, bounds=(0, 20))

    assert res.success, ("Line search for alpha failed!\n", str(res))

    # get the result
    alpha = res.x

    return alpha


def cleaned_up_coordinates(
    xmat: NDArrayLike2D,
    lmat: NDArrayLike2D,
    umat: NDArrayLike2D,
    chi_dct: SignedVolumeContraints = None,
    pla_dct: SignedVolumeContraints = None,
    conv_=None,
    max_dist_err: float = 0.2,
    grad_thresh: float = 0.2,
    maxiter: int | None = None,
    chi_flip: bool = True,
    dim4: bool = True,
    log: bool = False,
):
    """Clean up coordinates by conjugate-gradients error minimization.

    :param xmat: the initial guess coordinates to be cleaned up
    :param lmat: lower-bound distance matrix
    :param umat: upper-bound distance matrix
    :param chi_dct: chirality constraints; the keys are tuples of four atoms,
        the values are lower and upper bounds on the four-point signed volume
        of these atoms
    :param pla_dct: planarity constraints; the keys are tuples of four atoms,
        the values are lower and upper bounds on the four-point signed volume
        of these atoms
    :param conv_: a callable convergence checker function of xmat, err, and
        grad which returns True if the geometry is converged
    :param max_dist_err: maximum distance error for the geometry to be
        considered converged
    :param grad_thresh: maximum gradient norm for the geometry to be
        considered converged
    :param maxiter: maximum number of iterations; default is three times the
        number of coordinates
    :param chi_flip: whether or not to invert the structure if more than half
        of the chiralities are reversed
    :param dim4: whether or not to include a fourth dimension, for allowing
        chiralities to flip as they correct themselves
    """
    xmat = numpy.array(xmat)

    # Make the coordinates four-dimensional, if they aren't already
    natms, ndims = numpy.shape(xmat)
    if ndims < 4 and dim4:
        xmat = numpy.hstack([xmat, numpy.zeros((natms, 4 - ndims))])

    # If less than half of the chiralities have correct sign, invert the
    # geometry
    if chi_flip and chi_dct:
        current_vols = numpy.array([volume(xmat, idxs) for idxs in chi_dct.keys()])
        target_vols = numpy.array(list(map(numpy.average, chi_dct.values())))
        comparison = numpy.sign(current_vols) == numpy.sign(target_vols)
        fraction = numpy.average(comparison)
        if fraction < 0.5:
            xmat *= -1.0

    maxiter = int(numpy.size(xmat) * 5 if maxiter is None else maxiter)

    err_ = error_function_(
        lmat, umat, chi_dct=chi_dct, pla_dct=pla_dct, wdim4=1.0, log=log
    )
    grad_ = error_function_gradient_(
        lmat, umat, chi_dct=chi_dct, pla_dct=pla_dct, wdim4=1.0
    )
    conv_ = (
        default_convergence_checker_(lmat, umat, max_dist_err, grad_thresh)
        if conv_ is None
        else conv_
    )

    xmat, conv = minimize_error(xmat, err_, grad_, conv_, maxiter)
    return xmat, conv


def default_convergence_checker_(
    lmat: NDArrayLike2D, umat: NDArrayLike2D, max_dist_err=0.2, grad_thresh=0.2
) -> Callable[[NDArrayLike2D, float, NDArrayLike2D], bool]:
    """Check for default convergence."""
    conv1_ = distance_convergence_checker_(lmat, umat, max_dist_err)
    conv2_ = gradient_convergence_checker_(grad_thresh)

    def _is_converged(xmat, err, grad):
        return conv1_(xmat, err, grad) and conv2_(xmat, err, grad)

    return _is_converged


def distance_convergence_checker_(
    lmat: NDArrayLike2D, umat: NDArrayLike2D, max_dist_err: float = 0.2
) -> numpy.ndarray:
    """Convergence checker based on the maximum distance error."""

    def _is_converged(xmat, err, grad):
        assert err or not err
        assert numpy.shape(xmat) == numpy.shape(grad)
        dmat = distance_matrix_from_coordinates(xmat)
        dist_err_dct = greatest_distance_errors(dmat, lmat, umat, count=1)
        (dist_err,) = dist_err_dct.values()
        return dist_err <= max_dist_err

    return _is_converged


def planarity_convergence_checker_(
    pla_dct: dict, max_vol_err: float = 0.2
) -> Callable[[NDArrayLike2D, float, NDArrayLike2D], bool]:
    """Convergence checker based on the maximum planarity error."""

    def _is_converged(xmat, err, grad):
        assert err or not err
        assert numpy.shape(xmat) == numpy.shape(grad)
        vols = numpy.array([volume(xmat, idxs) for idxs in pla_dct.keys()])
        lvols, uvols = map(numpy.array, zip(*pla_dct.values(), strict=True))
        lmax = numpy.amax((lvols - vols) * (vols < lvols))
        umax = numpy.amax((vols - uvols) * (vols > uvols))
        return max(lmax, umax) <= max_vol_err

    return _is_converged


def gradient_convergence_checker_(
    thresh: float = 1e-1,
) -> Callable[[NDArrayLike2D, float, NDArrayLike2D], bool]:
    """Maximum gradient convergence checker."""

    def _is_converged(xmat, err, grad):
        assert numpy.shape(xmat) == numpy.shape(grad)
        grad_max = numpy.amax(numpy.abs(grad))
        logging.info(f"\tError: {err:f}")
        logging.info(f"\tMax gradient: {grad_max:f}")
        logging.info("\n")
        return grad_max < thresh

    return _is_converged


def minimize_error(
    xmat: NDArrayLike2D,
    err_: Callable[[NDArrayLike2D], float],
    grad_,
    conv_,
    maxiter: int | None = None,
) -> bool:
    """Do conjugate-gradients error minimization.

    :param err_: a callable error function of xmat
    :param grad_: a callable error gradient function of xmat
    :param conv_: a callable convergence checker function of xmat, err_(xmat),
        and grad_(xmat) which returns True if the geometry is converged
    :param maxiter: maximum number of iterations; default is three times the
        number of coordinates
    :returns: the optimized coordinates and a boolean which is True if
        converged and False if not
    """
    maxiter = numpy.size(xmat) * 5 if maxiter is None else maxiter

    sd0 = None
    cd0 = None
    logging.info(f"Initial error: {err_(xmat):f}")

    converged = False

    for niter in range(maxiter):
        logging.info(f"Iteration {niter:d}")

        # 1. Calculate the steepest direction
        sd1 = -grad_(xmat)

        # 2-3. Determine the conjugate direction
        if sd0 is None:
            cd1 = sd1
        else:
            # 2. Compute beta
            beta = min(0.0, polak_ribiere_beta(sd1, sd0))

            # 3. determine step direction
            cd1 = sd1 + beta * cd0

        # 4. Perform a line search
        alpha = line_search_alpha(err_, xmat, cd1)
        logging.info(f"{niter} alpha:")
        logging.info(alpha)

        # 5. Check convergence
        if conv_(xmat, err_(xmat), sd1):
            converged = True
            break

        # 6. Take the step
        xmat += alpha * cd1

        sd0 = sd1
        cd0 = cd1

    logging.info(f"Niter: {niter:d}")
    logging.info(f"Converged: {('Yes' if converged else 'No'):s}")
    logging.info("\n")

    return xmat, converged
