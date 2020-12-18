""" implements geometry purification for the distance geometry algorithm

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
"""
import itertools
import numpy
import scipy.optimize
from automol.embed._dgeom import distance_matrix_from_coordinates


def distance_error_function_(lmat, umat):
    """ the distance error function
    """

    def _error_function(xmat):
        natms = len(xmat)

        dmat = distance_matrix_from_coordinates(xmat)

        err = sum(
            max(0, ((dmat[i, j]/umat[i, j])**2 - 1)**2) +
            max(0, (2*lmat[i, j]**2/(lmat[i, j]**2 + dmat[i, j]**2) - 1)**2)
            for i, j in itertools.combinations(range(natms), 2))

        return err

    return _error_function


def fourth_dimension_error_function_(weight):
    """ the fourth-dimention error function
    """

    def _error_function(xmat):
        err = weight * numpy.linalg.norm(xmat[:, 3]) ** 2
        return err

    return _error_function


def distance_error_function_gradient_(lmat, umat):
    """ the gradient of the distance error function

    (Currently determined numerically. We can work out an analytic formula when
    the time comes...)
    """

    erf_ = distance_error_function_(lmat, umat)

    def _error_function_gradient(xmat):
        grad = central_difference(erf_, xmat)
        return grad

    return _error_function_gradient


def fourth_dimension_error_function_gradient_(weight):
    """ the gradient of the fourth-dimension error function

    (Currently determined numerically. We can work out an analytic formula when
    the time comes...)
    """

    erf_ = fourth_dimension_error_function_(weight)

    def _error_function_gradient(xmat):
        grad = central_difference(erf_, xmat)
        return grad

    return _error_function_gradient


def error_function_(lmat, umat, weight_4d=1.):
    """ the total error function
    """
    erf1_ = distance_error_function_(lmat, umat)
    erf2_ = fourth_dimension_error_function_(weight=weight_4d)

    def _error_function(xmat):
        return erf1_(xmat) + erf2_(xmat)

    return _error_function


def error_function_gradient_(lmat, umat, weight_4d=1.):
    """ the total error function
    """
    grad1_ = distance_error_function_gradient_(lmat, umat)
    grad2_ = fourth_dimension_error_function_gradient_(weight=weight_4d)

    def _error_function_gradient(xmat):
        return grad1_(xmat) + grad2_(xmat)

    return _error_function_gradient


def polak_ribiere_beta(sd1, sd0):
    """ calculate the Polak-Ribiere Beta coefficient
    """
    return numpy.vdot(sd1, sd1-sd0) / numpy.vdot(sd0, sd0)


def line_search_alpha(fun_, sd1, cd1):
    """ perform a line search to determine the alpha coefficient
    """

    # define the objective function
    def _function_of_alpha(alpha):
        return fun_(sd1 + alpha*cd1)

    # do the line search and make sure it worked
    res = scipy.optimize.minimize_scalar(_function_of_alpha)
    assert res.success, ("Line search for alpha failed!\n", str(res))

    # get the result
    alpha = res.x

    return alpha


def cleaned_up_coordinates(xmat, lmat, umat, thresh=1e-2, maxiter=None):
    """ clean up coordinates by conjugate-gradients error minimization
    """
    fun_ = error_function_(lmat, umat, weight_4d=1.)
    grad_ = error_function_gradient_(lmat, umat, weight_4d=1.)

    maxiter = numpy.size(xmat) * 3
    print(maxiter)

    sd0 = None
    cd0 = None
    print('function:', fun_(xmat))

    converged = False

    for _ in range(100):
        # 1. Calculate the steepest direction
        sd1 = -grad_(xmat)

        # 2-3. Determine the conjugate direction
        if sd0 is None:
            cd1 = sd1
        else:
            # 2. Cumpute beta
            beta = min(0., polak_ribiere_beta(sd1, sd0))

            # 3. determine step direction
            cd1 = sd1 + beta * cd0

        # 4. Perform a line search
        alpha = line_search_alpha(fun_, xmat, cd1)

        # 5. Take the step
        xmat += alpha*cd1

        sd0 = sd1
        cd0 = cd1

        grad_max = numpy.amax(numpy.abs(sd1))
        if grad_max < thresh:
            converged = True
            break

        print('function:', fun_(xmat))
        print('gradient:', grad_(xmat))

    return xmat, converged


if __name__ == '__main__':
    from automol.embed._findif import central_difference
    import automol
    ICH = automol.smiles.inchi('C1CCC2CC(CCC3C4CCC5CC4C53)CC2C1')

    # 1. Generate distance bounds matrices, L and U
    GRA = automol.inchi.graph(ICH)
    GRA = automol.graph.explicit(GRA)
    KEYS = sorted(automol.graph.atom_keys(GRA))
    LMAT, UMAT = automol.graph.embed.distance_bounds_matrices(GRA, KEYS)
    XMAT = automol.graph.embed.sample_raw_distance_coordinates(GRA, KEYS,
                                                               dim4=True)

    XMAT, CONVERGED = cleaned_up_coordinates(XMAT, LMAT, UMAT)
    print(CONVERGED)

    SYMS = list(map(automol.graph.atom_symbols(GRA).__getitem__, KEYS))
    XYZS = XMAT[:, :3]
    GEO = automol.create.geom.from_data(SYMS, XYZS, angstrom=True)
    print(automol.geom.string(GEO))
