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
import numpy
import scipy.optimize
from automol.embed._dgeom import distance_matrix_from_coordinates

X = numpy.newaxis


def error_function_(lmat, umat, wdim4=1., leps=0.1, ueps=0.1):
    """ the embedding error function
    """
    triu = numpy.triu_indices_from(lmat)

    def _function(xmat):
        dmat = distance_matrix_from_coordinates(xmat)

        utf = ((dmat**2-umat**2) / (ueps**2+umat**2))[triu]
        ltf = ((lmat**2-dmat**2) / (leps**2+dmat**2))[triu]

        utf *= (utf > 0.)
        ltf *= (ltf > 0.)

        dist_err = numpy.vdot(utf, utf) + numpy.vdot(ltf, ltf)
        dim4_err = wdim4 * numpy.vdot(xmat[:, 3], xmat[:, 3])

        return dist_err + dim4_err

    return _function


def error_function_gradient_(lmat, umat, wdim4=1., leps=0.1, ueps=0.1):
    """ the embedding error function gradient
    """
    # natms = len(lmat)

    def _gradient(xmat):
        dmat = distance_matrix_from_coordinates(xmat)

        utf = (dmat**2-umat**2) / (ueps**2+umat**2)
        ltf = (lmat**2-dmat**2) / (leps**2+dmat**2)
        utg = (+4.*utf/(ueps**2+umat**2))*(utf > 0.)
        ltg = (-4.*ltf/(leps**2+dmat**2)**2)*(leps**2+lmat**2)*(ltf > 0.)
        utg = utg[:, :, X]
        ltg = ltg[:, :, X]
        xmx = xmat[:, X, :] - xmat[X, :, :]

        dist_grad = numpy.sum(xmx*(ltg + utg), axis=1)
        dim4_grad = numpy.hstack([numpy.zeros_like(xmat[:, :3]),
                                  2*wdim4*xmat[:, 3:]])

        return dist_grad + dim4_grad

    return _gradient


def error_function_numerical_gradient_(lmat, umat, wdim4=1.):
    """ the gradient of the distance error function

    (For testing purposes only; Used to check the analytic gradient formula.)
    """

    erf_ = error_function_(lmat, umat, wdim4=wdim4)

    def _gradient(xmat):
        grad = central_difference(erf_, xmat, npts=11)
        return grad

    return _gradient


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


def cleaned_up_coordinates(xmat, lmat, umat, thresh=1e-1, maxiter=None):
    """ clean up coordinates by conjugate-gradients error minimization
    """
    fun_ = error_function_(lmat, umat, wdim4=1.)
    grad_ = error_function_gradient_(lmat, umat, wdim4=1.)

    maxiter = numpy.size(xmat) * 3
    print(maxiter)

    sd0 = None
    cd0 = None
    print('function:', fun_(xmat))

    converged = False

    for niter in range(maxiter):
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

    print('niter:', niter)

    return xmat, converged


if __name__ == '__main__':
    numpy.random.seed(5)
    from automol.embed._findif import central_difference
    import automol
    ICH = automol.smiles.inchi('CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O')
    # ICH = automol.smiles.inchi('C1CCC2CC(CCC3C4CCC5CC4C53)CC2C1')

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

    GRA2 = automol.geom.connectivity_graph(GEO)
    print(GRA == GRA2)
