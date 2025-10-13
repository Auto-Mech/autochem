""" evaluate finite difference derivatives to arbitrary order
"""
import numpy
import scipy.misc


def central_difference(fun, arg, step=0.01, nder=1, npts=None):
    """differentiate a function using central differences
    :param fun: the function
    :type fun: typing.Callable
    :param arg: the point at which to evaluate the derivative
    :type arg: float or numpy.ndarrray
    :param step: step size, or a grid of step sizes corresponding to `x`
    :type step: float or numpy.ndarray
    :param nder: return the nth derivative
    :type nder: int
    :param npts: the number of grid points, default `nder` + `1` + `nder % 2`
    :type npts: int
    """
    if npts is None:
        npts = nder + 1 + nder % 2
    if numpy.ndim(step) == 0:
        step = float(step) * numpy.ones_like(arg)

    weights = scipy.misc.central_diff_weights(Np=npts, ndiv=nder)

    def derivative(index):
        darg = numpy.zeros_like(arg)
        darg[index] = step[index]
        grid = [numpy.array(arg) + (k - npts//2) * darg for k in range(npts)]
        vals = map(fun, grid)
        return (sum(numpy.multiply(w, v) for w, v in zip(weights, vals))
                / (step[index] ** nder))

    der = tuple(map(derivative, numpy.ndindex(numpy.shape(arg))))
    shape = numpy.shape(arg) + numpy.shape(fun(arg))
    return numpy.reshape(der, shape)
