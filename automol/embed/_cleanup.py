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

See Havel, T. F.; "Distance Geometry: Theory, Algorithms, and Chemical
Applications"; Encyclopedia of Computational Chemistry (2002) for details.

The error function formulas are found on page 11 of this paper. The gradient
formulas I derived myself. On pages 12-13, the use of four-dimensional
coordinates to improve convergence is described.
"""
import numpy
import scipy.optimize
<<<<<<< HEAD
=======
from automol.embed._dgeom import sample_raw_distance_coordinates
>>>>>>> Cleans up example
from automol.embed._dgeom import distance_matrix_from_coordinates
from automol.embed._findif import central_difference

X = numpy.newaxis


def volume(xmat, idxs):
<<<<<<< HEAD
    """ calculate tetrahedral volume
=======
    """ calculate signed tetrahedral volume for a tetrad of atoms

    for a tetrad of four atoms (1, 2, 3, 4) around a central atom, the signed
    volume formula of this tetrahedral pyramid is given by

        d12 . (d13 x d14)

    where dij = rj - ri, . is the dot product, and x is the cross product
>>>>>>> Cleans up example
    """
    idxs = list(idxs)
    xyzs = xmat[:, :3][idxs]
    d12 = xyzs[1] - xyzs[0]
    d13 = xyzs[2] - xyzs[0]
    d14 = xyzs[3] - xyzs[0]
    vol = numpy.dot(d12, numpy.cross(d13, d14))
    return vol


def volume_gradient(xmat, idxs):
    """ calculate the tetrahedral volume gradient for a tetrad of atoms
    """
    idxs = list(idxs)
    xyzs = xmat[:, :3][idxs]

    grad = numpy.zeros_like(xmat)
    grad[idxs[0], :3] = (numpy.cross(xyzs[1], xyzs[3]-xyzs[2]) -
                         numpy.cross(xyzs[2], xyzs[3]))
    grad[idxs[1], :3] = +numpy.cross(xyzs[2]-xyzs[0], xyzs[3]-xyzs[0])
    grad[idxs[2], :3] = -numpy.cross(xyzs[1]-xyzs[0], xyzs[3]-xyzs[0])
    grad[idxs[3], :3] = +numpy.cross(xyzs[1]-xyzs[0], xyzs[2]-xyzs[0])

    return grad


def error_function_(lmat, umat, chip_dct=None, wdist=1., wchip=1., wdim4=1.,
                    leps=0.1, ueps=0.1):
    """ the embedding error function

    :param lmat: lower-bound distance matrix
    :param umat: upper-bound distance matrix
    :param chip_dct: chirality and planarity constraints; the keys are tuples
        of four atoms, the values are lower and upper bounds on the four-point
        signed volume of these atoms
    :param wdist: weight on the distance constraint
    :param wchip: weight on the chirality/planarity constraint
    :param wdim4: weight on the fourth dimension constraint
    :param leps: denominator epsilon for lower bound distances
    :param ueps: denominator epsilon for upper bound distances
    """
    triu = numpy.triu_indices_from(lmat)
    chip_dct = {} if chip_dct is None else chip_dct

    def _function(xmat):
        dmat = distance_matrix_from_coordinates(xmat)

        # distance error (equation 61 in the paper referenced above)
        ltf = ((lmat**2-dmat**2) / (leps**2+dmat**2))[triu]
        utf = ((dmat**2-umat**2) / (ueps**2+umat**2))[triu]
        ltf *= (ltf > 0.)
        utf *= (utf > 0.)
        dist_err = wdist * (numpy.vdot(utf, utf) + numpy.vdot(ltf, ltf))

        # chirality/planarity error (equation 62 in the paper referenced above)
<<<<<<< HEAD
        chip_err = wchip * sum(max(0, volume(xmat, idxs) - uvol)**2 +
                               max(0, lvol - volume(xmat, idxs))**2
                               for idxs, (lvol, uvol) in chip_dct.items())

        vols = numpy.array([volume(xmat, idxs) for idxs in chip_dct.keys()])
        lvols, uvols = map(numpy.array, zip(*chip_dct.values()))
        ltv = (lvols - vols) * (vols < lvols)
        utv = (vols - uvols) * (vols > uvols)
        alt_chip_err = wchip * (numpy.vdot(ltv, ltv) + numpy.vdot(utv, utv))
        print(alt_chip_err)
        print('chip err diff', chip_err - alt_chip_err)

        # fourth-dimension error
        dim4_err = wdim4 * numpy.vdot(xmat[:, 3], xmat[:, 3])
=======
        if chip_dct:
            vols = numpy.array(
                [volume(xmat, idxs) for idxs in chip_dct.keys()])
            lvols, uvols = map(numpy.array, zip(*chip_dct.values()))
            ltv = (lvols - vols) * (vols < lvols)
            utv = (vols - uvols) * (vols > uvols)
            chip_err = wchip * (numpy.vdot(ltv, ltv) + numpy.vdot(utv, utv))
        else:
            chip_err = 0.

        # fourth-dimension error
        if numpy.shape(xmat)[1] == 4:
            dim4_err = wdim4 * numpy.vdot(xmat[:, 3], xmat[:, 3])
        else:
            dim4_err = 0.
>>>>>>> Cleans up example

        return dist_err + chip_err + dim4_err

    return _function


<<<<<<< HEAD
def error_function_numerical_gradient_(lmat, umat, chip_dct=None,
                                       wdist=1., wchip=1., wdim4=1.,
                                       leps=0.1, ueps=0.1):
    """ the gradient of the distance error function

    (For testing purposes only; Used to check the analytic gradient formula.)
    """

    erf_ = error_function_(lmat, umat, chip_dct,
                           wdist=wdist, wchip=wchip, wdim4=wdim4,
                           leps=leps, ueps=ueps)

    def _gradient(xmat):
        grad = central_difference(erf_, xmat, npts=11)
        return grad

    return _gradient


=======
>>>>>>> Cleans up example
def error_function_gradient_(lmat, umat, chip_dct=None,
                             wdist=1., wchip=1., wdim4=1., leps=0.1, ueps=0.1):
    """ the embedding error function gradient

    :param lmat: lower-bound distance matrix
    :param umat: upper-bound distance matrix
    :param chip_dct: chirality and planarity constraints; the keys are tuples
        of four atoms, the values are lower and upper bounds on the four-point
        signed volume of these atoms
    :param wdist: weight on the distance constraint
    :param wchip: weight on the chirality/planarity constraint
    :param wdim4: weight on the fourth dimension constraint
    :param leps: denominator epsilon for lower bound distances
    :param ueps: denominator epsilon for upper bound distances
    """
    # natms = len(lmat)

    def _gradient(xmat):
        dmat = distance_matrix_from_coordinates(xmat)

        # distance error gradient
        utf = (dmat**2-umat**2) / (ueps**2+umat**2)
        ltf = (lmat**2-dmat**2) / (leps**2+dmat**2)
        utg = (+4.*utf/(ueps**2+umat**2))*(utf > 0.)
        ltg = (-4.*ltf/(leps**2+dmat**2)**2)*(leps**2+lmat**2)*(ltf > 0.)
        utg = utg[:, :, X]
        ltg = ltg[:, :, X]
        xmx = xmat[:, X, :] - xmat[X, :, :]
        dist_grad = numpy.sum(xmx*(ltg + utg), axis=1)
        dist_grad *= wdist

        # chirality/planarity error gradient
<<<<<<< HEAD
        vols = numpy.array([volume(xmat, idxs) for idxs in chip_dct.keys()])
        vol_grads = numpy.array(
            [volume_gradient(xmat, idxs) for idxs in chip_dct.keys()])
        lvols, uvols = map(numpy.array, zip(*chip_dct.values()))
        ltv = (lvols - vols) * (vols < lvols)
        utv = (vols - uvols) * (vols > uvols)
        ltg = -2. * ltv * vol_grads
        utg = +2. * utv * vol_grads
        chip_grad = numpy.sum(ltg+utg, axis=0)
        chip_grad *= wchip

        # fourth-dimension error gradient
        dim3_zeros = numpy.zeros_like(xmat[:, :3])
        dim4_grad = numpy.hstack([dim3_zeros, 2*xmat[:, 3:]])
        dim4_grad *= wdim4
=======
        if chip_dct:
            vols = numpy.array(
                [volume(xmat, idxs) for idxs in chip_dct.keys()])
            vol_grads = numpy.array(
                [volume_gradient(xmat, idxs) for idxs in chip_dct.keys()])
            lvols, uvols = map(numpy.array, zip(*chip_dct.values()))
            ltv = (lvols - vols) * (vols < lvols)
            utv = (vols - uvols) * (vols > uvols)
            ltg = -2. * ltv[:, X, X] * vol_grads
            utg = +2. * utv[:, X, X] * vol_grads
            chip_grad = numpy.sum(ltg+utg, axis=0)
            chip_grad *= wchip
        else:
            chip_grad = numpy.zeros_like(xmat)

        # fourth-dimension error gradient
        if numpy.shape(xmat)[1] == 4:
            dim3_zeros = numpy.zeros_like(xmat[:, :3])
            dim4_grad = numpy.hstack([dim3_zeros, 2*xmat[:, 3:]])
            dim4_grad *= wdim4
        else:
            dim4_grad = numpy.zeros_like(xmat)
>>>>>>> Cleans up example

        return dist_grad + chip_grad + dim4_grad

    return _gradient


<<<<<<< HEAD
=======
def error_function_numerical_gradient_(lmat, umat, chip_dct=None,
                                       wdist=1., wchip=1., wdim4=1.,
                                       leps=0.1, ueps=0.1):
    """ the gradient of the distance error function

    (For testing purposes only; Used to check the analytic gradient formula.)
    """

    erf_ = error_function_(lmat, umat, chip_dct,
                           wdist=wdist, wchip=wchip, wdim4=wdim4,
                           leps=leps, ueps=ueps)

    def _gradient(xmat):
        grad = central_difference(erf_, xmat, npts=11)
        return grad

    return _gradient


>>>>>>> Cleans up example
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


<<<<<<< HEAD
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
=======
def cleaned_up_coordinates(xmat, lmat, umat, chip_dct=None,
                           thresh=5e-2, maxiter=None,
                           # pre_thresh=5e-1, pre_maxiter=None,
                           chi_flip=False):
    """ clean up coordinates by conjugate-gradients error minimization

    :param xmat: the initial guess coordinates to be cleaned up
    :param lmat: lower-bound distance matrix
    :param umat: upper-bound distance matrix
    :param chip_dct: chirality and planarity constraints; the keys are tuples
        of four atoms, the values are lower and upper bounds on the four-point
        signed volume of these atoms
    :param thresh: convergence threshold, specifying the maximum gradient value
    :param maxiter: maximum number of iterations; default is three times the
        number of coordinates
    :param pre_thresh: convergence threshold for 4-dimensional pre-optimization
    :param pre_maxiter: maximum number of iterations for 4-dimensional
        pre-optimization; default is half of `maxiter`
    :chi_flip: whether or not to invert the structure if more than half of the
        chiralities are reversed
    """
    if chi_flip:
        raise NotImplementedError("Chirality flip not yet implemented!")

    # thresh1 = pre_thresh
    thresh2 = thresh

    maxiter2 = int(numpy.size(xmat) * 2 if maxiter is None else maxiter)
    # maxiter1 = int(3 if pre_maxiter is None else pre_maxiter)

    # fun1_ = error_function_(lmat, umat, chip_dct, wdim4=0.)
    fun2_ = error_function_(lmat, umat, chip_dct, wdim4=1.)
    # grad1_ = error_function_gradient_(lmat, umat, chip_dct, wdim4=0.)
    grad2_ = error_function_gradient_(lmat, umat, chip_dct, wdim4=1.)

    # xmat, _ = minimize_error(xmat, fun1_, grad1_, thresh1, maxiter1)
    xmat, conv = minimize_error(xmat, fun2_, grad2_, thresh2, maxiter2)
    return xmat, conv


def minimize_error(xmat, fun_, grad_, thresh=1e-1, maxiter=None):
    """ do conjugate-gradients error minimization

    :param fun_: a callable error function
    :param grad_: a callable error gradient function
    :param thresh: convergence threshold, specifying the maximum gradient value

    :returns: the optimized coordinates and a boolean which is True if
        converged and False if not
    """
    maxiter = numpy.size(xmat) * 3 if maxiter is None else maxiter

    sd0 = None
    cd0 = None
    print('Initial error:', fun_(xmat))
>>>>>>> Cleans up example

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

<<<<<<< HEAD
        print('function:', fun_(xmat))
        print('gradient:', grad_(xmat))

    print('niter:', niter)
=======
        print('Iteration {:d}'.format(niter))
        print('\tError:', fun_(xmat))
        print('\tMax gradient:', grad_max)
        print()

    print('Niter:', niter)
    print('Converged:', converged)
    print()
>>>>>>> Cleans up example

    return xmat, converged


if __name__ == '__main__':
<<<<<<< HEAD
    # numpy.random.seed(7)
    import automol

    # ICH = 'InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1'
    ICH = 'InChI=1S/CHBrClF/c2-1(3)4/h1H'
    GEO = automol.inchi.geometry(ICH)
    ICH = automol.geom.inchi(GEO)

    GRA = automol.geom.graph(GEO)
    print(automol.graph.string(GRA))

    STE_KEYS = automol.graph.atom_neighbor_keys(GRA)[0]
    CHIP_IDXS = automol.graph.stereo_sorted_atom_neighbor_keys(
        GRA, 0, STE_KEYS)

    CHIP_DCT = {CHIP_IDXS: (-999., -7.)}
    # CHIP_DCT = {CHIP_IDXS: (7., 999.)}

    KEYS = sorted(automol.graph.atom_keys(GRA))
    LMAT, UMAT = automol.graph.embed.distance_bounds_matrices(GRA, KEYS)
    XMAT = automol.graph.embed.sample_raw_distance_coordinates(GRA, KEYS,
                                                               dim4=True)

    print(CHIP_IDXS)
    print('volume', volume(XMAT, CHIP_IDXS))
    FUN_ = error_function_(LMAT, UMAT, CHIP_DCT, wdist=0., wchip=1., wdim4=0.)
    print(FUN_(XMAT))

    GRAD1_ = error_function_numerical_gradient_(
        LMAT, UMAT, CHIP_DCT, wdist=0., wchip=1., wdim4=0.)
    GRAD1 = GRAD1_(XMAT)

    GRAD2_ = error_function_gradient_(
        LMAT, UMAT, CHIP_DCT, wdist=0., wchip=1., wdim4=0.)
    GRAD2 = GRAD2_(XMAT)

    print(numpy.round(GRAD1, 2))
    print(numpy.round(GRAD2, 2))
    print(numpy.amax(numpy.abs(GRAD1-GRAD2)))

    # SYMS = list(map(automol.graph.atom_symbols(GRA).__getitem__, KEYS))
    # XYZS = XMAT[:, :3]
    # GEO = automol.create.geom.from_data(SYMS, XYZS, angstrom=True)
    # print(automol.geom.string(GEO))

    # OLD

    # # # morphine without stereo
    # # ICH = automol.smiles.inchi('CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O')

    # # morphine with stereo
    # ICH = ('InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)'
    #        '4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/'
    #        't10-,11+,13-,16-,17-/m0/s1')

    # # 1. Generate distance bounds matrices, L and U
    # # GRA = automol.inchi.graph(ICH)
    # GRA = automol.graph.from_string("""
    #     atoms:
    #       1: {symbol: C, implicit_hydrogen_valence: 3, stereo_parity: null}
    #       2: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: null}
    #       3: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: null}
    #       4: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: null}
    #       5: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: null}
    #       6: {symbol: C, implicit_hydrogen_valence: 2, stereo_parity: null}
    #       7: {symbol: C, implicit_hydrogen_valence: 2, stereo_parity: null}
    #       8: {symbol: C, implicit_hydrogen_valence: 2, stereo_parity: null}
    #       9: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
    #       10: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: true}
    #       11: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: true}
    #       12: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
    #       13: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: true}
    #       14: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
    #       15: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
    #       16: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: false}
    #       17: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: false}
    #       18: {symbol: N, implicit_hydrogen_valence: 0, stereo_parity: null}
    #       19: {symbol: O, implicit_hydrogen_valence: 1, stereo_parity: null}
    #       20: {symbol: O, implicit_hydrogen_valence: 1, stereo_parity: null}
    #       21: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
    #     bonds:
    #       1-18: {order: 1, stereo_parity: null}
    #       2-4: {order: 1, stereo_parity: null}
    #       2-9: {order: 1, stereo_parity: null}
    #       3-5: {order: 1, stereo_parity: null}
    #       3-10: {order: 1, stereo_parity: null}
    #       4-12: {order: 1, stereo_parity: null}
    #       5-13: {order: 1, stereo_parity: null}
    #       6-7: {order: 1, stereo_parity: null}
    #       6-17: {order: 1, stereo_parity: null}
    #       7-18: {order: 1, stereo_parity: null}
    #       8-9: {order: 1, stereo_parity: null}
    #       8-11: {order: 1, stereo_parity: null}
    #       9-14: {order: 1, stereo_parity: null}
    #       10-11: {order: 1, stereo_parity: null}
    #       10-17: {order: 1, stereo_parity: null}
    #       11-18: {order: 1, stereo_parity: null}
    #       12-15: {order: 1, stereo_parity: null}
    #       12-19: {order: 1, stereo_parity: null}
    #       13-16: {order: 1, stereo_parity: null}
    #       13-20: {order: 1, stereo_parity: null}
    #       14-15: {order: 1, stereo_parity: null}
    #       14-17: {order: 1, stereo_parity: null}
    #       15-21: {order: 1, stereo_parity: null}
    #       16-17: {order: 1, stereo_parity: null}
    #       16-21: {order: 1, stereo_parity: null}
    # """)
    # GRA = automol.graph.explicit(GRA)
    # KEYS = sorted(automol.graph.atom_keys(GRA))
    # LMAT, UMAT = automol.graph.embed.distance_bounds_matrices(GRA, KEYS)
    # XMAT = automol.graph.embed.sample_raw_distance_coordinates(GRA, KEYS,
    #                                                            dim4=True)

    # NGB_KEYS_DCT = automol.graph.atom_neighbor_keys(GRA)
    # STE_KEYS = automol.graph.atom_stereo_keys(GRA)
    # STE_NGB_KEYS_LST = [
    #     automol.graph.stereo_sorted_atom_neighbor_keys(
    #         GRA, KEY, NGB_KEYS_DCT[KEY])
    #     for KEY in STE_KEYS]
    # print(STE_KEYS)
    # print(STE_NGB_KEYS_LST)
    # STE_IDXS = list(map(KEYS.index, STE_KEYS))
    # STE_NGB_IDXS_LST = [tuple(map(KEYS.index, NGB_KEYS))
    #                     for NGB_KEYS in STE_NGB_KEYS_LST]
    # print(STE_IDXS)
    # print(STE_NGB_IDXS_LST)
    # STE_NGB_IDXS_DCT = dict(zip(STE_IDXS, STE_NGB_IDXS_LST))
    # print(STE_NGB_IDXS_DCT)

    # # XMAT, CONVERGED = cleaned_up_coordinates(XMAT, LMAT, UMAT)
    # # print(CONVERGED)

    # # SYMS = list(map(automol.graph.atom_symbols(GRA).__getitem__, KEYS))
    # # XYZS = XMAT[:, :3]
    # # GEO = automol.create.geom.from_data(SYMS, XYZS, angstrom=True)
    # # print(automol.geom.string(GEO))

    # # GRA2 = automol.geom.connectivity_graph(GEO)
    # # print(GRA == GRA2)
=======
    numpy.random.seed(1)
    import automol

    ICH = automol.smiles.inchi('O')  # water
    # ICH = automol.smiles.inchi('CO')  # methanol
    # ICH = automol.smiles.inchi('C1CCCCC1')  # hexane
    # ICH = automol.smiles.inchi('C1C2CC3CC1CC(C2)C3')  # adamantane
    GEO = automol.inchi.geometry(ICH)
    GRA = automol.geom.graph(GEO)

    KEYS = sorted(automol.graph.atom_keys(GRA))
    LMAT, UMAT = automol.graph.embed.distance_bounds_matrices(GRA, KEYS)
    P_DCT = automol.graph.embed.planarity_constraint_bounds(GRA, KEYS)
    XMAT = sample_raw_distance_coordinates(GRA, KEYS, dim4=True)

    print(numpy.shape(XMAT))
    # XMAT, CONV = cleaned_up_coordinates(XMAT, LMAT, UMAT, chip_dct=P_DCT)
    XMAT, CONV = cleaned_up_coordinates(XMAT, LMAT, UMAT, chip_dct={})
    print(CONV)

    SYMS = list(map(automol.graph.atom_symbols(GRA).__getitem__, KEYS))
    XYZS = XMAT[:, :3]
    GEO = automol.create.geom.from_data(SYMS, XYZS, angstrom=True)
    print(automol.geom.string(GEO))

    GRA = automol.graph.without_stereo_parities(GRA)
    GRA2 = automol.geom.connectivity_graph(GEO)
    print(GRA == GRA2)
    print(ICH)
>>>>>>> Cleans up example
