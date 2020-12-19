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
from automol.embed._dgeom import distance_matrix_from_coordinates
from automol.embed._findif import central_difference

X = numpy.newaxis


def volume(xmat, idxs):
    """ calculate chiral volume
    """
    idxs = list(idxs)
    xyzs = xmat[:, :3][idxs]
    d12 = xyzs[1] - xyzs[0]
    d13 = xyzs[2] - xyzs[0]
    d14 = xyzs[3] - xyzs[0]
    vol = numpy.dot(d12, numpy.cross(d13, d14))
    return vol


def volume_gradient(xmat, idxs, idx):
    """ calculate the chiral volume gradient for a single atom
    """
    idxs = list(idxs)

    grad = numpy.zeros_like(xmat)

    if idx in idxs:
        xyzs = xmat[:, :3][idxs]
        pos = idxs.index(idx)

        if pos == 0:
            vec = (numpy.cross(xyzs[1], xyzs[3]-xyzs[2]) -
                   numpy.cross(xyzs[2], xyzs[3]))
        elif pos == 1:
            vec = +numpy.cross(xyzs[2]-xyzs[0], xyzs[3]-xyzs[0])
        elif pos == 2:
            vec = -numpy.cross(xyzs[1]-xyzs[0], xyzs[3]-xyzs[0])
        elif pos == 3:
            vec = +numpy.cross(xyzs[1]-xyzs[0], xyzs[2]-xyzs[0])

        grad[idx, :3] = vec

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
        utf = ((dmat**2-umat**2) / (ueps**2+umat**2))[triu]
        ltf = ((lmat**2-dmat**2) / (leps**2+dmat**2))[triu]
        utf *= (utf > 0.)
        ltf *= (ltf > 0.)
        dist_err = wdist * (numpy.vdot(utf, utf) + numpy.vdot(ltf, ltf))

        # chirality/planarity error (equation 62 in the paper referenced above)
        chip_err = wchip * sum(max(0, volume(xmat, idxs) - uvol)**2 +
                               max(0, lvol - volume(xmat, idxs))**2
                               for idxs, (lvol, uvol) in chip_dct.items())

        # fourth-dimension error
        dim4_err = wdim4 * numpy.vdot(xmat[:, 3], xmat[:, 3])

        return dist_err + chip_err + dim4_err

    return _function


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

    GRAD_ = error_function_numerical_gradient_(
        LMAT, UMAT, CHIP_DCT, wdist=0., wchip=1., wdim4=0.)
    print(GRAD_(XMAT))

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
