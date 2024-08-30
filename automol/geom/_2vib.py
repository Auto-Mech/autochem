"""Vibrational analysis."""

from collections.abc import Sequence

import numpy
from qcelemental import constants as qcc

from .base import masses, rotational_normal_modes, translational_normal_modes

MatrixLike = Sequence[Sequence[float]] | numpy.ndarray


def normal_mode_projection(
    geo, trans: bool = False, rot: bool = False
) -> numpy.ndarray:
    """Get the matrix for projecting onto a subset of normal modes.

    Note: This must be applied to the *mass-weighted* normal modes!

    Uses SVD to calculate an orthonormal basis for the null space of the external normal
    modes, which is the space of internal normal modes.

    :param geo: The geometry
    :param trans: Keep translations, instead of removing them?
    :param rot: Keep rotations, instead of removing them?
    :return: The projection onto a subset normal modes
    """
    coos_lst = []
    if not trans:
        coos_lst.append(translational_normal_modes(geo, mass_weight=True))
    if not rot:
        coos_lst.append(rotational_normal_modes(geo, mass_weight=True))
    coos = numpy.hstack(coos_lst)
    dim = numpy.shape(coos)[-1]
    coo_basis, *_ = numpy.linalg.svd(coos, full_matrices=True)
    proj = coo_basis[:, dim:]
    return proj


def vibrational_analysis(
    geo, hess: MatrixLike, trans: bool = False, rot: bool = False, wavenum: bool = True
) -> tuple[tuple[float, ...], numpy.ndarray]:
    """Get vibrational modes and frequencies from the Hessian matrix.

    :param geo: The molecular geometry
    :param hess: The Hessian matrix
    :param trans: Keep translations, instead of removing them?
    :param rot: Keep rotations, instead of removing them?
    :param wavenum: Return frequencies (cm^-1), instead of force constants (a.u.)?
    :return: The vibrational frequencies (or force constants) and normal modes
    """
    # 1. Mass-weight the Hessian matrix
    mw_vec = numpy.sqrt(numpy.repeat(masses(geo), 3))
    hess_mw = hess / mw_vec[:, numpy.newaxis] / mw_vec[numpy.newaxis, :]

    # 2. Project onto the space of internal motions
    #       K = Qmwt Hmw Qmw = (PI)t Hmw (PI) = It Hint I
    proj = normal_mode_projection(geo, trans=trans, rot=rot)
    hess_proj = proj.T @ hess_mw @ proj

    # 2. Compute eigenvalues and eigenvectors of the mass-weighted Hessian matrix
    eig_vals, eig_vecs = numpy.linalg.eigh(hess_proj)

    # 3. Un-mass-weight the normal coordinates
    #       Qmw = PI (see above)    Q = Qmw / mw_vec
    norm_coos_mw = numpy.dot(proj, eig_vecs) / mw_vec[:, numpy.newaxis]
    norm_coos = norm_coos_mw / mw_vec[:, numpy.newaxis]
    if not wavenum:
        return eig_vals, norm_coos

    # 4. Get wavenumbers from a.u. force constants
    har2J = qcc.conversion_factor("hartree", "J")
    amu2kg = qcc.conversion_factor("atomic_mass_unit", "kg")
    bohr2m = qcc.conversion_factor("bohr", "meter")
    sol = qcc.get("speed of light in vacuum") * 100  # in cm / s
    to_inv_cm = numpy.sqrt(har2J / (amu2kg * bohr2m * bohr2m)) / (sol * 2 * numpy.pi)
    freqs = numpy.sqrt(numpy.complex_(eig_vals)) * to_inv_cm
    freqs = tuple(map(float, numpy.real(freqs) - numpy.imag(freqs)))

    return freqs, norm_coos
