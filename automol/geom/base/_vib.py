"""Vibrational analysis."""

from collections.abc import Sequence

import numpy
from qcelemental import constants as qcc

from ._0core import masses, projection_onto_internal_modes

MatrixLike = Sequence[Sequence[float]] | numpy.ndarray


def vibrational_analysis(
    geo, hess: MatrixLike, freq: bool = True
) -> tuple[tuple[float, ...], numpy.ndarray]:
    """Get vibrational modes and frequencies from the Hessian matrix.

    :param geo: The molecular geometry
    :param hess: The Hessian matrix
    :param freq: Return frequencies (cm^-1), instead of force constants (a.u.)?
    :return: The vibrational frequencies (or force constants) and normal modes
    """
    # 1. Mass-weight the Hessian matrix
    mw_vec = numpy.sqrt(numpy.repeat(masses(geo), 3))
    hess_mw = hess / mw_vec[:, numpy.newaxis] / mw_vec[numpy.newaxis, :]

    # 2. Project onto the space of internal motions
    #       K = Qmwt Hmw Qmw = (PI)t Hmw (PI) = It Hint I
    int_proj = projection_onto_internal_modes(geo)
    hess_int = int_proj.T @ hess_mw @ int_proj

    # 2. compute eigenvalues and eigenvectors of the mass-weighted Hessian matrix
    eig_vals, eig_vecs = numpy.linalg.eigh(hess_int)

    # 3. un-mass-weight the normal coordinates
    #       Qmw = PI (see above)    Q = Qmw / mw_vec
    norm_coos_mw = numpy.dot(int_proj, eig_vecs) / mw_vec[:, numpy.newaxis]
    norm_coos = norm_coos_mw / mw_vec[:, numpy.newaxis]
    if not freq:
        return eig_vals, norm_coos

    # 4. get wavenumbers from a.u. force constants
    har2J = qcc.conversion_factor("hartree", "J")
    amu2kg = qcc.conversion_factor("atomic_mass_unit", "kg")
    bohr2m = qcc.conversion_factor("bohr", "meter")
    sol = qcc.get("speed of light in vacuum") * 100  # in cm / s
    to_inv_cm = numpy.sqrt(har2J / (amu2kg * bohr2m * bohr2m)) / (sol * 2 * numpy.pi)
    freqs = numpy.sqrt(numpy.complex_(eig_vals)) * to_inv_cm
    freqs = tuple(map(float, numpy.real(freqs) - numpy.imag(freqs)))

    return freqs, norm_coos
