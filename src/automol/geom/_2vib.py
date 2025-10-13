"""Vibrational analysis."""

from collections.abc import Callable, Collection, Sequence

import numpy
from qcelemental import constants as qcc

from .. import util
from ..graph import base as graph_base
from ._1conv import graph
from .base import (
    coordinates,
    count,
    mass_weight_vector,
    masses,
    rotational_normal_modes,
    translational_normal_modes,
)

MatrixLike = Sequence[Sequence[float]] | numpy.ndarray


def vibrational_analysis(
    geo,
    hess: MatrixLike,
    trans: bool = False,
    rot: bool = False,
    tors: bool = True,
    gra: object | None = None,
    bkeys: Sequence[Collection[int]] | None = None,
    with_h_rotors: bool = True,
    with_ch_rotors: bool = True,
    wavenum: bool = True,
) -> tuple[tuple[float, ...], numpy.ndarray]:
    """Get vibrational modes and frequencies from the Hessian matrix.

    :param geo: The molecular geometry
    :param hess: The Hessian matrix
    :param trans: Keep translations, instead of removing them?
    :param rot: Keep rotations, instead of removing them?
    :param tors: Keep torsions, instead of removing them?
    :param gra: For torsions, a graph specifying connectivity
    :param bkeys: For torsions, specify the rotational bonds by index
    :param with_h_rotors: For torsions, include XH rotors?
    :param with_ch_rotors: For torsions, include CH rotors?
    :param wavenum: Return frequencies (cm^-1), instead of force constants (a.u.)?
    :return: The vibrational frequencies (or force constants) and normal modes
    """
    # 1. Mass-weight the Hessian matrix
    mw_vec = numpy.sqrt(numpy.repeat(masses(geo), 3))
    hess_mw = hess / mw_vec[:, numpy.newaxis] / mw_vec[numpy.newaxis, :]

    # 2. Project onto the space of internal motions
    #       K = Qmwt Hmw Qmw = (PI)t Hmw (PI) = It Hint I
    proj = normal_mode_projection(
        geo,
        trans=trans,
        rot=rot,
        tors=tors,
        gra=gra,
        bkeys=bkeys,
        with_h_rotors=with_h_rotors,
        with_ch_rotors=with_ch_rotors,
    )
    hess_proj = proj.T @ hess_mw @ proj

    # 2. Compute eigenvalues and eigenvectors of the mass-weighted Hessian matrix
    eig_vals, eig_vecs = numpy.linalg.eigh(hess_proj)

    # 3. Un-mass-weight the normal coordinates
    #       Qmw = PI (see above)    Q = Qmw / mw_vec
    norm_coos_mw = numpy.dot(proj, eig_vecs) / mw_vec[:, numpy.newaxis]
    norm_coos = norm_coos_mw / mw_vec[:, numpy.newaxis]
    norm_coos = norm_coos / numpy.linalg.norm(norm_coos, axis=0)
    if not wavenum:
        return eig_vals, norm_coos

    # 4. Get wavenumbers from a.u. force constants
    har2J = qcc.conversion_factor("hartree", "J")
    amu2kg = qcc.conversion_factor("atomic_mass_unit", "kg")
    bohr2m = qcc.conversion_factor("bohr", "meter")
    sol = qcc.get("speed of light in vacuum") * 100  # in cm / s
    to_inv_cm = numpy.sqrt(har2J / (amu2kg * bohr2m * bohr2m)) / (sol * 2 * numpy.pi)
    freqs = numpy.sqrt(numpy.complex128(eig_vals)) * to_inv_cm
    freqs = tuple(map(float, numpy.real(freqs) - numpy.imag(freqs)))

    return freqs, norm_coos


def normal_mode_projection(
    geo,
    trans: bool = False,
    rot: bool = False,
    tors: bool = True,
    gra: object | None = None,
    bkeys: Sequence[Collection[int]] | None = None,
    with_h_rotors: bool = True,
    with_ch_rotors: bool = True,
) -> numpy.ndarray:
    """Get the matrix for projecting onto a subset of normal modes.

    Note: This must be applied to the *mass-weighted* normal modes!

    Uses SVD to calculate an orthonormal basis for the null space of the external normal
    modes, which is the space of internal normal modes.

    :param geo: The geometry
    :param trans: Keep translations, instead of removing them?
    :param rot: Keep rotations, instead of removing them?
    :param tors: Keep torsions, instead of removing them?
    :param gra: For torsions, a graph specifying connectivity
    :param bkeys: For torsions, specify the rotational bonds by index
    :param with_h_rotors: For torsions, include XH rotors?
    :param with_ch_rotors: For torsions, include CH rotors?
    :return: The projection onto a subset normal modes
    """
    coos_lst = []
    if not trans:
        coos_lst.append(translational_normal_modes(geo, mass_weight=True))
    if not rot:
        coos_lst.append(rotational_normal_modes(geo, mass_weight=True))
    if not tors:
        coos_lst.append(
            torsional_normal_modes(
                geo,
                gra=gra,
                bkeys=bkeys,
                with_h_rotors=with_h_rotors,
                with_ch_rotors=with_ch_rotors,
            )
        )

    # If no modes are being project, return an identity matrix
    if not coos_lst:
        return numpy.eye(count(geo) * 3)

    coos = numpy.hstack(coos_lst)
    dim = numpy.shape(coos)[-1]
    coo_basis, *_ = numpy.linalg.svd(coos, full_matrices=True)
    proj = coo_basis[:, dim:]
    return proj


# Torsions
def torsional_normal_modes(
    geo,
    gra: object | None = None,
    bkeys: Sequence[Collection[int]] | None = None,
    mass_weight: bool = True,
    with_h_rotors: bool = True,
    with_ch_rotors: bool = True,
) -> numpy.ndarray:
    """Calculate the torsional normal modes at rotational bonds.

    :param geo: The geometry
    :param gra: A graph specifying the connectivity of the geometry
    :param bkeys: Specify the rotational bonds by index (ignores `with_x_rotors` flags);
        If `None`, they will be determined automatically
    :param mass_weight: Return mass-weighted normal modes?
    :param with_h_rotors: Include rotors neighbored by only hydrogens on one side?
    :param with_ch_rotors: Include rotors with C neighbored by only hydrogens?
    :return: The torsional normal modes, as a 3N x n_tors matrix
    """
    gra = graph(geo, stereo=False) if gra is None else gra
    bkeys = (
        graph_base.rotational_bond_keys(
            gra,
            with_h_rotors=with_h_rotors,
            with_ch_rotors=with_ch_rotors,
            with_rings_rotors=False,
        )
        if bkeys is None
        else bkeys
    )

    all_idxs = list(range(count(geo)))
    tors_coos = []
    for bkey in bkeys:
        axis_idxs = sorted(bkey)
        groups = graph_base.rotational_groups(gra, *axis_idxs)
        tors_ = torsional_motion_calculator_(geo, axis_idxs, groups)
        tors_coo = numpy.concatenate([tors_(idx) for idx in all_idxs])
        tors_coos.append(tors_coo)
    tors_coos = numpy.transpose(tors_coos)

    if mass_weight:
        tors_coos *= mass_weight_vector(geo)[:, numpy.newaxis]
    return tors_coos


def torsional_motion_calculator_(
    geo, axis_idxs: Sequence[int], groups: Sequence[Sequence[int]]
) -> Callable[[int], numpy.ndarray]:
    """Generate a torsional motion calculator.

    :param geo: The geometry
    :param axis_idxs: The indices of the rotational axis
    :param groups: The groups of the rotational axis
    """
    group1, group2 = groups
    xyz1, xyz2 = coordinates(geo, idxs=axis_idxs)
    axis1 = util.vector.unit_norm(numpy.subtract(xyz1, xyz2))
    axis2 = numpy.negative(axis1)

    def _torsional_motion(idx: int) -> numpy.ndarray:
        """Determine the torsional motion of an atom."""
        # If it is one of the axis atoms, it doesn't move
        if idx in axis_idxs:
            return numpy.zeros(3)
        # If it is in the first group, subtract and cross
        (xyz,) = coordinates(geo, idxs=(idx,))
        if idx in group1:
            dxyz = numpy.subtract(xyz, xyz1)
            return numpy.cross(dxyz, axis1)
        if idx in group2:
            dxyz = numpy.subtract(xyz, xyz2)
            return numpy.cross(dxyz, axis2)

        raise ValueError(f"{idx} not in {axis_idxs} {group1} {group2}")

    return _torsional_motion
