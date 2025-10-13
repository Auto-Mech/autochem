"""Functions for describing intermolecular interactions in a geometry
"""
import itertools

import numpy

from ...util import dict_
from ._0core import count, distance, symbols


def has_low_relative_repulsion_energy(
    geo, ref_geo, model: str = "exp6", thresh=40.0
) -> bool:
    """Identify whether a geometry has low repulsion energy relative to a reference
    geometry

    :param geo: A geometry
    :type geo: automol geom data structure
    :param ref_geo: A reference geometry to compare to
    :type ref_geo: automol geom data structure
    :param model: The model potential to use, "exp6" (default) or "lj_12_6"
    :type model: str, optional
    :param thresh: Threshold for excess repulsion energy (kcal/mol), defaults to 40.0
    :type thresh: float, optional
    :return: `True` if it does, `False` if it doesn't
    :rtype: bool
    """

    pot = total_repulsion_energy(geo, model=model)
    ref_pot = total_repulsion_energy(ref_geo, model=model)
    return bool((pot - ref_pot) <= thresh)


def total_repulsion_energy(geo, model: str = "exp6") -> float:
    """Calculate the geometry's total repulsion energy using a model potential

    :param geo: A geometry
    :type geo: automol geom data structure
    :param model: The model potential to use, "exp6" (default) or "lj_12_6"
    :type model: str, optional
    :return: The repulsion energy
    :rtype: float
    """
    idxs = range(count(geo))
    return sum(
        repulsion_energy(geo, i1, i2, model=model)
        for i1, i2 in itertools.product(idxs, repeat=2)
        if i1 != i2
    )


def repulsion_energy(geo, idx1: int, idx2: int, model: str = "exp6") -> float:
    """Calculate the repulsion energy between two atoms using a model potential

    :param geo: A geometry
    :type geo: automol geom data structure
    :param idx1: The index of the first atom
    :type idx1: int
    :param idx2: The index of the second atom
    :type idx2: int
    :param model: The model potential to use, "exp6" (default) or "lj_12_6"
    :type model: str, optional
    :return: The repulsion energy
    :rtype: float
    """
    # If we are comparing the atom to itself, return None
    if idx1 == idx2:
        return None

    symbs = symbols(geo)
    symb1 = symbs[idx1]
    symb2 = symbs[idx2]
    dist = distance(geo, idx1, idx2, angstrom=True)
    return _model_potential_energy(symb1, symb2, dist, model=model)


# helpers
LJ_DCT = {
    # eps[whatever], sig[ang] params
    ("H", "H"): (0.25, 1.0),
    ("H", "C"): (0.25, 1.0),
    ("H", "O"): (0.25, 1.0),
    ("C", "C"): (0.25, 1.0),
    ("C", "O"): (0.25, 1.0),
    ("O", "O"): (0.25, 1.0),
}

EXP6_DCT = {
    # A, B, C params E[kcal] R[Ang]; R cutoff
    ("H", "H"): (2.442e3, 3.74, 48.8, 1.0),
    ("H", "C"): (6.45e3, 3.67, 116.0, 1.0),
    ("H", "N"): (6.45e3, 3.67, 116.0, 1.0),
    ("H", "O"): (6.45e3, 3.67, 116.0, 1.0),
    ("C", "C"): (7.69e4, 3.6, 460.0, 0.8),
    ("C", "O"): (7.69e4, 3.6, 460.0, 0.8),
    ("O", "O"): (7.69e4, 3.6, 460.0, 0.8),
    ("Cl", "Cl"): (0.0, 3.6, 460.0, 0.8),
    ("Cl", "C"): (0.0, 3.6, 460.0, 0.8),
    ("Cl", "O"): (0.0, 3.6, 460.0, 0.8),
    ("Cl", "H"): (0.0, 3.6, 460.0, 0.8),
    ("N", "N"): (7.69e4, 3.6, 460.0, 0.8),
    ("N", "O"): (7.69e4, 3.6, 460.0, 0.8),
    ("N", "C"): (7.69e4, 3.6, 460.0, 0.8),
}


def _model_potential_energy(
    symb1: str, symb2: str, dist: float, model: str = "exp6"
) -> float:
    """Calculate an interatomic potential value according to a certian model, given the
    atomic symbols and the distances

    :param symb1: The first atomic symbol
    :type symb1: str
    :param symb2: The second atomic symbol
    :type symb2: str
    :param dist: The distance, in angstroms
    :type dist: float
    :param model: The model potential to use, "exp6" (default) or "lj_12_6"
    :type model: str, optional
    :return: The repulsion energy
    :rtype: float
    """
    assert model in ("exp6", "lj_12_6"), f"model {model} != exp6 or lj_12_6"

    param_dct = EXP6_DCT if model == "exp6" else LJ_DCT
    pot_ = _exp6_potential_energy if model == "exp6" else _lj_potential_energy

    params = dict_.value_by_unordered_key(param_dct, (symb1, symb2))
    pot_val = pot_(dist, *params)
    return pot_val


def _lj_potential_energy(rdist, eps, sig):
    """Calculate potential energy value of two interacting bodies
    assuming a 12-6 Lennard-Jones potential.

    :param rdist: distance between two interacting bodies (Bohr)
    :type rdist: float
    :param eps: Lennard-Jones epsilon parameter for interaction (_)
    :type eps: float
    :param eps: Lennard-Jones sigma parameter for interaction (_)
    :type eps: float
    :rtpe: float
    """
    return (4.0 * eps) * ((sig / rdist) ** 12 - (sig / rdist) ** 6)


def _exp6_potential_energy(rdist, apar, bpar, cpar, rcut):
    """Calculate potential energy value of two interacting bodies
    assuming a modified Buckingham potential.

    :param rdist: distance between two interacting bodies (Bohr)
    :type rdist: float
    :param apar: potential parameter A
    :type apar: float
    :param bpar: potential parameter B
    :type bpar: float
    :param cpar: potential parameter C
    :type cpar: float
    :param rcut: threshhold where interaction potential becomes constant
    :type rcut: float
    :rtpe: float
    """
    if rdist < rcut:
        pot_val = apar * numpy.exp(-1.0 * bpar * rcut) - (cpar / rcut**6)
    else:
        pot_val = apar * numpy.exp(-1.0 * bpar * rdist) - (cpar / rdist**6)
    return pot_val
