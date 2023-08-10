"""Structural heuristics
"""
from phydat import ptab


def bond_distance(symb1: str, symb2: str, angstrom: bool = True) -> float:
    """The heuristic bond distance between two atoms, based on their symbols

    Returns whichever is smaller of (a.) the sum of covalent radii, and (b.) the average
    vdw radius.

    :param symb1: The first atom symbol
    :type symb1: int
    :param symb2: The second atom symbol
    :type symb2: int
    :param angstrom: Return in angstroms intead of bohr?, defaults to True
    :type angstrom: bool, optional
    :return: The heuristic distance
    :rtype: float
    """
    rcov1 = ptab.covalent_radius(symb1, angstrom=angstrom)
    rcov2 = ptab.covalent_radius(symb2, angstrom=angstrom)
    rvdw1 = ptab.van_der_waals_radius(symb1, angstrom=angstrom)
    rvdw2 = ptab.van_der_waals_radius(symb2, angstrom=angstrom)

    cov_dist = rcov1 + rcov2
    vdw_dist = (rvdw1 + rvdw2) / 2.0

    return min(cov_dist, vdw_dist)


def bond_distance_limit(
    symb1: str, symb2: str, bdist_factor: float = None, angstrom: bool = True
) -> float:
    """The heuristic bond distance limit (largest possible bond distance) between two
    atoms, based on their symbols

    Returns `bdist_factor` times whichever is larger of of (a.) the sum of covalent
    radii, and (b.) the average vdw radius.

    :param symb1: The first atom symbol
    :type symb1: int
    :param symb2: The second atom symbol
    :type symb2: int
    :param bdist_factor: The multiplier on the distance limit, defaults to 1.05
    :type bdist_factor: float, optional
    :param angstrom: Return in angstroms intead of bohr?, defaults to True
    :type angstrom: bool, optional
    :return: The heuristic distance
    :rtype: float
    """
    bdist_factor = 1.05 if bdist_factor is None else bdist_factor

    rcov1 = ptab.covalent_radius(symb1, angstrom=angstrom)
    rcov2 = ptab.covalent_radius(symb2, angstrom=angstrom)
    rvdw1 = ptab.van_der_waals_radius(symb1, angstrom=angstrom)
    rvdw2 = ptab.van_der_waals_radius(symb2, angstrom=angstrom)

    cov_dist = rcov1 + rcov2
    vdw_dist = (rvdw1 + rvdw2) / 2.0

    return max(cov_dist, vdw_dist) * bdist_factor
