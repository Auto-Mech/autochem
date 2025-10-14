"""Interface to the MolSym code
"""
import qcelemental as qcel

import molsym
from molsym.symtext.symel import PointGroup

from .base import xyz_string


def point_group_from_geometry(geo, tol: float = 0.05) -> PointGroup:
    """Generate a MolSym symmetry object from a geometry

    :param geo: A geometry
    :type geo: automol geom data structure
    :param tol: Tolerance threshold (in Bohr) for symmetry detection
    :type tol: float, optional
    :return: A PointGroup object
    :rtype: PointGroup
    """
    qc_mol = qcel.models.Molecule.from_data(xyz_string(geo))
    ms_mol = molsym.Molecule.from_schema(qc_mol.dict())
    ms_mol.tol = tol
    ms_mol.translate(ms_mol.find_com())
    pg_str, *_ = molsym.find_point_group(ms_mol)
    pg_obj = PointGroup.from_string(pg_str)
    return pg_obj


def point_group_is_chiral(pg_obj: PointGroup) -> bool:
    """Based on its point group, is this molecule chiral?

    :param pg_obj: A MolSym point group symmetry object
    :type pg_obj: molsym.PointGroup
    :return: `True` if it is, `False` if it isn't
    :rtype: bool
    """
    return pg_obj.family in ("C", "D") and pg_obj.subfamily is None


def point_group_symmetry_number(pg_obj: PointGroup) -> int:
    """Get the external symmetry number of a molecule

    See Table II [here](https://cccbdb.nist.gov/thermo.asp) for a definition.

    :param pg_obj: A MolSym point group symmetry object
    :type pg_obj: molsym.PointGroup
    :return: The external symmetry number
    :rtype: int
    """
    sym_num = None

    if pg_obj.family == "C":
        if pg_obj.n == 0 or pg_obj.n is None:
            sym_num = 1
        else:
            sym_num = pg_obj.n
    elif pg_obj.family == "D":
        if pg_obj.n == 0:
            sym_num = 2
        else:
            sym_num = 2 * pg_obj.n
    elif pg_obj.family == "S":
        sym_num = pg_obj.n >> 1
    elif pg_obj.family == "T":
        sym_num = 12
    elif pg_obj.family == "O":
        sym_num = 24
    elif pg_obj.family == "I":
        sym_num = 60

    return sym_num
