"""Interface to the MolSym code
"""
import qcelemental as qcel
from automol.geom.base import xyz_string
from molsym.molecule import Molecule
from molsym.symtext.symtext import Symtext


def symtext_from_geometry(geo, tol: float=0.05) -> Symtext:
    """Generate a MolSym symmetry object from a geometry

    :param geo: A geometry
    :type geo: automol geom data structure
    :param tol: Tolerance threshold (in Bohr) for symmetry detection
    :type tol: float, optional
    :return: A MolSym symmetry object
    :rtype: Symtext
    """
    qcm = qcel.models.Molecule.from_data(xyz_string(geo))
    msm = Molecule.from_schema(qcm.dict())
    msm.tol = tol
    sym_obj = Symtext.from_molecule(msm)
    return sym_obj


def symtext_is_chiral(sym_obj: Symtext) -> bool:
    """Based on its point group, is this molecule chiral?

    :param sym_obj: A MolSym symmetry object
    :type sym_obj: Symtext
    :return: `True` if it is, `False` if it isn't
    :rtype: bool
    """
    return sym_obj.pg.family in ("C", "D") and sym_obj.pg.subfamily is None


def symtext_external_symmetry_number(sym_obj: Symtext) -> int:
    """Get the external symmetry number of a molecule

    See Table II [here](https://cccbdb.nist.gov/thermo.asp) for a definition.

    :param sym_obj: A MolSym symmetry object
    :type sym_obj: Symtext
    :return: The external symmetry number
    :rtype: int
    """
    return sym_obj.rotational_symmetry_number
