"""Interface to the MolSym code
"""
import qcelemental as qcel
from automol.geom.base import xyz_string
from molsym.molecule import Molecule
from molsym.symtext.symtext import Symtext


def symtext_from_geometry(geo) -> Symtext:
    """Generate a MolSym symmetry object from a geometry

    :param geo: A geometry
    :type geo: automol geom data structure
    :return: A MolSym symmetry object
    :rtype: Symtext
    """
    qcm = qcel.models.Molecule.from_data(xyz_string(geo))
    msm = Molecule.from_schema(qcm.dict())
    sym_obj = Symtext.from_molecule(msm)
    return sym_obj


def symtext_external_symmetry_number(sym_obj: Symtext) -> int:
    """Get the external symmetry number of a molecule

    See Table II [here](https://cccbdb.nist.gov/thermo.asp) for a definition.

    :param sym_obj: A MolSym symmetry object
    :type sym_obj: Symtext
    :return: The external symmetry number
    :rtype: int
    """
    return sym_obj.rotational_symmetry_number
