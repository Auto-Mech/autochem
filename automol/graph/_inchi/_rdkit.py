""" rdkit interface
"""
from rdkit import RDLogger
import rdkit.Chem as _rd_chem

_LOGGER = RDLogger.logger()
_LOGGER.setLevel(RDLogger.ERROR)


def from_molfile(mfl):
    """ rdkit molecule object from a mol block string
    """
    rdm = _rd_chem.rdmolfiles.MolFromMolBlock(mfl, removeHs=False)
    assert rdm is not None
    return rdm


def to_inchi_with_aux_info(rdm):
    """ InChI string from an rdkit molecule object
    """
    ich, ich_aux = _rd_chem.inchi.MolToInchiAndAuxInfo(rdm)
    return ich, ich_aux
