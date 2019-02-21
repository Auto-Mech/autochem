""" functions operating on SMILES strings
"""
from rdkit import RDLogger
import rdkit.Chem as _rd_chem

_LOGGER = RDLogger.logger()
_LOGGER.setLevel(RDLogger.ERROR)


def inchi(smi):
    """ InChI string from a SMILES string
    """
    rdm = _rdm_from_smiles(smi)
    ich = _rdm_to_inchi(rdm)
    return ich


def _rdm_from_smiles(smi):
    """ rdkit molecule object from a SMILES string
    """
    rdm = _rd_chem.MolFromSmiles(smi)
    assert rdm is not None
    return rdm


def _rdm_to_inchi(rdm):
    """ InChI string from an rdkit molecule object
    """
    ret = _rd_chem.inchi.MolToInchi(rdm)
    return ret
