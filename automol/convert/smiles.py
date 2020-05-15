""" smiles conversions
"""
import automol.convert.inchi
from automol.convert import _rdkit


def inchi(smi):
    """ SMILES => InChI
    """
    ich = automol.convert.inchi.object_to_hardcoded_inchi_by_key(
        'smiles', smi, comp=_compare)
    if ich is None:
        rdm = _rdkit.from_smiles(smi)
        ich = _rdkit.to_inchi(rdm)
    return ich


# helpers
def _compare(smi1, smi2):
    return _canonicalize(smi1) == _canonicalize(smi2)


def _canonicalize(smi):
    return _rdkit.to_smiles(_rdkit.from_smiles(smi))
