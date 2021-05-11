""" smiles conversions
"""

from automol.convert.inchi import object_to_hardcoded_inchi_by_key
from automol.convert import _rdkit


def inchi(smi):
    """ Convert a SMILES string into an InChI string.

        :param smi: SMILES string
        :type smi: str
        :rtype: str
    """

    ich = object_to_hardcoded_inchi_by_key(
        'smiles', smi, comp=_compare)
    if ich is None:
        rdm = _rdkit.from_smiles(smi)
        ich = _rdkit.to_inchi(rdm)
    return ich


# helpers
def _compare(smi1, smi2):
    """ Check if two SMILES strings are similar.

        :param smi1: SMILES string 1
        :type smi1: str
        :param smi2: SMILES string 2
        :type smi2: str
        :rtype: bool
    """
    return _canonicalize(smi1) == _canonicalize(smi2)


def _canonicalize(smi):
    """ Convert a SMILES string into its canonical form.

        :param smi: SMILES string
        :type smi: str
        :rtype: str
    """
    return _rdkit.to_smiles(_rdkit.from_smiles(smi))
