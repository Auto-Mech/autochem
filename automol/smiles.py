""" SMILES
"""

import automol.convert.smiles


def inchi(smi):
    """ Convert a SMILES string into an InChI string.

        :param smi: SMILES string
        :type smi: str
        :rtype: str
    """
    return automol.convert.smiles.inchi(smi)
