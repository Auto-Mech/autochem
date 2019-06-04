""" smiles conversions
"""
import automol.convert.inchi
from automol.convert import _rdkit

HARDCODED_SMILES_TO_INCHI_DCT = dict(
    map(reversed, automol.convert.inchi.HARDCODED_INCHI_TO_SMILES_DCT.items())
)


def inchi(smi):
    """ SMILES => InChI
    """

    # first, canonicalize the SMILES string
    rdm = _rdkit.from_smiles(smi)
    smi = _rdkit.to_smiles(rdm)

    if smi in HARDCODED_SMILES_TO_INCHI_DCT:
        ich = HARDCODED_SMILES_TO_INCHI_DCT[smi]
    else:
        ich = _rdkit.to_inchi(rdm)
    return ich
