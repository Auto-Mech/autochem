""" SMILES
"""
import automol.convert.smiles


def inchi(smi):
    """ SMILES => InChI
    """
    ich = automol.convert.smiles.inchi(smi)
    return ich
