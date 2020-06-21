""" SMILES
"""
import automol.convert.smiles


def inchi(smi):
    """ SMILES => InChI
    """
    # print('smi', smi)
    ich = automol.convert.smiles.inchi(smi)
    return ich
