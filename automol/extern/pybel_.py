""" pybel interface
"""

import pybel
from phydat import ptab
import automol.geom.base


# geometry
def to_geometry(pbm):
    """ Build an automol geometry data structure from a Pybel molecule object.

        :param pbm: Pybel molecule object
        :type pbm: Pybel object
        :rtype: automol geometry data structure
    """

    pbm.addh()
    pbm.make3D()
    nums = [atm.atomicnum for atm in pbm.atoms]
    symbs = list(map(ptab.to_symbol, nums))
    xyzs = tuple(tuple(atm.coords) for atm in pbm.atoms)
    geo = automol.geom.base.from_data(symbs, xyzs, angstrom=True)

    return geo


# inchi
def from_inchi(ich):
    """ Build a Pybel molecule object from an InChI string.

        :param ich: InChI string for a species
        :type ich: str
        :rtype: Pybel molecule object
    """
    return pybel.readstring('inchi', ich)


def to_inchi(pbm):
    """ Build an InChI string from a Pybel molecule object.

        :param pbm: Pybel molecule object
        :type pbm: Pybel object
        :rtype: str
    """
    return pbm.write('inchi').strip()


def from_smiles(smi):
    """ Build a Pybel molecule object from an SMILES string.

        :param smi: SMILES string for a species
        :type smi: str
        :rtype: Pybel molecule object
    """
    return pybel.readstring('smiles', smi)


def to_smiles(pbm):
    """ Build an SMILES string from a Pybel molecule object.

        :param pbm: Pybel molecule object
        :type pbm: Pybel object
        :rtype: str
    """
    return pbm.write('smiles').strip()
