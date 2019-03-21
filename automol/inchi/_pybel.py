""" pybel interface
"""
import pybel
from ..constructors.geom import from_data as _geom_from_data


def from_inchi(ich):
    """ pybel molecule object from an InChI string
    """
    pbm = pybel.readstring('inchi', ich)
    return pbm


def geometry(pbm):
    """ cartesian geometry from a pybel molecule object
    """
    pbm.addh()
    pbm.make3D()
    anbs = [atm.atomicnum for atm in pbm.atoms]
    syms = list(map(_atomic_symbol, anbs))
    xyzs = tuple(tuple(atm.coords) for atm in pbm.atoms)
    geo = _geom_from_data(syms, xyzs, angstrom=True)
    return geo


def _atomic_symbol(anum):
    """ convert atomic number to atomic symbol
    """
    anum2asymb = ['X', 'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE',
                  'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA',
                  'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
                  'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR',
                  'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN',
                  'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND',
                  'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB',
                  'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
                  'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH',
                  'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM',
                  'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS',
                  'RG', 'UUB', 'UUT', 'UUQ', 'UUP', 'UUH', 'UUS', 'UUO']
    return anum2asymb[anum].title()
