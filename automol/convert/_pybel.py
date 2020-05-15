""" pybel interface
"""
import pybel
from qcelemental import periodictable as pt
import automol.create


# geometry
def to_geometry(pbm):
    """ cartesian geometry from a pybel molecule object
    """
    pbm.addh()
    pbm.make3D()
    nums = [atm.atomicnum for atm in pbm.atoms]
    syms = list(map(pt.to_E, nums))
    xyzs = tuple(tuple(atm.coords) for atm in pbm.atoms)
    geo = automol.create.geom.from_data(syms, xyzs, angstrom=True)
    return geo


# inchi
def from_inchi(ich):
    """ pybel molecule object from an InChI string
    """
    pbm = pybel.readstring('inchi', ich)
    return pbm


def to_inchi(pbm):
    """ InChI string from a pybel molecule object
    """
    ich = pbm.write('inchi').strip()
    return ich
