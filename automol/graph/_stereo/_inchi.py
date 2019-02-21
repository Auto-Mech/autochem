""" stereo graph => InChI string conversion
"""
from .._inchi import atom_inchi_numbers
from .._inchi import stereo_inchi_from_coordinates
from ._stereo_ import _explicit_stereo
from ._stereo_ import _is_incomplete_or_higher_order
from ._intco import atom_stereo_coordinates as _atom_stereo_coordinates


def stereo_inchi(sgr):
    """ InChI string of this stereo graph
    """
    assert not _is_incomplete_or_higher_order(sgr)

    sgr = _explicit_stereo(sgr)
    atm_ich_num_dct = atom_inchi_numbers(sgr)
    atm_xyz_dct = _atom_stereo_coordinates(sgr, atm_ich_num_dct)
    ich = stereo_inchi_from_coordinates(sgr, atm_xyz_dct)
    return ich
