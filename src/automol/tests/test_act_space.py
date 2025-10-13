""" simple test to get active space
"""

from phydat import act_space


def test__():
    """ test read act_space dct
    """

    assert act_space.DCT.get(('InChI=1S/O2/c1-2', 3), None) == (4, 6, 1)
    assert act_space.DCT.get(('InChI=1S/O2/c1-2', 1), None) is None
