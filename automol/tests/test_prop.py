""" test automol.geom
"""

import numpy
from automol import prop


def test__dipole_moment():
    """ test prop.total_dipole_moment
    """

    ref_dip_mom = 3.74166

    dip_mom_vec = (1.00, 2.00, 3.00)
    dip_mom = prop.total_dipole_moment(dip_mom_vec)

    assert numpy.isclose(ref_dip_mom, dip_mom)


def test__polarizability():
    """ test prop.total_polarizability
    """

    ref_polar = 5.65

    polar_tensor = ((1.12, 2.09, 3.21),
                    (4.45, 5.86, 6.63),
                    (7.53, 8.24, 9.97))
    polar = prop.total_polarizability(polar_tensor)

    assert numpy.isclose(ref_polar, polar)


if __name__ == '__main__':
    test__dipole_moment()
    test__polarizability()
