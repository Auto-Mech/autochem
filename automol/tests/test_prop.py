""" test automol.prop
"""

import numpy
from automol import prop


# Wavefunction Props
DIP_MOM_VEC = (1.00, 2.00, 3.00)
POLAR_TENSOR = ((1.12, 2.09, 3.21),
                (4.45, 5.86, 6.63),
                (7.53, 8.24, 9.97))

# Frequency Props
FREQS = (100, 200, 300, 400, 500)
METHOD1, BASIS1 = ('b2plypd3', 'cc-pvtz')
METHOD2, BASIS2 = ('wb97xd', '6-31g*')
METHOD3, BASIS3 = ('b3lyp', '6-31g*')


def test__dipole_moment():
    """ test prop.total_dipole_moment
        test prop.total_polarizability
    """

    ref_dip_mom = 3.74166

    dip_mom = prop.total_dipole_moment(DIP_MOM_VEC)
    assert numpy.isclose(ref_dip_mom, dip_mom)

    ref_polar = 5.65

    polar = prop.total_polarizability(POLAR_TENSOR)
    assert numpy.isclose(ref_polar, polar)


def test__freqs():
    """ test prop.freq.anharm_by_scaling
    """

    ref_freqs1 = (102.92272, 203.95523, 303.94752, 403.15841, 501.72796)
    ref_azpve1 = 0.0034083028828633282

    freqs1, azpve1 = prop.freq.anharm_by_scaling(FREQS, METHOD1, BASIS1)
    assert numpy.allclose(freqs1, ref_freqs1)
    assert numpy.isclose(azpve1, ref_azpve1)

    ref_freqs2 = (101.53824, 200.41326, 298.24476, 395.38958, 492.01924)
    ref_azpve2 = 0.0034243108652987474

    freqs2, azpve2 = prop.freq.anharm_by_scaling(FREQS, METHOD2, BASIS2)
    assert numpy.allclose(freqs2, ref_freqs2)
    assert numpy.isclose(azpve2, ref_azpve2)

    ref_freqs3 = (100.0, 200.0, 300.0, 400.0, 500.0)
    ref_azpve3 = 0.0034172514395752504

    freqs3, azpve3 = prop.freq.anharm_by_scaling(FREQS, METHOD3, BASIS3)
    assert numpy.allclose(freqs3, ref_freqs3)
    assert numpy.isclose(azpve3, ref_azpve3)
