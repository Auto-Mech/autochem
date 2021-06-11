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

HARM_FREQS = (
    269.39, 325.04, 429.05, 818.35, 913.01, 1084.71, 1116.71,
    1169.23, 1307.13, 1408.92, 1434.56, 1459.29, 1527.27,
    1531.68, 1555.98, 3020.12, 3062.32, 3121.51, 3142.58,
    3160.32, 3841.83)
PROJ_FREQS = (
    428.01, 815.69, 912.95, 1084.7, 1116.7, 1169.2, 1306.78,
    1408.87, 1434.54, 1459.25, 1527.25, 1531.57, 1555.97,
    3020.11, 3062.31, 3121.51, 3142.57, 3160.32, 3841.82)
TORS_FREQS = (263.39, 320.923)


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


def test__freq_anharm():
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


def test__freq_scale():
    """ test automol.prop.freq.
    """

    ref_scale_factor = ([], 1.0423855517115035)
    scale_factor = prop.freq.rotor_scale_factor_from_harmonics(
        HARM_FREQS, PROJ_FREQS, TORS_FREQS)

    assert not scale_factor[0]
    assert numpy.isclose(scale_factor[1], ref_scale_factor[1])

    # need tests where idx remove list nonempty
    # nce to have 1DHRFA from large molecule
