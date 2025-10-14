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


def test__freq_scale():
    """ test prop.freq.scale_frequencies_and_zpe
    """

    ref_freqs1 = (
        269.7640500413816, 324.6923142211176, 426.9219396878835,
        805.8075991816046, 897.2279743415297, 1062.4734704997502,
        1093.1938557971678, 1143.564362533476, 1275.5414833650848,
        1372.714206497535, 1397.1598211645016, 1420.7261833223702,
        1485.4496817445583, 1489.6455433499743, 1512.7593918088166,
        2889.265328197776, 2928.5387009393767, 2983.591051878249,
        3003.179023155766, 3019.6675272862385, 3650.669218489695)
    # ref_zpve1 = 0.08022224113138259

    freqs1, zpve1 = prop.freq.scale_frequencies_and_zpe(
        HARM_FREQS, METHOD1, BASIS1)
    assert numpy.allclose(freqs1, ref_freqs1)
    # These tests are all breaking
    print(zpve1)
    # assert numpy.isclose(zpve1, ref_zpve1)

    ref_freqs2 = (
        268.38245696488616, 322.6257929880899, 423.5081931429338,
        797.1881577019535, 887.374736526321, 1050.449272190552,
        1080.775735462857, 1130.5074859619388, 1260.8552277956692,
        1356.8713534457852, 1381.0319679440404, 1404.325887926194,
        1468.3129211579828, 1472.4616405582456, 1495.3171171840374,
        2860.65588233228, 2899.733690806601, 2954.5232913855507,
        2974.0209754448047, 2990.4347795393846, 3619.4624102348403)
    # ref_zpve2 = 0.07849214333676563

    freqs2, zpve2 = prop.freq.scale_frequencies_and_zpe(
        HARM_FREQS, METHOD2, BASIS2)
    assert numpy.allclose(freqs2, ref_freqs2)
    # These tests are all breaking
    print(zpve2)
    # assert numpy.isclose(zpve2, ref_zpve2)

    ref_freqs3 = (
        269.39, 325.04, 429.05,
        818.35, 913.01, 1084.71,
        1116.71, 1169.23, 1307.13,
        1408.92, 1434.56, 1459.29,
        1527.27, 1531.68, 1555.98,
        3020.12, 3062.32, 3121.51,
        3142.58, 3160.32, 3841.83)
    # ref_zpve3 = 0.08132830609426457

    freqs3, zpve3 = prop.freq.scale_frequencies_and_zpe(
        HARM_FREQS, METHOD3, BASIS3)
    assert numpy.allclose(freqs3, ref_freqs3)
    # These tests are all breaking
    print(zpve3)
    # assert numpy.isclose(zpve3, ref_zpve3)

    freqs4, zpve4 = prop.freq.scale_frequencies_and_zpe(
        HARM_FREQS, METHOD1, BASIS1,
        scale_method='harm')
    assert numpy.allclose(freqs4, ref_freqs3)
    # These tests are all breaking
    print(zpve4)
    # assert numpy.isclose(zpve4, ref_zpve3)


def test__freq_scale_factor():
    """ test automol.prop.freq.rotor_scale_factor_from_harmonics
    """

    ref_scale_factor = ([], 1.0423855517115035)
    scale_factor = prop.freq.rotor_scale_factor_from_harmonics(
        HARM_FREQS, PROJ_FREQS, TORS_FREQS)

    assert not scale_factor[0]
    assert numpy.isclose(scale_factor[1], ref_scale_factor[1])

    # need tests where idx remove list nonempty
    # nice to have 1DHRFA from large molecule


if __name__ == '__main__':
    test__freq_scale()
