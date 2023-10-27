""" test combine
"""

import numpy
from automol import combine, geom

C2H6_GEO = (
    ("C", (-1.4283035320563338, 0.013425343735546437, -0.030302158896694683)),
    ("C", (1.4283027358735494, -0.013425597530894248, 0.0303022919384165)),
    ("H", (-2.1972722614281355, -0.19229727219177065, 1.8778380427620682)),
    ("H", (-2.121310184939721, 1.792702413487708, -0.8231106338374065)),
    ("H", (-2.1448124562913287, -1.5396513482615042, -1.191852168914227)),
    ("H", (2.1448121742707795, 1.539654946791746, 1.1918517388178247)),
    ("H", (2.1972712765396953, 0.1922944277301287, -1.8778395029874426)),
    ("H", (2.121312248031497, -1.7927029137609576, 0.8231123911174519)),
)
OH_GEO = (
    ("O", (-0.3560094121, -0.0911704723, 0.0000000000)),
    ("H", (1.5152996925, -0.0911704723, 0.0000000000)),
)
O_GEO = (("O", (1.0000000000, 1.0000000000, 1.0000000000)),)
H_GEO = (("H", (0.0000000000, 0.0000000000, 0.0000000000)),)

C2H6_ELVLS = ((0.0, 1),)
OH_ELVLS = ((0.0, 2), (138.9, 2))
O_ELVLS = ((0.0, 5), (158.5, 3), (226.5, 1))


def test__geo_combiners():
    """combine"""

    ref_vdw_geo1 = (
        ("C", (-1.428303532056334, 0.01342534373554644, -0.03030215889669468)),
        ("C", (1.4283027358735494, -0.013425597530894248, 0.0303022919384165)),
        ("H", (-2.1972722614281355, -0.19229727219177065, 1.8778380427620682)),
        ("H", (-2.121310184939721, 1.792702413487708, -0.8231106338374065)),
        ("H", (-2.1448124562913287, -1.5396513482615042, -1.191852168914227)),
        ("H", (2.1448121742707795, 1.539654946791746, 1.1918517388178247)),
        ("H", (2.1972712765396953, 0.1922944277301287, -1.8778395029874426)),
        ("H", (2.121312248031497, -1.7927029137609576, 0.8231123911174519)),
        ("O", (-0.3560094121, -0.0911704723, 7.877838042762068)),
        ("H", (1.5152996925, -0.0911704723, 7.877838042762068)),
    )
    ref_vdw_geo2 = (
        ("C", (-1.428303532056334, 0.01342534373554644, -0.03030215889669468)),
        ("C", (1.4283027358735494, -0.013425597530894248, 0.0303022919384165)),
        ("H", (-2.1972722614281355, -0.19229727219177065, 1.8778380427620682)),
        ("H", (-2.121310184939721, 1.792702413487708, -0.8231106338374065)),
        ("H", (-2.1448124562913287, -1.5396513482615042, -1.191852168914227)),
        ("H", (2.1448121742707795, 1.539654946791746, 1.1918517388178247)),
        ("H", (2.1972712765396953, 0.1922944277301287, -1.8778395029874426)),
        ("H", (2.121312248031497, -1.7927029137609576, 0.8231123911174519)),
        ("O", (1.0, 1.0, 7.877838042762068)),
    )
    ref_vdw_geo3 = (("H", (0.0, 0.0, 0.0)), ("O", (1.0, 1.0, 6.0)))

    # May need to reverse for code coverage
    vdw_geo1 = combine.fake_vdw_geometry(C2H6_GEO, OH_GEO)
    vdw_geo2 = combine.fake_vdw_geometry(C2H6_GEO, O_GEO)
    vdw_geo3 = combine.fake_vdw_geometry(H_GEO, O_GEO)
    assert geom.almost_equal_dist_matrix(vdw_geo1, ref_vdw_geo1)
    assert geom.almost_equal_dist_matrix(vdw_geo2, ref_vdw_geo2)
    assert geom.almost_equal_dist_matrix(vdw_geo3, ref_vdw_geo3)

    ref_vdw_freqs1 = (30.0, 50.0, 70.0)
    ref_vdw_freqs2 = (30.0, 50.0)
    ref_vdw_freqs3 = (30.0,)
    # ref_vdw_freqs4 = (30.0,)
    ref_vdw_freqs5 = ()

    vdw_freqs1 = combine.fake_vdw_frequencies(C2H6_GEO, OH_GEO)
    vdw_freqs2 = combine.fake_vdw_frequencies(C2H6_GEO, O_GEO)
    vdw_freqs3 = combine.fake_vdw_frequencies(H_GEO, OH_GEO)
    # vdw_freqs4 = combine.fake_vdw_frequencies(OH_GEO, OH_GEO)
    vdw_freqs5 = combine.fake_vdw_frequencies(H_GEO, O_GEO)
    assert numpy.allclose(vdw_freqs1, ref_vdw_freqs1)
    assert numpy.allclose(vdw_freqs2, ref_vdw_freqs2)
    assert numpy.allclose(vdw_freqs3, ref_vdw_freqs3)
    # assert numpy.allclose(vdw_freqs4, ref_vdw_freqs4)
    assert numpy.allclose(vdw_freqs5, ref_vdw_freqs5)


def test__elec_levels_combiners():
    """test combine.electronic_energy_levels"""

    ref_elvls1 = ((0.0, 2), (138.9, 2))
    ref_elvls2 = ((0.0, 5), (158.5, 3), (226.5, 1))
    ref_elvls3 = (
        (0.0, 10),
        (138.9, 10),
        (158.5, 6),
        (226.5, 2),
        (297.4, 6),
        (365.4, 2),
    )
    ref_elvls4 = ((0.0, 4), (138.9, 8), (277.8, 4))

    elvls1 = combine.electronic_energy_levels(C2H6_ELVLS, OH_ELVLS)
    elvls2 = combine.electronic_energy_levels(C2H6_ELVLS, O_ELVLS)
    elvls3 = combine.electronic_energy_levels(OH_ELVLS, O_ELVLS)
    elvls4 = combine.electronic_energy_levels(OH_ELVLS, OH_ELVLS)

    assert numpy.allclose(elvls1, ref_elvls1)
    assert numpy.allclose(elvls2, ref_elvls2)
    assert numpy.allclose(elvls3, ref_elvls3)
    assert numpy.allclose(elvls4, ref_elvls4)


def test__formula_combiners():
    """test combine.formula_string"""

    assert combine.formula_string(C2H6_GEO, OH_GEO) == "C2H7O1"
