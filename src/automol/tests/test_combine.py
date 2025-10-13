""" test combine
"""

import numpy
import pytest
from automol import combine, geom, graph

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


@pytest.mark.parametrize(
    "geo1,geo2,nfreqs1,nfreqs2",
    [
        (C2H6_GEO, OH_GEO, 18, 1),
        (C2H6_GEO, O_GEO, 18, 0),
        (OH_GEO, OH_GEO, 1, 1),
        (H_GEO, O_GEO, 0, 0),
    ],
)
def test__geo_combiners(geo1, geo2, nfreqs1, nfreqs2):
    """combine"""
    geos = (geo1, geo2)
    gras = list(map(geom.graph, geos))
    vdw_geo = combine.fake_vdw_geometry(*geos)
    vdw_gra = geom.graph(vdw_geo)
    assert vdw_gra == graph.union_from_sequence(gras, shift_keys=True)

    vdw_freqs = combine.fake_vdw_frequencies(*geos)
    print(vdw_freqs)
    num_atms = geom.count(vdw_geo)
    num_ext_modes = 5 if geom.is_linear(vdw_geo) else 6
    # Check that the total number of frequencies adds up to 3N-6 (or 3N-5 for linear)
    assert len(vdw_freqs) + nfreqs1 + nfreqs2 == 3 * num_atms - num_ext_modes


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
