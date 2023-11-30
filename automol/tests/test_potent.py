""" test automol.data.potent module
"""

from automol.data import potent

# Potentials
POT1_DCT = {
    (0.00000000,): 0.00,
    (0.52359878,): 0.77,
    (1.04719755,): 1.62,
    (1.57079633,): 0.90,
    (2.09439510,): 0.01,
    (2.61799388,): 0.61,
    (3.14159265,): 1.50,
    (3.66519143,): 0.99,
    (4.18879020,): 0.32,
    (4.71238898,): 0.89,
    (5.23598776,): 1.52,
    (5.75958653,): 0.74,
}
POT1_GEO_DCT = {
    (0.0,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (0.016867, 0.751572, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (0.52359878,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-0.402803, 0.43847, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (1.04719755,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-0.822473, 0.125368, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (1.57079633,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-1.242143, -0.187733, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (2.0943951,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-1.661813, -0.500835, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (2.61799388,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-2.081483, -0.813937, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (3.14159265,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-2.501153, -1.127039, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (3.66519143,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-2.920823, -1.440141, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (4.1887902,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-3.340493, -1.753243, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (4.71238898,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-3.760163, -2.066345, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (5.23598776,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-4.179833, -2.379447, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
    (5.75958653,): (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-4.599503, -2.692549, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    ),
}

POT2_DCT = {
    (0.00000000,): 0.00,
    (0.52359878,): 1.74,
    (1.04719755,): 3.58,
    (1.57079633,): 1.68,
    (2.09439510,): 0.01,
    (2.61799388,): 1.75,
    (3.14159265,): 3.59,
    (3.66519143,): 1.69,
    (4.18879020,): 0.02,
    (4.71238898,): 1.72,
    (5.23598776,): 3.60,
    (5.75958653,): 1.60,
}
POT3_DCT = {
    (0.00000000,): 0.00,
    (0.52359878,): None,
    (1.04719755,): 3.58,
    (1.57079633,): None,
    (2.09439510,): 0.01,
    (2.61799388,): 1.75,
    (3.14159265,): 3.59,
    (3.66519143,): 1.69,
    (4.18879020,): 0.02,
    (4.71238898,): 1.72,
    (5.23598776,): 3.60,
    (5.75958653,): None,
}
POT4_DCT = {
    (0.00000000,): None,
    (0.52359878,): None,
    (1.04719755,): None,
    (1.57079633,): None,
    (2.09439510,): None,
    (2.61799388,): None,
    (3.14159265,): None,
    (3.66519143,): None,
    (4.18879020,): None,
    (4.71238898,): None,
    (5.23598776,): None,
    (5.75958653,): None,
}
POT5_DCT = {
    (1.0, 0.1): 1.1,
    (1.0, 0.2): 1.2,
    (2.0, 0.1): 1.3,
    (2.0, 0.2): 1.4,
    (3.0, 0.1): 1.5,
    (3.0, 0.2): 1.6,
    (4.0, 0.1): 1.7,
    (4.0, 0.2): 1.8,
}
POT6_DCT = {
    (0.00000000,): 0.001,
    (0.52359878,): 0.77,
    (1.04719755,): 1.62,
    (1.57079633,): -4.0,
    (2.09439510,): 45.0,
    (2.61799388,): 55.0,
    (3.14159265,): 1.50,
    (3.66519143,): 0.99,
    (4.18879020,): 600.1,
    (4.71238898,): 0.89,
    (5.23598776,): -6.0,
    (5.75958653,): 0.74,
}


def test__potential():
    """test potent.from_dict, potent.scaled, potent.dict_"""
    pot1 = potent.from_dict(POT1_DCT, ["D5"])
    assert potent.value(pot1, D5=0.523) == 0.77
    assert potent.value(pot1, 0.523) == 0.77

    pot5 = potent.from_dict(POT5_DCT, ["D5", "D8"])
    assert potent.value(pot5, D5=3.0, D8=0.2) == 1.6
    assert potent.value(pot5, 3.0, 0.2) == 1.6

    # Test scaling
    ref_pot1_scaled = potent.from_dict(
        {
            (0.0,): 0.0,
            (0.52359878,): 0.9625,
            (1.04719755,): 2.0250,
            (1.57079633,): 1.1250,
            (2.0943951,): 0.0125,
            (2.61799388,): 0.7625,
            (3.14159265,): 1.8750,
            (3.66519143,): 1.2375,
            (4.1887902,): 0.4000,
            (4.71238898,): 1.1125,
            (5.23598776,): 1.900,
            (5.75958653,): 0.9250,
        },
        ["D5"],
    )

    pot1_scaled = potent.scale(pot1, 1.25)
    assert potent.almost_equal(pot1_scaled, ref_pot1_scaled)

    ref_pot1_idx_dct = {
        (0,): 0.00,
        (1,): 0.77,
        (2,): 1.62,
        (3,): 0.90,
        (4,): 0.01,
        (5,): 0.61,
        (6,): 1.50,
        (7,): 0.99,
        (8,): 0.32,
        (9,): 0.89,
        (10,): 1.52,
        (11,): 0.74,
    }
    assert potent.dict_(pot1, index=True) == ref_pot1_idx_dct

    ref_pot5_idx_dct = {
        (0, 0): 1.1,
        (0, 1): 1.2,
        (1, 0): 1.3,
        (1, 1): 1.4,
        (2, 0): 1.5,
        (2, 1): 1.6,
        (3, 0): 1.7,
        (3, 1): 1.8,
    }
    assert potent.dict_(pot5, index=True) == ref_pot5_idx_dct


def test__potential_with_geom():
    """test potential with geometries"""
    pot1 = potent.from_dict(POT1_DCT, aux_dct_dct={"geom": POT1_GEO_DCT})

    assert potent.keys(pot1) == ("energy", "geom")
    assert potent.value(pot1, 0.523) == 0.77
    assert potent.value(pot1, 0.523, key="geom") == (
        ("O", (0.016867, 0.751572, 0.0)),
        ("H", (-0.402803, 0.43847, 0.0)),
        ("H", (1.433949, -0.408156, 0.0)),
    )
    assert potent.dict_(pot1) == POT1_DCT
    assert potent.dict_(pot1, key="geom") == POT1_GEO_DCT

    pot1_without_geos = potent.from_dict(POT1_DCT)

    assert potent.keys(pot1_without_geos) == ("energy",)
    assert potent.value(pot1_without_geos, 0.523) == 0.77
    assert potent.value(pot1_without_geos, 0.523, key="geom") is None
    assert potent.dict_(pot1_without_geos) == POT1_DCT
    assert potent.dict_(pot1_without_geos, key="geom") is None


def test__has_defined_values():
    """test potent.has_defined_values and drop_null flag"""
    pot1 = potent.from_dict(POT1_DCT)
    assert potent.has_defined_values(pot1)

    pot4 = potent.from_dict(POT4_DCT)
    assert not potent.has_defined_values(pot4)

    ref_filt_pot3_dct = {
        (0.00000000,): 0.00,
        (1.04719755,): 3.58,
        (2.09439510,): 0.01,
        (2.61799388,): 1.75,
        (3.14159265,): 3.59,
        (3.66519143,): 1.69,
        (4.18879020,): 0.02,
        (4.71238898,): 1.72,
        (5.23598776,): 3.60,
    }
    pot3 = potent.from_dict(POT3_DCT)
    assert potent.dict_(pot3, drop_null=True) == ref_filt_pot3_dct


def test__clean():
    """test potent.clean"""
    ref_clean_pot6_idx_dct = {
        (0,): 0.0,
        (1,): 0.77,
        (2,): 1.62,
        (3,): 0.0,
        (4,): 45.0,
        (5,): 50.0,
        (6,): 1.5,
        (7,): 0.99,
        (9,): 0.89,
        (11,): 0.74,
    }

    pot6 = potent.from_dict(POT6_DCT)
    clean_pot6 = potent.clean(pot6)
    clean_pot6_idx_dct = potent.dict_(clean_pot6, index=True, drop_null=True)
    print(clean_pot6_idx_dct)
    assert clean_pot6_idx_dct == ref_clean_pot6_idx_dct


def test__dict__zero_start_coord():
    """test potent.dict with zero_start_coord flag"""
    ref_zero_pot5_dct = {
        (0.0, 0.0): 1.1,
        (0.0, 0.1): 1.2,
        (1.0, 0.0): 1.3,
        (1.0, 0.1): 1.4,
        (2.0, 0.0): 1.5,
        (2.0, 0.1): 1.6,
        (3.0, 0.0): 1.7,
        (3.0, 0.1): 1.8,
    }

    pot5 = potent.from_dict(POT5_DCT)
    assert potent.dict_(pot5) == POT5_DCT
    assert potent.dict_(pot5, zero_start_coord=True) == ref_zero_pot5_dct


if __name__ == "__main__":
    # test__potential()
    # test__potential_with_geom()
    # test__has_defined_values()
    test__clean()
    test__dict__zero_start_coord()
