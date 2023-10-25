""" test pot_
"""

import numpy
from automol import pot as pot_


# Potentials
POT1 = {
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
POT2 = {
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
POT3 = {
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
POT4 = {
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
POT5 = {
    (1.0, 0.1): 1.1,
    (1.0, 0.2): 1.2,
    (2.0, 0.1): 1.3,
    (2.0, 0.2): 1.4,
    (3.0, 0.1): 1.5,
    (3.0, 0.2): 1.6,
    (4.0, 0.1): 1.7,
    (4.0, 0.2): 1.8,
}
PCOORDS1 = (1.00, 2.00, 3.00, 4.00)
PCOORDS2 = (0.10, 0.20)
SYM_NUM1 = 1
SYM_NUM2 = 3
SCALE_COEFF = 1.25
NUM_TORS = 3
# SCAN_INCREMENT = 0.5
SCAN_INCREMENT = 0.523599


def test__transform_potential():
    """test pot_.scale
    test pot_.truncate
    """

    # Test scaling
    ref_pot_scaled = {
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
    }

    pot_scaled = pot_.scale(POT1, SCALE_COEFF)

    assert numpy.allclose(list(pot_scaled.keys()), list(ref_pot_scaled.keys()))
    for key, val in pot_scaled.items():
        assert numpy.isclose(val, ref_pot_scaled[key])

    ref_idx_pot1 = {
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
    ref_idx_pot2 = {
        (0, 0): 1.1,
        (0, 1): 1.2,
        (1, 0): 1.3,
        (1, 1): 1.4,
        (2, 0): 1.5,
        (2, 1): 1.6,
        (3, 0): 1.7,
        (3, 1): 1.8,
    }
    idx_pot1 = pot_.by_index(POT1)
    idx_pot2 = pot_.by_index(POT5)

    assert numpy.allclose(list(idx_pot1.keys()), list(ref_idx_pot1.keys()))
    assert numpy.allclose(list(idx_pot2.keys()), list(ref_idx_pot2.keys()))
    for key, val in idx_pot1.items():
        assert numpy.isclose(val, ref_idx_pot1[key])
    for key, val in idx_pot2.items():
        assert numpy.isclose(val, ref_idx_pot2[key])


def test__empty_terms_in_potential():
    """test pot_.is_nonempty
    test pot_.remove_empty_terms
    """

    assert pot_.has_defined_values(POT1)
    assert not pot_.has_defined_values(POT4)

    ref_filt_pot = {
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

    filt_pot = pot_.remove_empty_terms(POT3)
    assert numpy.allclose(
        tuple(filt_pot.keys()), tuple(ref_filt_pot.keys()), atol=1.0e-2
    )
    assert numpy.allclose(
        tuple(filt_pot.values()), tuple(ref_filt_pot.values()), atol=1.0e-2
    )
