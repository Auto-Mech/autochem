"""
    test automol.rotor
"""

import numpy
import automol


# Potentials
BAD_POT = {
    (0.00000000,): 0.00, (0.52359878, 1.0): 0.77, (1.04719755,): 1.62
}
POT1 = {
    (0.00000000,): 0.00, (0.52359878,): 0.77, (1.04719755,): 1.62,
    (1.57079633,): 0.90, (2.09439510,): 0.01, (2.61799388,): 0.61,
    (3.14159265,): 1.50, (3.66519143,): 0.99, (4.18879020,): 0.32,
    (4.71238898,): 0.89, (5.23598776,): 1.52, (5.75958653,): 0.74
}
POT2 = {
    (0.00000000,): 0.00, (0.52359878,): 1.74, (1.04719755,): 3.58,
    (1.57079633,): 1.68, (2.09439510,): 0.01, (2.61799388,): 1.75,
    (3.14159265,): 3.59, (3.66519143,): 1.69, (4.18879020,): 0.02,
    (4.71238898,): 1.72, (5.23598776,): 3.60, (5.75958653,): 1.60
}
POT3 = {
    (1.0, 0.1): 1.1, (1.0, 0.2): 1.2,
    (2.0, 0.1): 1.3, (2.0, 0.2): 1.4,
    (3.0, 0.1): 1.5, (3.0, 0.2): 1.6,
    (4.0, 0.1): 1.7, (4.0, 0.2): 1.8
}
PCOORDS1 = (1.00, 2.00, 3.00, 4.00)
PCOORDS2 = (0.10, 0.20)
SYM_NUM1 = 1
SYM_NUM2 = 3
SCALE_COEFF = 1.25
NUM_TORS = 3
# SCAN_INCREMENT = 0.5
SCAN_INCREMENT = 0.523599

# Scaling potential


def test__valid_potential():
    """ test automol.rotor.pot.valid
    """

    assert automol.rotor.pot.valid(POT1)
    assert not automol.rotor.pot.valid(BAD_POT)


def test__build_potential():
    """ test automol.rotor.pot.points
        test automol.rotor.pot.coords
    """

    ref_grid_pts = ((0, 0), (0, 1),
                    (1, 0), (1, 1),
                    (2, 0), (2, 1),
                    (3, 0), (3, 1))
    ref_grid_coords = ((1.0, 0.1), (1.0, 0.2),
                       (2.0, 0.1), (2.0, 0.2),
                       (3.0, 0.1), (3.0, 0.2),
                       (4.0, 0.1), (4.0, 0.2))

    assert automol.rotor.pot.points((PCOORDS1, PCOORDS2)) == ref_grid_pts
    assert automol.rotor.pot.coords((PCOORDS1, PCOORDS2)) == ref_grid_coords


def test__transform_potential():
    """ test automol.rotor.pot.scale
        test automol.rotor.pot.truncate
    """

    # Test scaling
    ref_pot_scaled = {(0.0,): 0.0,
                      (0.52359878,): 0.89350585047046,
                      (1.04719755,): 1.8798434776131758,
                      (1.57079633,): 1.0443574875628754,
                      (2.0943951,): 0.011603972084031949,
                      (2.61799388,): 0.7078422971259488,
                      (3.14159265,): 1.7405958126047922,
                      (3.66519143,): 1.1487932363191629,
                      (4.1887902,): 0.37132710668902236,
                      (4.71238898,): 1.0327535154788434,
                      (5.23598776,): 1.7638037567728562,
                      (5.75958653,): 0.8586939342183642}

    pot_scaled = automol.rotor.pot.scale(POT1, SCALE_COEFF, NUM_TORS)

    assert numpy.allclose(list(pot_scaled.keys()), list(ref_pot_scaled.keys()))
    for key, val in pot_scaled.items():
        assert numpy.isclose(val, ref_pot_scaled[key])

    # Test truncating
    ref_pot_trunc1 = {(0.00000000,): 0.00, (0.52359878,): 0.77,
                      (1.04719755,): 1.62, (1.57079633,): 0.90,
                      (2.09439510,): 0.01, (2.61799388,): 0.61,
                      (3.14159265,): 1.50, (3.66519143,): 0.99,
                      (4.18879020,): 0.32, (4.71238898,): 0.89,
                      (5.23598776,): 1.52, (5.75958653,): 0.74}
    ref_pot_trunc2 = {(0.00000000,): 0.00, (0.52359878,): 1.74,
                      (1.04719755,): 3.58, (1.57079633,): 1.68}

    pot_trunc1 = automol.rotor.pot.truncate(POT1, SYM_NUM1)
    pot_trunc2 = automol.rotor.pot.truncate(POT2, SYM_NUM2)

    assert numpy.allclose(list(pot_trunc1.keys()), list(ref_pot_trunc1.keys()))
    assert numpy.allclose(list(pot_trunc2.keys()), list(ref_pot_trunc2.keys()))
    for key, val in pot_trunc1.items():
        assert numpy.isclose(val, ref_pot_trunc1[key])
    for key, val in pot_trunc2.items():
        assert numpy.isclose(val, ref_pot_trunc2[key])

    ref_idx_pot1 = {(0,): 0.00, (1,): 0.77, (2,): 1.62, (3,): 0.90,
                    (4,): 0.01, (5,): 0.61, (6,): 1.50, (7,): 0.99,
                    (8,): 0.32, (9,): 0.89, (10,): 1.52, (11,): 0.74}
    ref_idx_pot2 = {(0, 0): 1.1, (0, 1): 1.2,
                    (1, 0): 1.3, (1, 1): 1.4,
                    (2, 0): 1.5, (2, 1): 1.6,
                    (3, 0): 1.7, (3, 1): 1.8}
    idx_pot1 = automol.rotor.pot.by_index(POT1)
    idx_pot2 = automol.rotor.pot.by_index(POT3)

    assert numpy.allclose(list(idx_pot1.keys()), list(ref_idx_pot1.keys()))
    assert numpy.allclose(list(idx_pot2.keys()), list(ref_idx_pot2.keys()))
    for key, val in idx_pot1.items():
        assert numpy.isclose(val, ref_idx_pot1[key])
    for key, val in idx_pot2.items():
        assert numpy.isclose(val, ref_idx_pot2[key])


if __name__ == '__main__':
    test__valid_potential()
    test__build_potential()
    test__transform_potential()
