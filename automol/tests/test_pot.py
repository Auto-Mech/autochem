""" test automol.intmol
"""

import numpy
import automol.pot


PROP_GEO1 = (
    ('C', (0.300662110143701, -0.45146519971545696, 0.5208682143568093)),
    ('H', (2.399149254920804, -0.3946288556825769, 0.536136790057131)),
    ('H', (-0.3993711916135264, 0.3895083230039884, -1.270731253238608)),
    ('H', (-0.4141144236399435, 0.6958184996927619, 2.125824819960968)),
    ('C', (-0.6303238583870289, -3.1678122617517217, 0.7839727870645417)),
    ('C', (0.3195573718790659, -4.813937448240943, -1.4016053045083412)),
    ('H', (-2.7315133185155234, -3.1653742416065107, 0.8049951600535731)),
    ('H', (0.06149708083384194, -3.947966968295672, 2.6083048214650737)),
    ('H', (-0.3722645227873339, -4.03378323264784, -3.2259372098405787)),
    ('C', (-0.6114269677077436, -7.530284973071136, -1.138499771441792)),
    ('H', (2.4207468066852305, -4.816374666386387, -1.4226284218604937)),
    ('H', (-2.70991407147779, -7.587122702273266, -1.153766809660938)),
    ('H', (0.08860820960266315, -8.371257730640474, 0.6530993359718258)),
    ('H', (0.10334927997136532, -8.677568509584962, -2.7434566228993194)))

PROP_GEO2 = (
    ('C', (0.5878206563471623, -0.19038043073027486, 1.0790169922039132)),
    ('H', (2.527962725012789, -0.5333845133705444, 1.8007522712353692)),
    ('H', (0.6391634644881258, 1.3170668967727561, -0.37719809900217977)),
    ('H', (-0.6121913407605588, 0.48983567918922677, 2.6593019822570736)),
    ('C', (-0.5345745789805177, -2.625274265384683, 0.0020410359294053986)),
    ('C', (1.07484339826323, -3.8096380986390974, -2.1041260920708122)),
    ('H', (-2.468920238734059, -2.229786545618703, -0.7158942226940402)),
    ('H', (-0.6907480986136927, -4.025052808686176, 1.561490570717891)),
    ('H', (0.11890751572678643, -5.5711110776432635, -2.7358332265801693)),
    ('C', (1.4145114676277863, -2.0893978516196934, -4.400589986011689)),
    ('H', (2.965904849365761, -4.327978946599484, -1.3499854232110824)),
    ('H', (-0.4512464979315055, -1.5355703545342925, -5.184753837615598)),
    ('H', (2.498511168547862, -0.36848779079717586, -3.891669130446078)),
    ('H', (2.4953702095902845, -3.1096562386285123, -5.8807421730423775)))


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


def test__grid():
    """ test automol.pot.grid
    """


def test__valid_potential():
    """ test automol.pot.valid
    """

    assert automol.pot.valid(POT1)
    assert not automol.pot.valid(BAD_POT)


def test__build_potential():
    """ test automol.pot.points
        test automol.pot.coords
    """

    ref_1dgrid_pts = ((0,), (1,), (2,), (3,))
    ref_1dgrid_coords = ((1.0,), (2.0,), (3.0,), (4.0,))

    ref_2dgrid_pts = ((0, 0), (0, 1),
                      (1, 0), (1, 1),
                      (2, 0), (2, 1),
                      (3, 0), (3, 1))
    ref_2dgrid_coords = ((1.0, 0.1), (1.0, 0.2),
                         (2.0, 0.1), (2.0, 0.2),
                         (3.0, 0.1), (3.0, 0.2),
                         (4.0, 0.1), (4.0, 0.2))

    assert automol.pot.points((PCOORDS1,)) == ref_1dgrid_pts
    assert automol.pot.coords((PCOORDS1,)) == ref_1dgrid_coords
    assert automol.pot.points((PCOORDS1, PCOORDS2)) == ref_2dgrid_pts
    assert automol.pot.coords((PCOORDS1, PCOORDS2)) == ref_2dgrid_coords


def test__transform_potential():
    """ test automol.pot.scale
        test automol.pot.truncate
    """

    # Test scaling
    ref_pot_scaled = {(0.0,): 0.0,
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
                      (5.75958653,): 0.9250}

    pot_scaled = automol.pot.scale(POT1, SCALE_COEFF)

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

    pot_trunc1 = automol.pot.truncate(POT1, SYM_NUM1)
    pot_trunc2 = automol.pot.truncate(POT2, SYM_NUM2)

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
    idx_pot1 = automol.pot.by_index(POT1)
    idx_pot2 = automol.pot.by_index(POT3)

    assert numpy.allclose(list(idx_pot1.keys()), list(ref_idx_pot1.keys()))
    assert numpy.allclose(list(idx_pot2.keys()), list(ref_idx_pot2.keys()))
    for key, val in idx_pot1.items():
        assert numpy.isclose(val, ref_idx_pot1[key])
    for key, val in idx_pot2.items():
        assert numpy.isclose(val, ref_idx_pot2[key])


def test__fit_potential():
    """ test automol.pot.fit_1d_potential
    """

    # full potential
    init_pot1 = {(0,): 0.000, (1,): 1.299, (2,): 3.085, (3,): 2.780,
                 (4,): 2.045, (5,): 3.052, (6,): 3.949, (7,): 2.655,
                 (8,): 1.480, (9,): 2.358, (10,): 2.948, (11,): 1.289}
    # one negative
    init_pot2 = {(0,): 0.000, (1,): 1.299, (2,): 3.085, (3,): 2.780,
                 (4,): 2.045, (5,): -10.0, (6,): 3.949, (7,): 2.655,
                 (8,): 1.480, (9,): 2.358, (10,): 2.948, (11,): 1.289}
    # two negatives
    init_pot3 = {(0,): 0.000, (1,): 1.299, (2,): 3.085, (3,): 2.780,
                 (4,): 2.045, (5,): -10.0, (6,): 3.949, (7,): 2.655,
                 (8,): 1.480, (9,): 2.358, (10,): 2.948, (11,): -10.0}
    # need one with max threshold closed

    pot1 = automol.pot.fit_1d_potential(
        init_pot1, min_thresh=-0.0001, max_thresh=50.0)
    pot2 = automol.pot.fit_1d_potential(
        init_pot2, min_thresh=-0.0001, max_thresh=50.0)
    pot3 = automol.pot.fit_1d_potential(
        init_pot3, min_thresh=-0.0001, max_thresh=50.0)

    ref_pot1 = {(0,): 0.000, (1,): 1.299, (2,): 3.085, (3,): 2.780,
                (4,): 2.045, (5,): 3.052, (6,): 3.949, (7,): 2.655,
                (8,): 1.480, (9,): 2.358, (10,): 2.948, (11,): 1.289}
    ref_pot2 = {(0,): 0.000, (1,): 1.299, (2,): 3.085, (3,): 2.780,
                (4,): 2.045, (5,): 3.091, (6,): 3.949, (7,): 2.655,
                (8,): 1.480, (9,): 2.358, (10,): 2.948, (11,): 1.289}
    ref_pot3 = {(0,): 0.000, (1,): 1.299, (2,): 3.085, (3,): 2.780,
                (4,): 2.045, (5,): 3.089, (6,): 3.949, (7,): 2.655,
                (8,): 1.480, (9,): 2.358, (10,): 2.948, (11,): 2.273}

    assert tuple(pot1.keys()) == tuple(ref_pot1.keys())
    assert numpy.allclose(
        tuple(pot1.values()), tuple(ref_pot1.values()), atol=1.0e-2)
    assert tuple(pot2.keys()) == tuple(ref_pot2.keys())
    assert numpy.allclose(
        tuple(pot2.values()), tuple(ref_pot2.values()), atol=1.0e-2)
    assert tuple(pot3.keys()) == tuple(ref_pot3.keys())
    assert numpy.allclose(
        tuple(pot3.values()), tuple(ref_pot3.values()), atol=1.0e-2)


def test__intmol():
    """ test pot.low_repulsion_struct
        test pot.intramol_interaction_potential_sum
    """

    assert automol.pot.low_repulsion_struct(
        PROP_GEO1, PROP_GEO2, thresh=40.0, potential='exp6')
    assert automol.pot.low_repulsion_struct(
        PROP_GEO1, PROP_GEO2, thresh=40.0, potential='lj_12_6')
