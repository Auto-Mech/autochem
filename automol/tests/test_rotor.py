"""
    test automol.rotor
"""

import numpy
import automol


SPC_ZMA = (
    (('C', (None, None, None), (None, None, None)),
     ('C', (0, None, None), ('R1', None, None)),
     ('H', (0, 1, None), ('R2', 'A2', None)),
     ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
     ('C', (1, 0, 2), ('R4', 'A4', 'D4')),
     ('O', (1, 0, 4), ('R5', 'A5', 'D5')),
     ('H', (1, 0, 4), ('R6', 'A6', 'D6')),
     ('H', (4, 1, 0), ('R7', 'A7', 'D7')),
     ('H', (4, 1, 7), ('R8', 'A8', 'D8')),
     ('H', (4, 1, 7), ('R9', 'A9', 'D9')),
     ('O', (5, 1, 0), ('R10', 'A10', 'D10')),
     ('H', (10, 5, 1), ('R11', 'A11', 'D11'))),
    {'R1': 2.81758165305762,
     'R2': 2.04865209260883, 'A2': 2.12613042418220,
     'R3': 2.05413229837265, 'A3': 2.08003802396629, 'D3': 3.03605084839244,
     'R4': 2.89052508150029, 'A4': 1.96723215842214, 'D4': 4.37781983475288,
     'R5': 2.70646575688070, 'A5': 1.94415017906451, 'D5': 4.24894995877338,
     'R6': 2.07227366917705, 'A6': 1.92579804197979, 'D6': 2.15147085959191,
     'R7': 2.06717140863831, 'A7': 1.92170873554237, 'D7': 3.17812064950478,
     'R8': 2.06641551818813, 'A8': 1.91917102680997, 'D8': 4.18693317446226,
     'R9': 2.06792729908850, 'A9': 1.92728855316099, 'D9': 2.09678271280992,
     'R10': 2.6947494549028, 'A10': 1.8812781834199, 'D10': 5.0954346459906,
     'R11': 1.8324674238564, 'A11': 1.7685054791315, 'D11': 1.6560504247670})
SPC_TORS_NAMES = ('D4', 'D7', 'D10', 'D11')
SPC_FRM_KEY = frozenset({})
SPC_BRK_KEY = frozenset({})

SPC2_ZMA = (
    (('C', (None, None, None), (None, None, None)),
     ('C', (0, None, None), ('R1', None, None)),
     ('H', (0, 1, None), ('R2', 'A2', None)),
     ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
     ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
     ('H', (1, 0, 2), ('R5', 'A5', 'D5')),
     ('H', (1, 0, 5), ('R6', 'A6', 'D6')),
     ('H', (1, 0, 5), ('R7', 'A7', 'D7'))),
    {'R1': 2.85738,
     'R2': 2.06752, 'A2': 1.9297756473450902,
     'R3': 2.06752, 'A3': 1.9297756473450902, 'D3': 2.0943951023931953,
     'R4': 2.06752, 'A4': 1.9297756473450902, 'D4': 4.1887902047863905,
     'R5': 2.06752, 'A5': 1.9297756473450902, 'D5': 5.235987755982989,
     'R6': 2.06752, 'A6': 1.9297756473450902, 'D6': 4.1887902047863905,
     'R7': 2.06752, 'A7': 1.9297756473450902, 'D7': 2.0943951023931953})
SPC2_NAME = 'D5'
SPC2_TRANS = ('', frozenset({}), frozenset({}))

TS_ZMA = (
    (('C', (None, None, None), (None, None, None)),
     ('C', (0, None, None), ('R1', None, None)),
     ('H', (0, 1, None), ('R2', 'A2', None)),
     ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
     ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
     ('O', (1, 0, 2), ('R5', 'A5', 'D5')),
     ('H', (1, 0, 5), ('R6', 'A6', 'D6')),
     ('H', (1, 0, 5), ('R7', 'A7', 'D7')),
     ('H', (5, 1, 0), ('R8', 'A8', 'D8')),
     ('X', (7, 1, 0), ('R9', 'A9', 'D9')),
     ('C', (7, 9, 1), ('R10', 'A10', 'D10')),
     ('H', (10, 7, 9), ('R11', 'A11', 'D11')),
     ('H', (10, 11, 7), ('R12', 'A12', 'D12')),
     ('H', (10, 11, 12), ('R13', 'A13', 'D13'))),
    {'R1': 2.86047843610551,
     'R2': 2.06868318953868, 'A2': 1.92742643417190,
     'R3': 2.07189572395196, 'A3': 1.94002247538354, 'D3': 4.17755377506205,
     'R4': 2.07454134052760, 'A4': 1.93385797246550, 'D4': 2.07811118047208,
     'R5': 2.63465616411330, 'A5': 2.00629786306952, 'D5': 3.13352050579931,
     'R6': 2.06509270990031, 'A6': 1.97849825874376, 'D6': 4.11924756074067,
     'R7': 2.41053464563400, 'A7': 1.85996596792382, 'D7': 2.11338428465489,
     'R8': 1.82415262890444, 'A8': 1.88251213120108, 'D8': 0.89541150213840,
     'R11': 2.0599904493615, 'R12': 2.0581007232361, 'A12': 1.9886403670271,
     'R13': 2.0584786684612, 'A13': 1.9847552641121, 'D13': 2.3199789082134,
     'D12': 1.9638130584174, 'D10': 3.2229703752932, 'A11': 1.8333968207207,
     'A10': 1.5739204661559, 'D11': 2.3565871892740, 'R10': 2.7297093882238,
     'D9': 3.14159265358979, 'A9': 1.57079632679489, 'R9': 1.88972612545782})
TS_TORS_NAMES = ('D5', 'D8', 'D11')
TS_TRANS = ('', frozenset({7, 10}), frozenset({1, 7}))


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


TORS_DCT = {
    'D12': {
        'group1': (5, 6, 7),
        'group2': (1, 2),
        'axis': (3, 4),
        'span': 6.28,
        'sym_num': 1,
    },
    'D10': {
        'group1': (1, 2, 3, 4, 5, 6, 7, 8),
        'group2': (11, 12),
        'axis': (9, 10),
        'span': 3.14,
        'sym_num': 2
    }
}


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


def test__tors_inf():
    """ test automol.rotor.tors.symmetry_number
        test.automol.rotor.tors.grid
    """

    ref_sym_num1 = (1, 3, 1, 1)
    ref_sym_num2 = 3
    ref_sym_num3 = (3, 1, 1)

    sym_num1 = tuple(
        automol.rotor.tors.symmetry_number(
            SPC_ZMA, name,
            frm_bnd_keys=SPC_FRM_KEY, brk_bnd_keys=SPC_BRK_KEY)
        for name in SPC_TORS_NAMES)
    sym_num2 = automol.rotor.tors.symmetry_number(
        SPC2_ZMA, SPC2_NAME,
        frm_bnd_keys=SPC_FRM_KEY, brk_bnd_keys=SPC_BRK_KEY)
    sym_num3 = tuple(
        automol.rotor.tors.symmetry_number(
            TS_ZMA, name,
            frm_bnd_keys=SPC_FRM_KEY, brk_bnd_keys=SPC_BRK_KEY)
        for name in TS_TORS_NAMES)

    assert numpy.allclose(sym_num1, ref_sym_num1)
    assert numpy.isclose(sym_num2, ref_sym_num2)
    assert numpy.allclose(sym_num3, ref_sym_num3)

    ref_span1 = (6.283185307179586, 2.0943951023931953,
                 6.283185307179586, 6.283185307179586)
    ref_span2 = 2.0943951023931953
    ref_span3 = (2.0943951023931953, 6.283185307179586, 6.283185307179586)

    span1 = tuple(
        automol.rotor.tors.span(
            SPC_ZMA, name,
            frm_bnd_keys=SPC_FRM_KEY, brk_bnd_keys=SPC_BRK_KEY)
        for name in SPC_TORS_NAMES)
    span2 = automol.rotor.tors.span(
        SPC2_ZMA, SPC2_NAME,
        frm_bnd_keys=SPC_FRM_KEY, brk_bnd_keys=SPC_BRK_KEY)
    span3 = tuple(
        automol.rotor.tors.span(
            TS_ZMA, name,
            frm_bnd_keys=SPC_FRM_KEY, brk_bnd_keys=SPC_BRK_KEY)
        for name in TS_TORS_NAMES)

    assert numpy.allclose(span1, ref_span1)
    assert numpy.isclose(span2, ref_span2)
    assert numpy.allclose(span3, ref_span3)

    ref_grid1 = (
        (4.377819834752880, 4.9014185899510245, 5.425017345149168,
         5.948616100347312, 6.4722148555454570, 6.9958136107436015,
         7.519412365941745, 8.0430111211398900, 8.566609876338035,
         9.090208631536179, 9.6138073867343220, 10.137406141932466),
        (3.178120649504780, 3.70171935030251160, 4.225318051100244,
         4.748916751897975),
        (5.095434645990600, 5.6190334011887440, 6.1426321563868880,
         6.666230911585032, 7.1898296667831770, 7.7134284219813210,
         8.237027177179465, 8.7606259323776090, 9.2842246875757550,
         9.807823442773898, 10.331422197972042, 10.855020953170186),
        (1.656050424767000, 2.1796491799651445, 2.7032479351632883,
         3.226846690361433, 3.7504454455595770, 4.2740442007577215,
         4.797642955955865, 5.3212417111540100, 5.8448404663521540,
         6.368439221550298, 6.8920379767484420, 7.4156367319465860)
    )
    ref_grid2 = (
        (3.13352050579931, 3.6571192065970415, 4.180717907394773,
         4.704316608192505),
        (0.8954115021384, 1.4190102573365442, 1.9426090125346884,
         2.4662077677328327, 2.989806522930977, 3.513405278129121,
         4.037004033327266, 4.56060278852541, 5.084201543723554,
         5.607800298921698, 6.131399054119843, 6.654997809317987),
        (2.356587189274, 2.880185944472144, 3.4037846996702883,
         3.9273834548684325, 4.450982210066577, 4.974580965264721,
         5.498179720462865, 6.0217784756610095, 6.545377230859154,
         7.068975986057298, 7.592574741255442, 8.116173496453586)
    )

    grid1 = tuple(
        automol.rotor.tors.grid(
            SPC_ZMA, name, SCAN_INCREMENT,
            frm_bnd_keys=SPC_FRM_KEY, brk_bnd_keys=SPC_BRK_KEY)
        for name in SPC_TORS_NAMES)
    grid2 = tuple(
        automol.rotor.tors.grid(
            TS_ZMA, name, SCAN_INCREMENT,
            frm_bnd_keys=SPC_FRM_KEY, brk_bnd_keys=SPC_BRK_KEY)
        for name in TS_TORS_NAMES)

    for row, ref_row in zip(grid1, ref_grid1):
        assert numpy.allclose(row, ref_row)
    for row, ref_row in zip(grid2, ref_grid2):
        assert numpy.allclose(row, ref_row)


def test__tors_idx():
    """ test automol.rotor.tors.axis
        test.automol.rotor.tors.groups
    """

    axis = automol.rotor.tors.torsional_axis(zma, tors_name)
    group = automol.rotor.tors.torsional_groups()


def test__tors_string():
    """ test.automol.rotor.tors.string
        test.automol.rotor.tors.from_string
    """

    tors_str = automol.rotor.tors.string(TORS_DCT)
    print(tors_str)
    tors_dct = automol.rotor.tors.from_string(tors_str)
    print(tors_dct)


def test__tors_build():
    """ test automol.rotor.tors
    """

    sym_num = automol.rotor.tors.symmetry_number(TORS_DCT)
    grid = automol.rotor.tors.grid(TORS_DCT)
    span = automol.rotor.tors.span(TORS_DCT)
    axis = automol.rotor.tors.axis(TORS_DCT)
    groups = automol.rotor.tors.rotational_groups(TORS_DCT)
    geo = automol.rotor.tors.geometry(TORS_DCT)
    zma = automol.rotor.tors.zmatrix(TORS_DCT)

if __name__ == '__main__':
    test__valid_potential()
    test__build_potential()
    test__transform_potential()
    test__tors_inf()
    test__tors_string()

