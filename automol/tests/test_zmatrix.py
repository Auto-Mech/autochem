""" test automol.zmatrix
"""
from automol import zmatrix
import automol


CH4O_ZMA = (
    (('H', (None, None, None), (None, None, None)),
     ('O', (0, None, None), ('R1', None, None)),
     ('C', (1, 0, None), ('R2', 'A2', None)),
     ('H', (2, 1, 0), ('R3', 'A3', 'D3')),
     ('H', (2, 1, 0), ('R3', 'A3', 'D4')),
     ('H', (2, 1, 0), ('R3', 'A3', 'D5'))),
    {'R1': 1.70075351,
     'R2': 2.64561657, 'A2': 1.74532925,
     'R3': 2.07869873, 'A3': 1.83259571,
     'D3': 1.04719755, 'D4': -1.04719755, 'D5': 3.1415926})

# CH4O_ZMA with R1 replaced with R7 and R2 with R8
CH4O_ZMA_2 = (
    (('H', (None, None, None), (None, None, None)),
     ('O', (0, None, None), ('R7', None, None)),
     ('C', (1, 0, None), ('R8', 'A2', None)),
     ('H', (2, 1, 0), ('R3', 'A3', 'D3')),
     ('H', (2, 1, 0), ('R3', 'A3', 'D4')),
     ('H', (2, 1, 0), ('R3', 'A3', 'D5'))),
    {'R7': 1.70075351,
     'R8': 2.64561657, 'A2': 1.74532925,
     'R3': 2.07869873, 'A3': 1.83259571,
     'D3': 1.04719755, 'D4': -1.04719755, 'D5': 3.1415926})

# CH4O_ZMA with R1 replaced with R7 and R2 with R8
CH4O_ZMA_3 = (
    (('H', (None, None, None), (None, None, None)),
     ('O', (0, None, None), ('R1', None, None)),
     ('C', (1, 0, None), ('R2', 'A2', None)),
     ('H', (0, 2, 4), ('R3', 'A3', 'D3')),
     ('H', (1, 2, 3), ('R3', 'A3', 'D4')),
     ('H', (2, 1, 0), ('R3', 'A3', 'D5'))),
    {'R1': 1.70075351,
     'R2': 2.64561657, 'A2': 1.74532925,
     'R3': 2.07869873, 'A3': 1.83259571,
     'D3': 1.04719755, 'D4': -1.04719755, 'D5': 3.1415926})

CH4O2_ZMA = (
    (('C', (None, None, None), (None, None, None)),
     ('O', (0, None, None), ('R1', None, None)),
     ('O', (0, 1, None), ('R2', 'A2', None)),
     ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
     ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
     ('H', (1, 0, 2), ('R5', 'A5', 'D5')),
     ('H', (2, 0, 1), ('R6', 'A6', 'D6'))),
    {'R1': 2.65933,
     'R2': 2.65933, 'A2': 1.90743,
     'R3': 2.06844, 'A3': 1.93366, 'D3': 4.1477,
     'R4': 2.06548, 'A4': 1.89469, 'D4': 2.06369,
     'R5': 1.83126, 'A5': 1.86751, 'D5': 1.44253,
     'R6': 1.83126, 'A6': 1.86751, 'D6': 4.84065})

C2O1H5_ZMA = ((('C', (None, None, None), (None, None, None)),
               ('C', (0, None, None), ('R1', None, None)),
               ('H', (0, 1, None), ('R2', 'A2', None)),
               ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
               ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
               ('O', (1, 0, 2), ('R5', 'A5', 'D5')),
               ('H', (1, 0, 5), ('R6', 'A6', 'D6')),
               ('H', (1, 0, 5), ('R7', 'A7', 'D7'))),
              {'R1': 2.87358,
               'R2': 2.09893, 'A2': 1.9205,
               'R3': 2.09835, 'A3': 1.9149, 'D3': 4.19162,
               'R4': 2.09893, 'A4': 1.9205, 'D4': 2.09994,
               'R5': 2.63652, 'A5': 1.9174, 'D5': 5.23321,
               'R6': 2.09831, 'A6': 1.9153, 'D6': 2.09638,
               'R7': 2.09831, 'A7': 1.9153, 'D7': 4.18680})
TS_ZMA = ((('C', (None, None, None), (None, None, None)),
           ('O', (0, None, None), ('R1', None, None)),
           ('H', (0, 1, None), ('R2', 'A2', None)),
           ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
           ('O', (1, 0, 2), ('R4', 'A4', 'D4')),
           ('H', (4, 1, 0), ('R5', 'A5', 'D5'))),
          {'R1': 2.27321, 'R2': 2.09434, 'A2': 2.1331937716650295,
           'R3': 2.09434, 'A3': 2.1331937716650295, 'D3': 3.141592653589793,
           'R5': 1.84616, 'A5': 1.4835298641951802, 'D4': 1.5707963267948966,
           'D5': 1.5707963267948966, 'R4': 3.0, 'A4': 1.4835298641951802})

CH4O2_ZMA_STR = """
C
O 1    R1
O 1    R2 2    A2
H 1    R3 2    A3 3    D3
H 1    R4 2    A4 3    D4
H 2    R5 1    A5 3    D5
H 3    R6 1    A6 2    D6

R1   =   1.407257
R2   =   1.407257
A2   = 109.287689
R3   =   1.094571
A3   = 110.790557
D3   = 237.645705
R4   =   1.093005
A4   = 108.557740
D4   = 118.240727
R5   =   0.969061
A5   = 107.000441
D5   =  82.650881
R6   =   0.969061
A6   = 107.000441
D6   = 277.348815
"""


def test__from_data():
    """ test getters
    """
    zma = automol.zmatrix.from_data(
        syms=zmatrix.symbols(CH4O2_ZMA),
        key_mat=zmatrix.key_matrix(CH4O2_ZMA),
        name_mat=zmatrix.name_matrix(CH4O2_ZMA),
        val_dct=zmatrix.values(CH4O2_ZMA),
    )
    assert zma == CH4O2_ZMA


def test__is_valid():
    """ test zmatrix.is_valid
    """
    zma1 = (
        (('H', (), ()),
         ('O', (0,), ('r1',)),
         ('C', (1, 0), ('r2', 'a1')),
         ('H', (2, 1, 0), ('r3', 'a2', 'd1')),
         ('H', (2, 1, 0), ('r3', 'a2', 'd2')),
         ('H', (2, 1, 0), ('r3', 'a2', 'd3'))),
        {'r1': 1.70075351,
         'r2': 2.64561657, 'a1': 1.74532925,
         'r3': 2.07869873, 'a2': 1.83259571,
         'd1': 1.04719755, 'd2': -1.04719755, 'd3': 3.1415926})

    zma2 = (
        (('H', (None, None, None), (None, None, None)),
         ('O', (0, None, None), ('r1', None, None)),
         ('C', (1, 0, None), ('r2', 'a1', None)),
         ('H', (2, 1, 0), ('r3', 'a2', 'd1')),
         ('H', (2, 1, 0), ('r3', 'a2', 'd2')),
         ('H', (2, 1, 0), ('r3', 'a2', 'd3'))),
        {'r1': 1.70075351,
         'r2': 2.64561657, 'a1': 1.74532925,
         'r3': 2.07869873, 'a2': 1.83259571,
         'd1': 1.04719755, 'd2': -1.04719755})

    assert zmatrix.is_valid(CH4O_ZMA)
    assert not zmatrix.is_valid(zma1)
    assert not zmatrix.is_valid(zma2)


def test__is_standard_form():
    """ test zmatrix.is_standard_form
    """
    nonst_zma = (
        (('H', (None, None, None), (None, None, None)),
         ('O', (0, None, None), ('r1', None, None)),
         ('C', (1, 0, None), ('r2', 'a1', None)),
         ('H', (2, 1, 0), ('r3', 'a2', 'd1')),
         ('H', (2, 1, 0), ('r3', 'a2', 'd2')),
         ('H', (2, 1, 0), ('r3', 'a2', 'd3'))),
        {'r1': 1.70075351,
         'r2': 2.64561657, 'a1': 1.74532925,
         'r3': 2.07869873, 'a2': 1.83259571,
         'd1': 1.04719755, 'd2': -1.04719755, 'd3': 3.1415926})
    zma = zmatrix.standard_form(nonst_zma)

    assert zmatrix.is_valid(nonst_zma)
    assert not zmatrix.is_standard_form(nonst_zma)
    assert zmatrix.is_standard_form(zma)
    assert zmatrix.almost_equal(zma, CH4O_ZMA)


def test__coordinate_key_matrix():
    """ test zmatrix.coordinate_key_matrix
    """
    coo_key_mat = ((None, None, None),
                   ((1, 0), None, None),
                   ((2, 1), (2, 1, 0), None),
                   ((3, 2), (3, 2, 1), (3, 2, 1, 0)),
                   ((4, 2), (4, 2, 1), (4, 2, 1, 0)),
                   ((5, 2), (5, 2, 1), (5, 2, 1, 0)))

    assert (zmatrix.coordinate_key_matrix(CH4O_ZMA)
            == coo_key_mat)

    one_indexed_coo_key_mat = ((None, None, None),
                               ((2, 1), None, None),
                               ((3, 2), (3, 2, 1), None),
                               ((4, 3), (4, 3, 2), (4, 3, 2, 1)),
                               ((5, 3), (5, 3, 2), (5, 3, 2, 1)),
                               ((6, 3), (6, 3, 2), (6, 3, 2, 1)))

    assert (zmatrix.coordinate_key_matrix(CH4O_ZMA, shift=1)
            == one_indexed_coo_key_mat)


def test__coordinates():
    """ test zmatrix.coordinates
    """
    coo_dct = {
        'R1': ((1, 0),),
        'R2': ((2, 1),), 'A2': ((2, 1, 0),),
        'R3': ((3, 2), (4, 2), (5, 2)),
        'A3': ((3, 2, 1), (4, 2, 1), (5, 2, 1)),
        'D3': ((3, 2, 1, 0),), 'D4': ((4, 2, 1, 0),), 'D5': ((5, 2, 1, 0),)}
    assert zmatrix.coordinates(CH4O_ZMA) == coo_dct


def test__bond_idxs():
    """ test bond indices
    """
    idxs = (2, 1)
    assert zmatrix.bond_idxs(CH4O_ZMA, 'R2') == idxs


def test__get_babs1():
    """ test get_babs1
    """
    dist_name = 'R4'
    assert zmatrix.get_babs1(TS_ZMA, dist_name) == 'D4'


def test__distance_names():
    """ test zmatrix.distance_names
    """
    assert zmatrix.distance_names(CH4O_ZMA) == ('R1', 'R2', 'R3')


def test__central_angle_names():
    """ test zmatrix.angle_names
    """
    assert zmatrix.central_angle_names(CH4O_ZMA) == ('A2', 'A3')


def test__dihedral_angle_names():
    """ test zmatrix.dihedral_names
    """
    assert zmatrix.dihedral_angle_names(CH4O_ZMA) == ('D3', 'D4', 'D5')


def test__angle_names():
    """ test zmatrix.angle_names
    """
    assert zmatrix.angle_names(CH4O_ZMA) == ('A2', 'A3', 'D3', 'D4', 'D5')


def test__from_string():
    """ test zmatrix.from_vmat_string
    """
    zma = zmatrix.from_string(CH4O2_ZMA_STR)
    assert zmatrix.almost_equal(zma, CH4O2_ZMA)


def test__string():
    """ test zmatrix.vmat_string
    """
    zma = zmatrix.from_string(zmatrix.string(CH4O_ZMA))
    assert zmatrix.almost_equal(zma, CH4O_ZMA)


def test__set_values():
    """ test zmatrix.set_values
    """
    val = zmatrix.values(CH4O_ZMA)['D3']

    zma = zmatrix.set_values(CH4O_ZMA, {'D3': val + 1e-6})
    assert zmatrix.almost_equal(zma, CH4O_ZMA)

    zma = zmatrix.set_values(CH4O_ZMA, {'D3': val + 1e-1})
    assert not zmatrix.almost_equal(zma, CH4O_ZMA)


def test__set_names():
    """ test zmatrix.set_names
    """
    name_dct = {'R1': 'R7',
                'R2': 'R8'}
    zma = zmatrix.set_names(CH4O_ZMA, name_dct)
    assert zmatrix.almost_equal(zma, CH4O_ZMA_2)


def test__set_keys():
    """ test zmatrix.set_keys
    """
    key_dct = {3: (0, 2, 4),
               4: (1, 2, 3)}
    zma = zmatrix.set_keys(CH4O_ZMA, key_dct)
    assert zmatrix.almost_equal(zma, CH4O_ZMA_3)


def test__join():
    """ test zmatrix.join
    """
    zma1 = ((('C', (None, None, None), (None, None, None)),
             ('O', (0, None, None), ('R1', None, None)),
             ('H', (0, 1, None), ('R2', 'A2', None)),
             ('H', (0, 1, 2), ('R3', 'A3', 'D3'))),
            {'R1': 2.31422,
             'R2': 2.08191, 'A2': 2.13342,
             'R3': 2.08191, 'A3': 2.13342, 'D3': 3.14159})
    zma2 = ((('X', (None, None, None), (None, None, None)),), {})
    zma3 = ((('N', (None, None, None), (None, None, None)),
             ('O', (0, None, None), ('R1', None, None)),
             ('O', (0, 1, None), ('R2', 'A2', None))),
            {'R1': 2.69082, 'R2': 2.69082, 'A2': 1.89019})

    zma1 = zmatrix.standard_form(zma1)
    zma3 = zmatrix.standard_form(zma3, shift=zmatrix.count(zma1)+1)

    join_key_mat1 = ((3, 0, 1),)
    join_name_mat1 = (('JR4', 'JA4', 'JD4'),)
    join_val_dct1 = {'JR4': 3.78, 'JA4': 1.57, 'JD4': 3.14}
    zma = zmatrix.join(zma1, zma2, join_key_mat1, join_name_mat1,
                       join_val_dct1)

    join_key_mat2 = ((3, 4, 0),
                     (None, 3, 4),
                     (None, None, 3))
    join_name_mat2 = (('JR5', 'JA5', 'JD5'),
                      (None, 'JA6', 'JD6'),
                      (None, None, 'JD7'))
    val_dct2 = {'JR5': 2.38, 'JA5': 1.53, 'JD5': 3.32,
                'JA6': 1.83, 'JD6': 1.39, 'JD7': 0.}
    zma = zmatrix.join(zma, zma3, join_key_mat2, join_name_mat2, val_dct2)
    assert zma == (
        (('C', (None, None, None), (None, None, None)),
         ('O', (0, None, None), ('R1', None, None)),
         ('H', (0, 1, None), ('R2', 'A2', None)),
         ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
         ('X', (3, 0, 1), ('JR4', 'JA4', 'JD4')),
         ('N', (3, 4, 0), ('JR5', 'JA5', 'JD5')),
         ('O', (5, 3, 4), ('R6', 'JA6', 'JD6')),
         ('O', (5, 6, 3), ('R7', 'A7', 'JD7'))),
        {'R1': 2.31422,
         'R2': 2.08191, 'A2': 2.13342,
         'R3': 2.08191, 'A3': 2.13342, 'D3': 3.14159,
         'JR4': 3.78, 'JA4': 1.57, 'JD4': 3.14,
         'JR5': 2.38, 'JA5': 1.53, 'JD5': 3.32,
         'R6': 2.69082, 'JA6': 1.83, 'JD6': 1.39,
         'R7': 2.69082, 'A7': 1.89019, 'JD7': 0.0})


def test__torsional_symmetry_numbers():
    """ test zmatrix.torsional_symmetry_numbers
    """
    assert zmatrix.torsional_symmetry_numbers(CH4O_ZMA, ('D3',)) == (3,)
    assert zmatrix.torsional_symmetry_numbers(CH4O2_ZMA, ('D5', 'D6')) == (
        1, 1)


def test__samples():
    """ test zmatrix.samples
    """
    tors_names = ['D3']
    tors_range_vals = zmatrix.torsional_sampling_ranges(CH4O_ZMA, tors_names)
    tors_range_dct = dict(zip(tors_names, tors_range_vals))
    zmas = zmatrix.samples(CH4O_ZMA, 9, tors_range_dct)
    assert len(zmas) == 9

    tors_names = ['D5', 'D6']
    tors_range_vals = zmatrix.torsional_sampling_ranges(CH4O2_ZMA, tors_names)
    tors_range_dct = dict(zip(tors_names, tors_range_vals))
    zmas = zmatrix.samples(CH4O2_ZMA, 7, tors_range_dct)
    assert len(zmas) == 7


def test__ts__addition():
    """ test zmatrix.ts.addition
    """
    rct_zmas = [
        ((('O', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None))),
         {'R1': 2.15608}),
        ((('H', (None, None, None), (None, None, None)),),
         {}),
    ]
    prd_zmas = [
        ((('O', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None))),
         {'R1': 2.48959, 'R2': 1.86213, 'A2': 1.9084302705931997})
    ]
    ts_zma, dist_name, tors_names = zmatrix.ts.addition(rct_zmas, prd_zmas)
    assert dist_name == 'R2'
    assert not tors_names
    assert zmatrix.almost_equal(
        ts_zma,
        ((('O', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None)),
          ('H', (1, 0, None), ('R2', 'A2', None))),
         {'R1': 2.15608, 'A2': 1.4835298641951802, 'R2': 3.0})
    )

    rct_zmas = list(reversed(rct_zmas))
    ts_zma, dist_name, tors_names = zmatrix.ts.addition(rct_zmas, prd_zmas)
    assert dist_name == 'R1'
    assert not tors_names
    assert zmatrix.almost_equal(
        ts_zma,
        ((('H', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None)),
          ('O', (1, 0, None), ('R2', 'A2', None))),
         {'R2': 2.15608, 'R1': 3.0, 'A2': 1.4835298641951802})
    )

    rct_zmas = [
        ((('O', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('x1', None, None)),
          ('H', (0, 1, None), ('x2', 'y2', None))),
         {'x1': 2.48959, 'x2': 1.86213, 'y2': 1.9084302705931997}),
        ((('C', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('x1', None, None)),
          ('H', (0, 1, None), ('x2', 'y2', None)),
          ('H', (0, 1, 2), ('x3', 'y3', 'z3'))),
         {'x1': 2.045,
          'x2': 2.045, 'y2': 2.0943,
          'x3': 2.045, 'y3': 2.0943, 'z3': 3.1415}),
    ]
    prd_zmas = [
        ((('C', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
          ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
          ('O', (1, 0, 2), ('R5', 'A5', 'D5')),
          ('H', (5, 1, 0), ('R6', 'A6', 'D6'))),
         {'R1': 2.687,
          'R2': 2.065, 'A2': 1.894,
          'R3': 2.067, 'A3': 1.920, 'D3': 4.202,
          'R4': 2.067, 'A4': 1.920, 'D4': 2.080,
          'R5': 2.750, 'A5': 1.842, 'D5': 3.140,
          'R6': 1.840, 'A6': 1.680, 'D6': 2.055})
    ]
    ts_zma, dist_name, tors_names = zmatrix.ts.addition(rct_zmas, prd_zmas)
    assert dist_name == 'R3'
    assert tors_names == ('D4',)
    assert zmatrix.almost_equal(
        ts_zma,
        ((('O', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('C', (1, 0, 2), ('R3', 'A3', 'D3')),
          ('H', (3, 1, 0), ('R4', 'A4', 'D4')),
          ('H', (3, 4, 1), ('R5', 'A5', 'D5')),
          ('H', (3, 4, 5), ('R6', 'A6', 'D6'))),
         {'R1': 2.48959,
          'R2': 1.86213, 'A2': 1.9084302705931997,
          'R4': 2.045, 'R5': 2.045,
          'A5': 2.0943, 'R6': 2.045,
          'A6': 2.0943, 'D6': 3.1415,
          'D3': 1.4835298641951802,
          'R3': 3.0,
          'D5': 1.4835298641951802,
          'D4': 1.4835298641951802,
          'A3': 1.4835298641951802, 'A4': 1.4835298641951802})
    )


def test__ts__hydrogen_abstraction():
    """ test zmatrix.ts.hydrogen_abstraction
    """
    rct_zmas = [
        ((('C', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
          ('H', (0, 1, 2), ('R4', 'A4', 'D4'))),
         {'R1': 2.063,
          'R2': 2.063, 'A2': 1.9106,
          'R3': 2.063, 'A3': 1.9106, 'D3': 2.0943,
          'R4': 2.063, 'A4': 1.9106, 'D4': 4.1887}),
        ((('H', (None, None, None), (None, None, None)),),
         {}),
    ]
    prd_zmas = [
        ((('C', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3'))),
         {'R1': 2.045,
          'R2': 2.045, 'A2': 2.0943,
          'R3': 2.045, 'A3': 2.0943, 'D3': 3.1415}),
        ((('H', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('R1', None, None))),
         {'R1': 1.31906}),
    ]
    ts_zma, dist_name, frm_key, brk_key, tors_names = (
        zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas))
    assert dist_name == 'R6'
    assert frm_key == frozenset({4, 5})
    assert brk_key == frozenset({0, 4})
    assert not tors_names
    assert zmatrix.almost_equal(
        ts_zma,
        ((('C', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
          ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
          ('X', (4, 0, 1), ('R5', 'A5', 'D5')),
          ('H', (4, 5, 0), ('R6', 'A6', 'D6'))),
         {'R1': 2.063,
          'R2': 2.063, 'A2': 1.9106,
          'R3': 2.063, 'A3': 1.9106, 'D3': 2.0943,
          'R4': 2.063, 'A4': 1.9106, 'D4': 4.1887,
          'A5': 1.5707963267948966, 'D5': 3.141592653589793,
          'R5': 1.8897261254578281, 'A6': 1.4835298641951802,
          'R6': 3.0, 'D6': 3.0543261909900767})
    )

    # example 2
    rct_zmas = [
        ((('O', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None))),
         {'R1': 1.83114, 'R2': 1.83115, 'A2': 1.814758449638664}),
        ((('C', (None, None, None), (None, None, None)),
          ('X', (0, None, None), ('R1', None, None)),
          ('C', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
          ('H', (2, 0, 1), ('R4', 'A4', 'D4')),
          ('H', (2, 0, 4), ('R5', 'A5', 'D5'))),
         {'R1': 1,
          'R2': 2.45333, 'A2': 1.570796,
          'R3': 2.01549, 'A3': 1.570796, 'D3': 3.14159,
          'R4': 2.04888, 'A4': 2.106804, 'D4': 0.00000,
          'R5': 2.04888, 'A5': 2.106804, 'D5': 3.14159}),
    ]
    prd_zmas = [
        ((('O', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('R1', None, None))),
         {'R1': 1.84779}),
        ((('C', (None, None, None), (None, None, None)),
          ('C', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
          ('H', (1, 0, 2), ('R4', 'A4', 'D4')),
          ('H', (1, 0, 4), ('R5', 'A5', 'D5'))),
         {'R1': 2.52427,
          'R2': 2.05128, 'A2': 2.112721059539136,
          'R3': 2.05128, 'A3': 2.112721059539136, 'D3': 3.141592653589793,
          'R4': 2.05128, 'A4': 2.112721059539136, 'D4': 3.141592653589793,
          'R5': 2.05128, 'A5': 2.112721059539136, 'D5': 3.141592653589793})
    ]
    ts_zma, dist_name, frm_key, brk_key, tors_names = (
        zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas))
    assert dist_name == 'R4'
    assert frm_key == frozenset({2, 3})
    assert brk_key == frozenset({0, 2})
    assert tors_names == ('D5',)
    assert zmatrix.almost_equal(
        ts_zma,
        ((('O', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('X', (2, 0, 1), ('R3', 'A3', 'D3')),
          ('C', (2, 3, 0), ('R4', 'A4', 'D4')),
          ('X', (4, 2, 3), ('R5', 'A5', 'D5')),
          ('C', (4, 5, 2), ('R6', 'A6', 'D6')),
          ('H', (4, 5, 6), ('R7', 'A7', 'D7')),
          ('H', (6, 4, 5), ('R8', 'A8', 'D8')),
          ('H', (6, 4, 8), ('R9', 'A9', 'D9'))),
         {'R1': 1.83114,
          'R2': 1.83115, 'A2': 1.814758449638664,
          'A3': 1.5707963267948966, 'D3': 3.141592653589793,
          'R3': 1.8897261254578281, 'R5': 1,
          'R6': 2.45333, 'A6': 1.570796,
          'R7': 2.01549, 'A7': 1.570796, 'D7': 3.14159,
          'R8': 2.04888, 'A8': 2.106804, 'D8': 0.0,
          'R9': 2.04888, 'A9': 2.106804, 'D9': 3.14159,
          'R4': 3.0, 'D6': 1.4835298641951802, 'D4': 3.0543261909900767,
          'A4': 1.4835298641951802, 'A5': 1.4835298641951802,
          'D5': 1.4835298641951802})
    )


def test__ts__hydrogen_migration():
    """ test zmatrix.ts.hydrogen_migration
    """

    # hydrogen migration
    print('hydrogen migration')
    rct_zmas1 = [
        ((('C', (None, None, None), (None, None, None)),
          ('C', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
          ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
          ('O', (1, 0, 2), ('R5', 'A5', 'D5')),
          ('H', (1, 0, 5), ('R6', 'A6', 'D6')),
          ('H', (1, 0, 5), ('R7', 'A7', 'D7'))),
         {'R1': 2.87358,
          'R2': 2.09893, 'A2': 1.9205,
          'R3': 2.09835, 'A3': 1.9149, 'D3': 4.19162,
          'R4': 2.09893, 'A4': 1.9205, 'D4': 2.09994,
          'R5': 2.63652, 'A5': 1.9174, 'D5': 5.23321,
          'R6': 2.09831, 'A6': 1.9153, 'D6': 2.09638,
          'R7': 2.09831, 'A7': 1.9153, 'D7': 4.18680})
    ]

    prd_zmas1 = [
        ((('C', (None, None, None), (None, None, None)),
          ('C', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
          ('O', (1, 0, 2), ('R4', 'A4', 'D4')),
          ('H', (1, 0, 4), ('R5', 'A5', 'D5')),
          ('H', (1, 0, 4), ('R6', 'A6', 'D6')),
          ('H', (4, 1, 0), ('R7', 'A7', 'D7'))),
         {'R1': 2.86956,
          'R2': 2.0985, 'A2': 1.9202,
          'R3': 2.09904, 'A3': 1.9244, 'D3': 4.1802,
          'R4': 2.64409, 'A4': 1.906, 'D4': 3.15423,
          'R5': 2.10184, 'A5': 1.8926, 'D5': 2.0903,
          'R6': 2.10265, 'A6': 1.9011, 'D6': 4.1848,
          'R7': 1.87536, 'A7': 1.8507, 'D7': 3.1216})
    ]

    ts_zma, dist_name, frm_key, brk_key, tors_names = (
        zmatrix.ts.hydrogen_migration(rct_zmas1, prd_zmas1))
    print(zmatrix.string(ts_zma))
    print(dist_name)
    print(tors_names)

    print('\nproton migration')

    rct_zmas2 = [
        ((('N', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None)),
          ('O', (0, 1, None), ('R2', 'A2', None)),
          ('H', (1, 0, 2), ('R3', 'A3', 'D3'))),
         {'R1': 2.48392,
          'R2': 2.25378, 'A2': 2.0944649155632753,
          'R3': 1.82641, 'A3': 2.0947267149510744, 'D3': 3.141575200297273})
    ]

    prd_zmas2 = [
        ((('N', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None)),
          ('O', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3'))),
         {'R1': 2.4839,
          'R2': 2.25363, 'A2': 2.0943951023931953,
          'R3': 1.94544, 'A3': 2.0943252892231157, 'D3': 3.1514363105710412})
    ]

    ts_zma, dist_name, frm_key, brk_key, tors_names = (
        zmatrix.ts.hydrogen_migration(rct_zmas2, prd_zmas2))
    print(zmatrix.string(ts_zma))
    print(dist_name)
    print(tors_names)


def test__ts__elimination():
    """ test zmatrix.ts.elimination
    """
    rct_zmas = [
        ((('C', (None, None, None), (None, None, None)),
          ('C', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
          ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
          ('O', (1, 0, 2), ('R5', 'A5', 'D5')),
          ('H', (1, 0, 5), ('R6', 'A6', 'D6')),
          ('H', (1, 0, 5), ('R7', 'A7', 'D7')),
          ('O', (5, 1, 0), ('R8', 'A8', 'D8'))),
         {'R1': 2.88178,
          'R2': 2.09836, 'A2': 1.9249236320245462,
          'R3': 2.09877, 'A3': 1.9198796304862826, 'D3': 4.181093302785096,
          'R4': 2.09829, 'A4': 1.914207310417301, 'D4': 2.0904157516986484,
          'R5': 2.65619, 'A5': 1.9326205340258409, 'D5': 1.0871254484944721,
          'R6': 2.10699, 'A6': 1.933894624379797, 'D6': 4.1123273302565195,
          'R7': 2.09849, 'A7': 1.8998083440883475, 'D7': 2.0737304040495825,
          'R8': 2.49456, 'A8': 1.858164788135763, 'D8': 4.852695998952513})
    ]

    prd_zmas = [
        ((('C', (None, None, None), (None, None, None)),
          ('C', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
          ('H', (1, 0, 2), ('R4', 'A4', 'D4')),
          ('H', (1, 0, 4), ('R5', 'A5', 'D5'))),
         {'R1': 2.51566,
          'R2': 2.05124, 'A2': 2.098950411740901,
          'R3': 2.05124, 'A3': 2.098950411740901, 'D3': 3.141592653589793,
          'R4': 2.05124, 'A4': 2.098950411740901, 'D4': 6.283185307179586,
          'R5': 2.05124, 'A5': 2.098950411740901, 'D5': 3.141592653589793}),
        ((('O', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None)),
          ('H', (1, 0, None), ('R2', 'A2', None))),
         {'R1': 2.48688,
          'R2': 1.87131, 'A2': 1.824043601259274})
    ]

    ts_zma, dist_name, dist_name2, tors_names = (
        zmatrix.ts.concerted_unimolecular_elimination(rct_zmas, prd_zmas))
    print(zmatrix.string(ts_zma))
    print(dist_name)
    print(dist_name2)
    print(tors_names)
    print(automol.geom.string(automol.zmatrix.geometry(ts_zma)))


if __name__ == '__main__':
    # test__from_data()
    # test__string()
    # test__set_keys()
    # test__join()
    # test__ts__hydrogen_abstraction()
    # test__is_standard_form()
    # test__join()
    # test__from_string()
    # test__ts__addition()
    test__ts__hydrogen_abstraction()
    # test__ts__hydrogen_migration()
    # test__ts__elimination()
    # test__bond_idxs()
    # test__get_babs1()
