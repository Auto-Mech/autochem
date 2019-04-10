""" test the automol.zmatrix module
"""
import automol
from automol import zmatrix

CH4O_ZMA = (
    (('H', (None, None, None), (None, None, None)),
     ('O', (0, None, None), ('d1', None, None)),
     ('C', (1, 0, None), ('d2', 'a1', None)),
     ('H', (2, 1, 0), ('d3', 'a2', 't1')),
     ('H', (2, 1, 0), ('d3', 'a2', 't2')),
     ('H', (2, 1, 0), ('d3', 'a2', 't3'))),
    {'d1': 1.70075351, 'd2': 2.64561657, 'a1': 1.74532925,
     'd3': 2.07869873, 'a2': 1.83259571, 't1': 1.04719755,
     't2': -1.04719755, 't3': 3.1415926})

CH4O_GRA = (
    {0: ('H', 0, None), 1: ('O', 0, None), 2: ('C', 0, None),
     3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
     frozenset({2, 3}): (1, None), frozenset({2, 4}): (1, None),
     frozenset({2, 5}): (1, None)})

CH4O_ZMA_STR = """
H
O  1    d1
C  2    d2 1    a1
H  3    d3 2    a2 1    t1
H  3    d3 2    a2 1    t2
H  3    d3 2    a2 1    t3

d1    =    0.900000
d2    =    1.400000
a1    =  100.000000
d3    =    1.100000
a2    =  105.000000
t1    =   60.000000
t2    =  -60.000000
t3    =  179.999997
"""

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
     'R3': 2.06844, 'A3': 1.93366, 'D3': 4.14770,
     'R4': 2.06548, 'A4': 1.89469, 'D4': 2.06369,
     'R5': 1.83126, 'A5': 1.86751, 'D5': 1.44253,
     'R6': 1.83126, 'A6': 1.86751, 'D6': 4.84065})

H2O_ZMA_STR = """Geometry (in Angstrom), charge = 0, multiplicity = 1:

    O
    H             1          R1
    H             1          R2      2          A2

    A2        =   99.9891574659
    R1        =    0.9894851960
    R2        =    0.9894852414
"""

H2O_ZMA = (
    (('O', (None, None, None), (None, None, None)),
     ('H', (0, None, None), ('R1', None, None)),
     ('H', (0, 1, None), ('R2', 'A2', None))),
    {'A2': 1.7451400140750248, 'R1': 1.8698560256349597,
     'R2': 1.869856111428526})

CH4O_VAR_ZMA = (
    (('H', (None, None, None), (None, None, None)),
     ('O', (0, None, None), ('d1', None, None)),
     ('C', (1, 0, None), ('d2', 'a1', None)),
     ('H', (2, 1, 0), ('d3', 'a2', 't1')),
     ('H', (2, 1, 0), ('d3', 'a2', 't2')),
     ('H', (2, 1, 0), ('d3', 'a2', 't3'))),
    {'d1': 1.70075351, 'a1': 1.74532925, 'a2': 1.83259571,
     't1': 1.04719755, 't2': -1.04719755})

CH4O_VAR_ZMA_STR = """
H
O 1    d1
C 2    d2 1    a1
H 3    d3 2    a2 1    t1
H 3    d3 2    a2 1    t2
H 3    d3 2    a2 1    t3

d1   =   0.900000
a1   = 100.000000
a2   = 105.000000
t1   =  60.000000
t2   = -60.000000
"""


def test__from_data():
    """ test zmatrix.from_data
    """
    zma = zmatrix.from_data(
        symbols=zmatrix.symbols(CH4O_ZMA),
        key_matrix=zmatrix.key_matrix(CH4O_ZMA),
        name_matrix=zmatrix.name_matrix(CH4O_ZMA),
        values=zmatrix.values(CH4O_ZMA),
    )
    assert zmatrix.almost_equal(zma, CH4O_ZMA)

    zma = zmatrix.from_data(
        symbols=zmatrix.symbols(CH4O_VAR_ZMA),
        key_matrix=zmatrix.key_matrix(CH4O_VAR_ZMA),
        name_matrix=zmatrix.name_matrix(CH4O_VAR_ZMA),
        values=zmatrix.values(CH4O_VAR_ZMA),
        complete=False
    )
    assert zmatrix.almost_equal(zma, CH4O_VAR_ZMA)


def test__coordinate_key_matrix():
    """ test zmatrix.coordinate_key_matrix
    """
    coo_key_mat = ((None, None, None),
                   ((1, 0), None, None),
                   ((2, 1), (2, 1, 0), None),
                   ((3, 2), (3, 2, 1), (3, 2, 1, 0)),
                   ((4, 2), (4, 2, 1), (4, 2, 1, 0)),
                   ((5, 2), (5, 2, 1), (5, 2, 1, 0)))

    assert zmatrix.coordinate_key_matrix(CH4O_ZMA) == coo_key_mat
    assert zmatrix.coordinate_key_matrix(CH4O_VAR_ZMA) == coo_key_mat

    one_indexed_coo_key_mat = ((None, None, None),
                               ((2, 1), None, None),
                               ((3, 2), (3, 2, 1), None),
                               ((4, 3), (4, 3, 2), (4, 3, 2, 1)),
                               ((5, 3), (5, 3, 2), (5, 3, 2, 1)),
                               ((6, 3), (6, 3, 2), (6, 3, 2, 1)))

    assert (zmatrix.coordinate_key_matrix(CH4O_ZMA, one_indexed=True)
            == one_indexed_coo_key_mat)
    assert (zmatrix.coordinate_key_matrix(CH4O_VAR_ZMA, one_indexed=True)
            == one_indexed_coo_key_mat)


def test__distance_names():
    """ test zmatrix.distance_names
    """
    assert zmatrix.distance_names(CH4O_ZMA) == ('d1', 'd2', 'd3')
    assert zmatrix.distance_names(CH4O_VAR_ZMA) == ('d1', 'd2', 'd3')
    assert zmatrix.distance_value_names(CH4O_ZMA) == ('d1', 'd2', 'd3')
    assert zmatrix.distance_value_names(CH4O_VAR_ZMA) == ('d1',)
    assert zmatrix.distance_variable_names(CH4O_ZMA) == ()
    assert zmatrix.distance_variable_names(CH4O_VAR_ZMA) == ('d2', 'd3')


def test__central_angle_names():
    """ test zmatrix.angle_names
    """
    assert zmatrix.central_angle_names(CH4O_ZMA) == ('a1', 'a2')
    assert zmatrix.central_angle_names(CH4O_VAR_ZMA) == ('a1', 'a2')
    assert zmatrix.central_angle_value_names(CH4O_ZMA) == ('a1', 'a2')
    assert zmatrix.central_angle_value_names(CH4O_VAR_ZMA) == ('a1', 'a2')
    assert zmatrix.central_angle_variable_names(CH4O_ZMA) == ()
    assert zmatrix.central_angle_variable_names(CH4O_VAR_ZMA) == ()


def test__dihedral_angle_names():
    """ test zmatrix.dihedral_names
    """
    assert zmatrix.dihedral_angle_names(CH4O_ZMA) == ('t1', 't2', 't3')
    assert zmatrix.dihedral_angle_names(CH4O_VAR_ZMA) == ('t1', 't2', 't3')
    assert zmatrix.dihedral_angle_value_names(CH4O_ZMA) == ('t1', 't2', 't3')
    assert zmatrix.dihedral_angle_value_names(CH4O_VAR_ZMA) == ('t1', 't2')
    assert zmatrix.dihedral_angle_variable_names(CH4O_ZMA) == ()
    assert zmatrix.dihedral_angle_variable_names(CH4O_VAR_ZMA) == ('t3',)


def test__angle_names():
    """ test zmatrix.angle_names
    """
    assert zmatrix.angle_names(CH4O_ZMA) == ('a1', 'a2', 't1', 't2', 't3')
    assert zmatrix.angle_names(CH4O_VAR_ZMA) == ('a1', 'a2', 't1', 't2', 't3')
    assert zmatrix.angle_value_names(CH4O_ZMA) == (
        'a1', 'a2', 't1', 't2', 't3')
    assert zmatrix.angle_value_names(CH4O_VAR_ZMA) == ('a1', 'a2', 't1', 't2')
    assert zmatrix.angle_variable_names(CH4O_ZMA) == ()
    assert zmatrix.angle_variable_names(CH4O_VAR_ZMA) == ('t3',)


def test__coordinates():
    """ test zmatrix.coordinates
    """
    coo_dct = {
        'd1': ((1, 0),), 'd2': ((2, 1),), 'd3': ((3, 2), (4, 2), (5, 2)),
        'a1': ((2, 1, 0),), 'a2': ((3, 2, 1), (4, 2, 1), (5, 2, 1)),
        't1': ((3, 2, 1, 0),), 't2': ((4, 2, 1, 0),), 't3': ((5, 2, 1, 0),)
    }
    assert zmatrix.coordinates(CH4O_ZMA) == coo_dct
    assert zmatrix.coordinates(CH4O_VAR_ZMA) == coo_dct

    one_indexed_coo_dct = {
        'd1': ((2, 1),), 'd2': ((3, 2),), 'd3': ((4, 3), (5, 3), (6, 3)),
        'a1': ((3, 2, 1),), 'a2': ((4, 3, 2), (5, 3, 2), (6, 3, 2)),
        't1': ((4, 3, 2, 1),), 't2': ((5, 3, 2, 1),), 't3': ((6, 3, 2, 1),)}
    assert (zmatrix.coordinates(CH4O_ZMA, one_indexed=True) ==
            one_indexed_coo_dct)
    assert (zmatrix.coordinates(CH4O_VAR_ZMA, one_indexed=True) ==
            one_indexed_coo_dct)


def test__set_names():
    """ test zmatrix.set_names
    """
    zma = zmatrix.set_names(CH4O_ZMA, {'d3': 'dCH'})
    assert zmatrix.distance_names(zma) == ('d1', 'd2', 'dCH')
    zma = zmatrix.set_names(CH4O_VAR_ZMA, {'d3': 'dCH'})
    assert zmatrix.distance_names(zma) == ('d1', 'd2', 'dCH')


def test__set_values():
    """ test zmatrix.set_values
    """
    val = zmatrix.values(CH4O_ZMA)['t1']

    zma = zmatrix.set_values(CH4O_ZMA, {'t1': val + 1e-6})
    assert zmatrix.almost_equal(zma, CH4O_ZMA)

    zma = zmatrix.set_values(CH4O_ZMA, {'t1': val + 1e-1})
    assert not zmatrix.almost_equal(zma, CH4O_ZMA)

    val = zmatrix.values(CH4O_VAR_ZMA)['t1']

    zma = zmatrix.set_values(CH4O_VAR_ZMA, {'t1': val + 1e-6})
    assert zmatrix.almost_equal(zma, CH4O_VAR_ZMA)

    zma = zmatrix.set_values(CH4O_VAR_ZMA, {'t1': val + 1e-1})
    assert not zmatrix.almost_equal(zma, CH4O_VAR_ZMA)


def test__is_valid():
    """ test zmatrix.is_valid
    """
    assert zmatrix.is_valid(CH4O_ZMA)
    assert not zmatrix.is_valid(zmatrix.geometry(CH4O_ZMA))
    assert not zmatrix.is_valid(CH4O_VAR_ZMA)
    assert zmatrix.is_valid(CH4O_VAR_ZMA, complete=False)


def test__is_complete():
    """ test zmatrix.is_complete
    """
    assert zmatrix.is_complete(CH4O_ZMA)
    assert not zmatrix.is_complete(CH4O_VAR_ZMA)


def test__is_standard_form():
    """ test zmatrix.is_standard_form
    """
    assert not zmatrix.is_standard_form(CH4O_ZMA)
    assert not zmatrix.is_standard_form(CH4O_VAR_ZMA)
    assert zmatrix.is_standard_form(zmatrix.standard_form(CH4O_ZMA))
    assert zmatrix.is_standard_form(zmatrix.standard_form(CH4O_VAR_ZMA))
    print(zmatrix.standard_form(CH4O_ZMA))
    print(zmatrix.standard_form(CH4O_VAR_ZMA))


def test__from_zmat_string():
    """ test zmatrix.from_zmat_string
    """
    zma = zmatrix.from_zmat_string(CH4O_ZMA_STR)
    assert zmatrix.almost_equal(zma, CH4O_ZMA)

    zma = zmatrix.from_zmat_string(H2O_ZMA_STR)
    assert zmatrix.almost_equal(zma, H2O_ZMA)

    zma = zmatrix.from_zmat_string(CH4O_VAR_ZMA_STR, complete=False)
    assert zmatrix.almost_equal(zma, CH4O_VAR_ZMA)


def test__zmat_string():
    """ test zmatrix.zmat_string
    """
    zma = zmatrix.from_zmat_string(zmatrix.zmat_string(CH4O_ZMA))
    assert zmatrix.almost_equal(zma, CH4O_ZMA)
    zma = zmatrix.from_zmat_string(zmatrix.zmat_string(CH4O_VAR_ZMA),
                                   complete=False)
    assert zmatrix.almost_equal(zma, CH4O_VAR_ZMA)


def test__geometry():
    """ test zmatrix.geometry
    """
    ref_zma = automol.geom.zmatrix(zmatrix.geometry(CH4O_ZMA))
    zma = automol.geom.zmatrix(zmatrix.geometry(ref_zma))
    assert zmatrix.almost_equal(zma, ref_zma)


def test__connectivity_graph():
    """ test zmatrix.connectivity_graph
    """
    gra = zmatrix.connectivity_graph(CH4O_ZMA)
    assert gra == CH4O_GRA


def test__tors__symmetry_numbers():
    """ test zmatrix.tors.symmetry_numbers
    """
    assert zmatrix.tors.symmetry_numbers(CH4O_ZMA, ('t1',)) == (3,)
    assert zmatrix.tors.symmetry_numbers(CH4O2_ZMA, ('D5', 'D6')) == (1, 1)


def test__tors__samples():
    """ test zmatrix.tors.samples
    """
    zmas = zmatrix.tors.samples(CH4O_ZMA, 9, ('t1',))
    assert len(zmas) == 9

    zmas = zmatrix.tors.samples(CH4O2_ZMA, 7, ('D5', 'D6'))
    assert len(zmas) == 7


if __name__ == '__main__':
    test__distance_names()
    test__central_angle_names()
    test__dihedral_angle_names()
    test__angle_names()
    # test__set_names()
    # test__set_values()
    test__from_zmat_string()
    # test__zmat_string()
    # test__geometry()
    # test__connectivity_graph()
    # test__coordinate_key_matrix()
    # test__tors__symmetry_numbers()
    # test__tors__samples()
    test__is_valid()
    # test__coordinates()
    # test__from_data()
    test__is_complete()
    test__is_standard_form()
