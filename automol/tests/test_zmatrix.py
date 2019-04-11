""" test the automol.zmatrix module
"""
import automol
from automol import zmatrix


CH4O_ZMA = (
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

CH4O_ZMA_STR = """
H
O 1    r1
C 2    r2 1    a1
H 3    r3 2    a2 1    d1
H 3    r3 2    a2 1    d2
H 3    r3 2    a2 1    d3

r1   =   0.900000
r2   =   1.400000
a1   = 100.000000
r3   =   1.100000
a2   = 105.000000
d1   =  60.000000
d2   = -60.000000
d3   = 179.999997
"""

CH4O2_ZMA = (
    (('C', (None, None, None), (None, None, None)),
     ('O', (0, None, None), ('r1', None, None)),
     ('O', (0, 1, None), ('r2', 'a1', None)),
     ('H', (0, 1, 2), ('r3', 'a2', 'd1')),
     ('H', (0, 1, 2), ('r4', 'a3', 'd2')),
     ('H', (1, 0, 2), ('r5', 'a4', 'd3')),
     ('H', (2, 0, 1), ('r6', 'a5', 'd4'))),
    {'r1': 2.65933,
     'r2': 2.65933, 'a1': 1.90743,
     'r3': 2.06844, 'a2': 1.93366, 'd1': 4.1477,
     'r4': 2.06548, 'a3': 1.89469, 'd2': 2.06369,
     'r5': 1.83126, 'a4': 1.86751, 'd3': 1.44253,
     'r6': 1.83126, 'a5': 1.86751, 'd4': 4.84065})


def test__from_data():
    """ test zmatrix.from_data
    """
    vma = zmatrix.v.from_data(
        symbols=zmatrix.v.symbols(zmatrix.var_(CH4O_ZMA)),
        key_matrix=zmatrix.v.key_matrix(zmatrix.var_(CH4O_ZMA)),
        name_matrix=zmatrix.v.name_matrix(zmatrix.var_(CH4O_ZMA)),
    )
    assert vma == zmatrix.var_(CH4O_ZMA)

    zma = zmatrix.from_data(
        symbols=zmatrix.symbols(CH4O_ZMA),
        key_matrix=zmatrix.key_matrix(CH4O_ZMA),
        name_matrix=zmatrix.name_matrix(CH4O_ZMA),
        values=zmatrix.values(CH4O_ZMA),
    )
    assert zma == CH4O_ZMA


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

    assert (zmatrix.coordinate_key_matrix(CH4O_ZMA, one_indexed=True)
            == one_indexed_coo_key_mat)


def test__coordinates():
    """ test zmatrix.v.coordinates
    """
    coo_dct = {
        'r1': ((1, 0),), 'r2': ((2, 1),), 'r3': ((3, 2), (4, 2), (5, 2)),
        'a1': ((2, 1, 0),), 'a2': ((3, 2, 1), (4, 2, 1), (5, 2, 1)),
        'd1': ((3, 2, 1, 0),), 'd2': ((4, 2, 1, 0),), 'd3': ((5, 2, 1, 0),)
    }
    assert zmatrix.coordinates(CH4O_ZMA) == coo_dct

    one_indexed_coo_dct = {
        'r1': ((2, 1),), 'r2': ((3, 2),), 'r3': ((4, 3), (5, 3), (6, 3)),
        'a1': ((3, 2, 1),), 'a2': ((4, 3, 2), (5, 3, 2), (6, 3, 2)),
        'd1': ((4, 3, 2, 1),), 'd2': ((5, 3, 2, 1),), 'd3': ((6, 3, 2, 1),)}
    assert (zmatrix.coordinates(CH4O_ZMA, one_indexed=True) ==
            one_indexed_coo_dct)


def test__distance_names():
    """ test zmatrix.distance_names
    """
    assert zmatrix.distance_names(CH4O_ZMA) == ('r1', 'r2', 'r3')


def test__central_angle_names():
    """ test zmatrix.angle_names
    """
    assert zmatrix.central_angle_names(CH4O_ZMA) == ('a1', 'a2')


def test__dihedral_angle_names():
    """ test zmatrix.dihedral_names
    """
    assert zmatrix.dihedral_angle_names(CH4O_ZMA) == ('d1', 'd2', 'd3')


def test__angle_names():
    """ test zmatrix.v.angle_names
    """
    assert zmatrix.angle_names(CH4O_ZMA) == ('a1', 'a2', 'd1', 'd2', 'd3')


def test__from_string():
    """ test zmatrix.v.from_vmat_string
    """
    vma = zmatrix.v.from_string(CH4O_ZMA_STR)
    assert vma == zmatrix.var_(CH4O_ZMA)

    zma = zmatrix.from_string(CH4O_ZMA_STR)
    assert zmatrix.almost_equal(zma, CH4O_ZMA)


def test__string():
    """ test zmatrix.v.vmat_string
    """
    vma = zmatrix.v.from_string(zmatrix.v.string(zmatrix.var_(CH4O_ZMA)))
    assert vma == zmatrix.var_(CH4O_ZMA)

    zma = zmatrix.from_string(zmatrix.string(CH4O_ZMA))
    assert zmatrix.almost_equal(zma, CH4O_ZMA)


def test__set_values():
    """ test zmatrix.set_values
    """
    val = zmatrix.values(CH4O_ZMA)['d1']

    zma = zmatrix.set_values(CH4O_ZMA, {'d1': val + 1e-6})
    assert zmatrix.almost_equal(zma, CH4O_ZMA)

    zma = zmatrix.set_values(CH4O_ZMA, {'d1': val + 1e-1})
    assert not zmatrix.almost_equal(zma, CH4O_ZMA)


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
         ('O', (0, None, None), ('R1', None, None)),
         ('C', (1, 0, None), ('R2', 'A2', None)),
         ('H', (2, 1, 0), ('R3', 'A3', 'D3')),
         ('H', (2, 1, 0), ('R3', 'A3', 'D4')),
         ('H', (2, 1, 0), ('R3', 'A3', 'D5'))),
        {'R1': 1.70075351,
         'R2': 2.64561657, 'A2': 1.74532925,
         'R3': 2.07869873, 'A3': 1.83259571,
         'D3': 1.04719755, 'D4': -1.04719755, 'D5': 3.1415926})
    zma = zmatrix.standard_form(nonst_zma)

    assert zmatrix.is_valid(nonst_zma)
    assert not zmatrix.is_standard_form(nonst_zma)
    assert zmatrix.is_standard_form(zma)
    assert zmatrix.almost_equal(zma, CH4O_ZMA)


def test__geometry():
    """ test zmatrix.geometry
    """
    ref_zma = (
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
    ref_zma = zmatrix.standard_form(ref_zma)

    zma = automol.geom.zmatrix(zmatrix.geometry(ref_zma))
    zma = zmatrix.standard_form(zma)
    assert zmatrix.almost_equal(zma, ref_zma)


def test__connectivity_graph():
    """ test zmatrix.connectivity_graph
    """
    ref_gra = (
        {0: ('H', 0, None), 1: ('O', 0, None), 2: ('C', 0, None),
         3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None)},
        {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
         frozenset({2, 3}): (1, None), frozenset({2, 4}): (1, None),
         frozenset({2, 5}): (1, None)})

    gra = zmatrix.connectivity_graph(CH4O_ZMA)
    assert gra == ref_gra


def test__tors__symmetry_numbers():
    """ test zmatrix.tors.symmetry_numbers
    """
    assert zmatrix.tors.symmetry_numbers(CH4O_ZMA, ('d1',)) == (3,)
    assert zmatrix.tors.symmetry_numbers(CH4O2_ZMA, ('d3', 'd4')) == (1, 1)


def test__tors__samples():
    """ test zmatrix.tors.samples
    """
    zmas = zmatrix.tors.samples(CH4O_ZMA, 9, ('d1',))
    assert len(zmas) == 9

    zmas = zmatrix.tors.samples(CH4O2_ZMA, 7, ('d3', 'd4'))
    assert len(zmas) == 7


if __name__ == '__main__':
    test__from_data()
    test__coordinate_key_matrix()
    test__coordinates()
    test__from_string()
    test__string()
    test__set_values()
    test__is_valid()
    test__is_standard_form()
    test__geometry()
    test__tors__symmetry_numbers()
    test__tors__samples()
