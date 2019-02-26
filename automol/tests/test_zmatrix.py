""" test the automol.zmatrix module
"""
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


def test__distance_names():
    """ test zmatrix.distance_names
    """
    assert zmatrix.distance_names(CH4O_ZMA) == ('d1', 'd2', 'd3')


def test__angle_names():
    """ test zmatrix.angle_names
    """
    assert zmatrix.angle_names(CH4O_ZMA) == ('a1', 'a2')


def test__torsion_names():
    """ test zmatrix.torsion_names
    """
    assert zmatrix.torsion_names(CH4O_ZMA) == ('t1', 't2', 't3')


def test__set_names():
    """ test zmatrix.set_names
    """
    zma = zmatrix.set_names(CH4O_ZMA, {'d3': 'dCH'})
    assert zmatrix.distance_names(zma) == ('d1', 'd2', 'dCH')


def test__set_values():
    """ test zmatrix.set_values
    """
    val = zmatrix.values(CH4O_ZMA)['t1']

    zma = zmatrix.set_values(CH4O_ZMA, {'t1': val + 1e-6})
    assert zmatrix.almost_equal(zma, CH4O_ZMA)

    zma = zmatrix.set_values(CH4O_ZMA, {'t1': val + 1e-1})
    assert not zmatrix.almost_equal(zma, CH4O_ZMA)


def test__from_zmat_string():
    """ test zmatrix.from_zmat_string
    """
    zma = zmatrix.from_zmat_string(CH4O_ZMA_STR)
    assert zmatrix.almost_equal(zma, CH4O_ZMA)


def test__zmat_string():
    """ test zmatrix.zmat_string
    """
    zma = zmatrix.from_zmat_string(zmatrix.zmat_string(CH4O_ZMA))
    assert zmatrix.almost_equal(zma, CH4O_ZMA)


if __name__ == '__main__':
    test__from_data()
    test__distance_names()
    test__angle_names()
    test__torsion_names()
    test__set_names()
    test__set_values()
    test__from_zmat_string()
    test__zmat_string()
