""" test the automol.zmatrix module
"""
from automol import zmatrix
from automol import geom

H_ZMA = (('H', (None, None, None), (None, None, None)),)

H_ZMA_STR = """charge: 0, mult: 2
H
"""

HO_ZMA = (('H', (None, None, None), (None, None, None)),
          ('O', (0, None, None), (0.47625948, None, None)))

HO_ZMA_STR = """charge: 0, mult: 2
H
O 1 0.9
"""

HO2_ZMA = (('H', (None, None, None), (None, None, None)),
           ('O', (0, None, None), (0.47625948, None, None)),
           ('O', (1, 0, None), (0.74084809, 1.83259571, None)))

HO2_ZMA_STR = """charge: 0, mult: 2
H
O 1 0.9
O 2 1.4 1 105.0
"""

H2O2_ZMA = (('H', (None, None, None), (None, None, None)),
            ('O', (0, None, None), (0.47625948, None, None)),
            ('O', (1, 0, None), (0.74084809, 1.83259571, None)),
            ('H', (2, 1, 0), (0.47625948, 1.83259571, 2.09439510)))

H2O2_ZMA_STR = """charge: 0, mult: 1
H
O 1 0.9
O 2 1.4 1 105.0
H 3 0.9 2 105.0 1 120.0
"""

CH4O_ZMA = (('H', (None, None, None), (None, None, None)),
            ('O', (0, None, None), (1.70075351, None, None)),
            ('C', (1, 0, None), (2.64561657, 1.74532925, None)),
            ('H', (2, 1, 0), (2.07869873, 1.83259571, 1.04719755)),
            ('H', (2, 1, 0), (2.07869873, 1.83259571, -1.04719755)),
            ('H', (2, 1, 0), (2.07869873, 1.83259571, 3.1415926)))

CH4O_VAR_DCT = {(1, 0): 'rOH', (3, 2): 'rCH', (4, 2): 'rCH', (5, 2): 'rCH',
                (3, 2, 1, 0): 'd1', (4, 2, 1, 0): 'd2', (5, 2, 1, 0): 'd3'}

CH4O_ZMA_STR = """charge: 0, mult: 1
H
O 1 rOH
C 2 1.4 1 100.0
H 3 rCH 2 105.0 1 d1
H 3 rCH 2 105.0 1 d2
H 3 rCH 2 105.0 1 d3

rOH = 0.9
rCH = 1.1
d1 = 60.0
d2 = -60.0
d3 = 180.0
"""


def test__from_data():
    """ test zmatrix.from_data
    """
    assert geom.is_valid(H_ZMA) is False
    assert geom.is_valid(HO_ZMA) is False
    assert geom.is_valid(HO2_ZMA) is False
    assert geom.is_valid(H2O2_ZMA) is False
    assert geom.is_valid(CH4O_ZMA) is False

    assert H_ZMA == zmatrix.from_data(
        symbols=zmatrix.symbols(H_ZMA),
        distance_column=zmatrix.distance_column(H_ZMA),
        angle_column=zmatrix.angle_column(H_ZMA),
        torsion_column=zmatrix.torsion_column(H_ZMA),
        angstroms=False,
        one_indexed=False,
    )
    assert HO_ZMA == zmatrix.from_data(
        symbols=zmatrix.symbols(HO_ZMA),
        distance_column=zmatrix.distance_column(HO_ZMA),
        angle_column=zmatrix.angle_column(HO_ZMA),
        torsion_column=zmatrix.torsion_column(HO_ZMA),
        angstroms=False,
        one_indexed=False,
    )
    assert HO2_ZMA == zmatrix.from_data(
        symbols=zmatrix.symbols(HO2_ZMA),
        distance_column=zmatrix.distance_column(HO2_ZMA),
        angle_column=zmatrix.angle_column(HO2_ZMA),
        torsion_column=zmatrix.torsion_column(HO2_ZMA),
        angstroms=False,
        one_indexed=False,
    )
    assert H2O2_ZMA == zmatrix.from_data(
        symbols=zmatrix.symbols(H2O2_ZMA),
        distance_column=zmatrix.distance_column(H2O2_ZMA),
        angle_column=zmatrix.angle_column(H2O2_ZMA),
        torsion_column=zmatrix.torsion_column(H2O2_ZMA),
        angstroms=False,
        one_indexed=False,
    )
    assert CH4O_ZMA == zmatrix.from_data(
        symbols=zmatrix.symbols(CH4O_ZMA),
        distance_column=zmatrix.distance_column(CH4O_ZMA),
        angle_column=zmatrix.angle_column(CH4O_ZMA),
        torsion_column=zmatrix.torsion_column(CH4O_ZMA),
        angstroms=False,
        one_indexed=False,
    )


def test__coordinate_matrix():
    """ test zmatrix.coordinate_matrix
    """
    assert zmatrix.coordinate_matrix(H_ZMA) == (
        ((None, None, None),)
    )
    assert zmatrix.coordinate_matrix(HO_ZMA) == (
        ((None, None, None),
         ((1, 0), None, None))
    )
    assert zmatrix.coordinate_matrix(HO2_ZMA) == (
        ((None, None, None),
         ((1, 0), None, None),
         ((2, 1), (2, 1, 0), None))
    )
    assert zmatrix.coordinate_matrix(H2O2_ZMA) == (
        ((None, None, None),
         ((1, 0), None, None),
         ((2, 1), (2, 1, 0), None),
         ((3, 2), (3, 2, 1), (3, 2, 1, 0)))
    )
    assert zmatrix.coordinate_matrix(CH4O_ZMA) == (
        ((None, None, None),
         ((1, 0), None, None),
         ((2, 1), (2, 1, 0), None),
         ((3, 2), (3, 2, 1), (3, 2, 1, 0)),
         ((4, 2), (4, 2, 1), (4, 2, 1, 0)),
         ((5, 2), (5, 2, 1), (5, 2, 1, 0)))
    )


def test__from_zmat_string():
    """ test zmatrix.from_zmat_string
    """
    zma, var_dct = zmatrix.from_zmat_string(CH4O_ZMA_STR)
    assert zmatrix.almost_equal(zma, CH4O_ZMA)
    assert var_dct == CH4O_VAR_DCT


def test__zmat_string():
    """ test zmatrix.zmat_string
    """
    zma_str = zmatrix.zmat_string(CH4O_ZMA, CH4O_VAR_DCT)
    zma, var_dct = zmatrix.from_zmat_string(zma_str)
    assert zmatrix.almost_equal(zma, CH4O_ZMA)
    assert var_dct == CH4O_VAR_DCT


if __name__ == '__main__':
    test__from_data()
    test__coordinate_matrix()
    test__from_zmat_string()
    test__zmat_string()
