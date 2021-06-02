
""" test automol.zmat
"""
import numpy
from automol import zmat

CH4O2_ZMA = (
    ('C', (None, None, None), (None, None, None), (None, None, None)),
    ('O', (0, None, None), ('R1', None, None), (2.659, None, None)),
    ('O', (0, 1, None), ('R2', 'A2', None), (2.659, 1.907, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'), (2.068, 1.933, 4.14)),
    ('H', (0, 1, 2), ('R4', 'A4', 'D4'), (2.065, 1.894, 2.06)),
    ('H', (1, 0, 2), ('R5', 'A5', 'D5'), (1.831, 1.867, 1.44)),
    ('H', (2, 0, 1), ('R6', 'A6', 'D6'), (1.831, 1.867, 4.84)))

CH4O2_ZMA_NO_NONES = (
    ('C', (), (), ()),
    ('O', (0,), ('R1',), (2.659,)),
    ('O', (0, 1,), ('R2', 'A2',), (2.659, 1.907,)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'), (2.068, 1.933, 4.14)),
    ('H', (0, 1, 2), ('R4', 'A4', 'D4'), (2.065, 1.894, 2.06)),
    ('H', (1, 0, 2), ('R5', 'A5', 'D5'), (1.831, 1.867, 1.44)),
    ('H', (2, 0, 1), ('R6', 'A6', 'D6'), (1.831, 1.867, 4.84)))

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
    """ test zmat.from_data
    """
    zma1 = zmat.from_data(
        symbs=zmat.symbols(CH4O2_ZMA),
        key_mat=zmat.key_matrix(CH4O2_ZMA),
        val_mat=zmat.value_matrix(CH4O2_ZMA),
        name_mat=zmat.name_matrix(CH4O2_ZMA),
    )
    assert zma1 == CH4O2_ZMA

    zma2 = zmat.from_data(
        symbs=zmat.symbols(CH4O2_ZMA_NO_NONES),
        key_mat=zmat.key_matrix(CH4O2_ZMA_NO_NONES),
        val_mat=zmat.value_matrix(CH4O2_ZMA_NO_NONES),
        name_mat=zmat.name_matrix(CH4O2_ZMA_NO_NONES),
    )
    assert zma2 == CH4O2_ZMA

    zma1 = list(map(list, zma1))
    zma2 = list(map(list, zma2))

    assert zmat.is_valid(zma1)
    assert zmat.is_valid(zma2)

    zma1[0] += [None]
    zma2[0][1] = zma2[0][1] + (None,)

    assert not zmat.is_valid(zma1)
    assert not zmat.is_valid(zma2)


def test__add_atom():
    """ test zmat.add_atom

    output below should look like this:
        C
        O  1    R2
        C  2    R3  1    A3
        H  3    R4  2    A4  1    D4
        H  3    R5  2    A5  1    D5
        H  3    R6  2    A6  1    D6

        R2   =   1.500000
        R3   =   1.500000
        R4   =   1.500000
        R5   =   1.500000
        R6   =   1.500000
        A3   = 120.000000
        A4   = 120.000000
        A5   = 120.000000
        A6   = 120.000000
        D4   =   0.000000
        D5   = 109.470000
        D6   =-109.470000
    """
    zma = ()
    zma = zmat.add_atom(zma, 'C', (), ())
    zma = zmat.add_atom(zma, 'O', (0,), (1.5,))
    zma = zmat.add_atom(zma, 'C', (1, 0), (1.5, 120.))
    zma = zmat.add_atom(zma, 'H', (2, 1, 0), (1.5, 120., 0.))
    zma = zmat.add_atom(zma, 'H', (2, 1, 0), (1.5, 120., +109.47))
    zma = zmat.add_atom(zma, 'H', (2, 1, 0), (1.5, 120., -109.47))
    print(zmat.string(zma))


def test__string():
    """ test zmat.string
    """
    zma = zmat.from_string(zmat.string(CH4O2_ZMA))
    assert zmat.almost_equal(zma, CH4O2_ZMA)


def test__distance():
    """ test zmat.distance
    """
    dist16 = zmat.distance(CH4O2_ZMA, 1, 6, angstrom=True)
    print(dist16)
    assert numpy.isclose(dist16, 2.6363)


if __name__ == '__main__':
    test__from_data()
    test__string()
    test__add_atom()
    test__distance()
