""" test automol.vmat
"""
from automol import vmat

CH4O_VMA = (
    ('H', (None, None, None), (None, None, None),),
    ('O', (0, None, None), ('R1', None, None),),
    ('C', (1, 0, None), ('R2', 'A2', None),),
    ('H', (2, 1, 0), ('R3', 'A3', 'D3'),),
    ('H', (2, 1, 0), ('R4', 'A4', 'D4'),),
    ('H', (2, 1, 0), ('R5', 'A5', 'D5'),))

CH4O_VMA_NO_NONES = (
    ('H', (), (),),
    ('O', (0,), ('R1',),),
    ('C', (1, 0,), ('R2', 'A2',),),
    ('H', (2, 1, 0), ('R3', 'A3', 'D3'),),
    ('H', (2, 1, 0), ('R4', 'A4', 'D4'),),
    ('H', (2, 1, 0), ('R5', 'A5', 'D5'),))

CH4O_VMA_STR = """
H
O 1    R1
C 2    R2 1    A2
H 3    R3 2    A3 1    D3
H 3    R4 2    A4 1    D4
H 3    R5 2    A5 1    D5
"""


def test__from_data():
    """ test vmat.from_data
    """
    vma1 = vmat.from_data(
        syms=vmat.symbols(CH4O_VMA),
        key_mat=vmat.key_matrix(CH4O_VMA),
        name_mat=vmat.name_matrix(CH4O_VMA),
    )
    assert vma1 == CH4O_VMA

    vma2 = vmat.from_data(
        syms=vmat.symbols(CH4O_VMA_NO_NONES),
        key_mat=vmat.key_matrix(CH4O_VMA_NO_NONES),
        name_mat=vmat.name_matrix(CH4O_VMA_NO_NONES),
    )
    assert vma2 == CH4O_VMA

    vma1 = list(map(list, vma1))
    vma2 = list(map(list, vma2))

    assert vmat.is_valid(vma1)
    assert vmat.is_valid(vma2)

    vma1[0] += [None]
    vma2[0][1] = vma2[0][1] + (None,)

    assert not vmat.is_valid(vma1)
    assert not vmat.is_valid(vma2)


def test__from_string():
    """ test vmat.from_string
    """
    vma = vmat.from_string(CH4O_VMA_STR)
    print(vma)
    assert vma == CH4O_VMA


def test__string():
    """ test vmat.string
    """
    vma = vmat.from_string(vmat.string(CH4O_VMA))
    assert vma == CH4O_VMA


if __name__ == '__main__':
    test__from_data()
