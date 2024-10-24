""" test automol.vmat
"""

import pytest
from automol import vmat


# VMAT Info for C#C[C@H](C)CO
C5H8O_VMA = (
    ('C', (None, None, None), (None, None, None)),
    ('C', (0, None, None), ('R1', None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
    ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
    ('C', (1, 0, 2), ('R5', 'A5', 'D5')),
    ('C', (1, 0, 5), ('R6', 'A6', 'D6')),
    ('H', (1, 0, 5), ('R7', 'A7', 'D7')),
    ('X', (5, 1, 0), ('R8', 'A8', 'D8')),
    ('C', (5, 8, 1), ('R9', 'A9', 'D9')),
    ('X', (9, 5, 8), ('R10', 'A10', 'D10')),
    ('H', (9, 10, 5), ('R11', 'A11', 'D11')),
    ('O', (6, 1, 0), ('R12', 'A12', 'D12')),
    ('H', (6, 1, 12), ('R13', 'A13', 'D13')),
    ('H', (6, 1, 12), ('R14', 'A14', 'D14')),
    ('H', (12, 6, 1), ('R15', 'A15', 'D15')))

C5H8O_VMA_STR = """
C
C  1    R1
H  1    R2  2    A2
H  1    R3  2    A3  3    D3
H  1    R4  2    A4  3    D4
C  2    R5  1    A5  3    D5
C  2    R6  1    A6  6    D6
H  2    R7  1    A7  6    D7
X  6    R8  2    A8  1    D8
C  6    R9  9    A9  2    D9
X  10   R10  6   A10  9   D10
H  10   R11  11   A11  6   D11
O  7   R12  2   A12  1   D12
H  7   R13  2   A13  13   D13
H  7   R14  2   A14  13   D14
H  13   R15  7   A15  2   D15"""


BAD_C5H8O_VMA = (
    ('C', (None, None, None)),
    ('C', (0, None, None)),
    ('H', (0, 1, None)),
    ('H', (0, 1, 2)),
    ('H', (0, 1, 2)),
    ('C', (1, 0, 2)),
    ('C', (1, 0, 5)),
    ('H', (1, 0, 5)),
    ('X', (5, 1, 0)),
    ('C', (5, 8, 1)),
    ('X', (9, 5, 8)),
    ('H', (9, 10, 5)),
    ('O', (6, 1, 0)),
    ('H', (6, 1, 12)),
    ('H', (6, 1, 12)),
    ('H', (12, 6, 1)))


def test__from_data():
    """ test vmat.from_data
    """
    vma = vmat.from_data(
        symbs=vmat.symbols(C5H8O_VMA),
        key_mat=vmat.key_matrix(C5H8O_VMA),
        name_mat=vmat.name_matrix(C5H8O_VMA),
    )
    assert vma == C5H8O_VMA


def test__atom_indices():
    """ test vmat.atom_indices
    """

    angle_idxs = vmat.atom_indices(C5H8O_VMA, 'C', match=False)
    assert angle_idxs == (2, 3, 4, 7, 8, 10, 11, 12, 13, 14, 15)


def test__coordinates():
    """ test vmat.angle_names
    """

    names = vmat.names(C5H8O_VMA)
    angle_names = vmat.angle_names(C5H8O_VMA)

    assert names == ('R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9',
                     'R10', 'R11', 'R12', 'R13', 'R14', 'R15', 'A2', 'A3',
                     'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12',
                     'A13', 'A14', 'A15', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8',
                     'D9', 'D10', 'D11', 'D12', 'D13', 'D14', 'D15')
    assert angle_names == ('A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8',
                           'A9', 'A10', 'A11', 'A12', 'A13', 'A14',
                           'A15', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8',
                           'D9', 'D10', 'D11', 'D12', 'D13', 'D14', 'D15')


def test__standardize():
    """ test vmat.standard_form
        test vmat.standard_name_matrix
    """

    nonstandard_vma = (
        ('C', (None, None, None), (None, None, None)),
        ('C', (0, None, None), ('R20', None, None)),
        ('H', (0, 1, None), ('R30', 'A30', None)),
        ('H', (0, 1, 2), ('R40', 'A40', 'D40')),
        ('H', (0, 1, 2), ('R50', 'A50', 'D50')),
        ('C', (1, 0, 2), ('R60', 'A60', 'D60')),
        ('C', (1, 0, 5), ('R70', 'A40', 'D70')),
        ('H', (1, 0, 5), ('R80', 'A80', 'D80')),
        ('X', (5, 1, 0), ('R90', 'A90', 'D90')),
        ('C', (5, 8, 1), ('R100', 'A100', 'D100')),
        ('X', (9, 5, 8), ('R110', 'A110', 'D110')),
        ('H', (9, 10, 5), ('R120', 'A120', 'D120')),
        ('O', (6, 1, 0), ('R130', 'A130', 'D130')),
        ('H', (6, 1, 12), ('R140', 'A140', 'D140')),
        ('H', (6, 1, 12), ('R150', 'A150', 'D150')),
        ('H', (12, 6, 1), ('R160', 'A160', 'D160')))

    assert C5H8O_VMA == vmat.standard_form(nonstandard_vma)
    assert vmat.is_standard_form(C5H8O_VMA)


def test__remove_atom():
    """ test vmat.remove_atom
    """

    # Remove a terminal hydrogen
    ref_vma = (
        ('C', (None, None, None), (None, None, None)),
        ('C', (0, None, None), ('R1', None, None)),
        ('H', (0, 1, None), ('R2', 'A2', None)),
        ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
        ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
        ('C', (1, 0, 2), ('R5', 'A5', 'D5')),
        ('C', (1, 0, 5), ('R6', 'A6', 'D6')),
        ('H', (1, 0, 5), ('R7', 'A7', 'D7')),
        ('X', (5, 1, 0), ('R8', 'A8', 'D8')),
        ('C', (5, 8, 1), ('R9', 'A9', 'D9')),
        ('X', (9, 5, 8), ('R10', 'A10', 'D10')),
        ('H', (9, 10, 5), ('R11', 'A11', 'D11')),
        ('O', (6, 1, 0), ('R12', 'A12', 'D12')),
        ('H', (6, 1, 12), ('R14', 'A14', 'D14')),
        ('H', (12, 6, 1), ('R15', 'A15', 'D15')))
    assert ref_vma == vmat.remove_atom(C5H8O_VMA, 13)

    # remove O atom from OH group; should break
    with pytest.raises(ValueError):
        _ = vmat.remove_atom(C5H8O_VMA, 12)


def test__from_string():
    """ test vmat.from_string
    """
    assert vmat.from_string(C5H8O_VMA_STR) == C5H8O_VMA


def test__string():
    """ test vmat.string
    """
    assert vmat.from_string(vmat.string(C5H8O_VMA)) == C5H8O_VMA


def test_valid():
    """ test vmat.is_valid
    """

    assert vmat.is_valid(C5H8O_VMA)
    assert not vmat.is_valid(BAD_C5H8O_VMA)


test__remove_atom()
test__atom_indices()
test__coordinates()
test__standardize()
