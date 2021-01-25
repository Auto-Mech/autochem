""" test automol.vmat
"""
from automol import vmat

CH4O_VMA = (('H', (None, None, None), (None, None, None)),
            ('O', (0, None, None), ('r1', None, None)),
            ('C', (1, 0, None), ('r2', 'a1', None)),
            ('H', (2, 1, 0), ('r3', 'a2', 'd1')),
            ('H', (2, 1, 0), ('r3', 'a2', 'd2')),
            ('H', (2, 1, 0), ('r3', 'a2', 'd3')))

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

CH4O_VMA_STR = """
H
O 1    r1
C 2    r2 1    a1
H 3    r3 2    a2 1    d1
H 3    r3 2    a2 1    d2
H 3    r3 2    a2 1    d3
"""


def test__from_data():
    """ test vmat.from_data
    """
    vma = vmat.from_data(
        symbs=vmat.symbols(CH4O_VMA),
        key_mat=vmat.key_matrix(CH4O_VMA),
        name_mat=vmat.name_matrix(CH4O_VMA),
    )
    assert vma == CH4O_VMA


def test__from_string():
    """ test vmat.from_string
    """
    vma = vmat.from_string(CH4O_VMA_STR)
    assert vma == CH4O_VMA


def test__string():
    """ test vmat.string
    """
    vma = vmat.from_string(vmat.string(CH4O_VMA))
    assert vma == CH4O_VMA
