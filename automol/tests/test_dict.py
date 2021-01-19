""" test automol.dict_ not utilized by other tests
"""

import numpy
from automol import dict_


DCT1 = {
    'ich': 'InChI=1S/N2/c1-2',
    'mult': 2
}
DCT2 = {
    'mult': 1,
    'charge': 0,
}
DCT3 = {
    'mult': 3,
    'charge': 0,
    'linear': True
}
DCT4 = {
    ('C', 'H'): 1.54
}


def test__build():
    """ test automol.dict_
    """

    # Combine dct functions
    ref_dct1 = {
        'ich': 'InChI=1S/N2/c1-2',
        'charge': 0,
        'mult': 1
    }

    dct1 = dict_.right_update(DCT1, DCT2)

    assert dct1 == ref_dct1

    ref_dct2 = {
        'ich': 'InChI=1S/N2/c1-2',
        'charge': 0,
        'mult': 3,
        'linear': True
    }
    dct2 = dict_.merge_sequence((DCT1, DCT2, DCT3))

    assert dct2 == ref_dct2


def test__read():
    """ test read
    """

    val1 = dict_.values_by_unordered_tuple(DCT4, ('C', 'H'))
    val2 = dict_.values_by_unordered_tuple(DCT4, ('H', 'C'))

    assert numpy.isclose(val1, val2, 1.54)
