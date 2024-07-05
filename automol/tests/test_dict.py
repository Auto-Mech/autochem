""" test automol.dict_ not utilized by other tests
"""

import numpy
from automol.util import dict_


# General
GDCT1 = {
    'key1': 'val1',
    'key2': 'val2',
    'key3': 'val3'
}
GDCT2 = {
    'key1': 'val1',
    'key4': 'val4'
}
GDCT3 = {
    'global': {
        'subkey1': 'subvalg1',
        'subkey3': 'subvalg3'
    },
    'key1': {
        'subkey1': 'subval1-1',
        'subkey2': 'subval1-2'
    },
    'key2': {
        'subkey1': 'subval2-1',
        'subkey2': 'subval2-2'
    }
}
GDCT4 = {
    'global': {
        'subkey1': 'subvalg1',
        'subkey3': 'subvalg3'
    }
}
GDCT5 = {'b': 2, 'd': 4, 'a': 1, 'e': 5, 'c': 3}


# Specific dictionaries
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
    ('C', 'H'): 1.54,
    ('C', 'O'): 1.42
}


def test__merge_sequence():
    """ test automol.dict_.merge_sequence
    """

    ref_dct2 = {
        'ich': 'InChI=1S/N2/c1-2',
        'charge': 0,
        'mult': 3,
        'linear': True
    }
    dct2 = dict_.merge_sequence((DCT1, DCT2, DCT3))

    assert dct2 == ref_dct2


def test__read():
    """ test dict_.values_by_unordered_tuple
        test dict_.values_in_multilevel
    """

    # Unordered unordered tuple reads
    val1 = dict_.value_by_unordered_key(DCT4, ('C', 'H'))
    val2 = dict_.value_by_unordered_key(DCT4, ('H', 'C'))
    assert numpy.isclose(1.54, val1, val2)

    val3 = dict_.value_by_unordered_key(DCT4, ('C', 'H'), fill_val=1.15)
    val4 = dict_.value_by_unordered_key(DCT4, ('C', 'N'), fill_val=1.15)
    assert numpy.isclose(1.54, val3)
    assert numpy.isclose(1.15, val4)

    # Multilevel dct reads
    val1 = dict_.values_in_multilevel_dct(
        GDCT3, 'key1', 'subkey2', fill_val='fill')
    val2 = dict_.values_in_multilevel_dct(
        GDCT3, 'key2', 'subkey1', fill_val='fill')
    val3 = dict_.values_in_multilevel_dct(
        GDCT3, 'key1', 'subkey3', fill_val='fill')
    val4 = dict_.values_in_multilevel_dct(
        GDCT3, 'key3', 'subkey1', fill_val='fill')

    assert val1 == 'subval1-2'
    assert val2 == 'subval2-1'
    assert val3 == 'fill'
    assert val4 == 'fill'


def test__sort():
    """ test dict_.keys_sorted_by_value
    """
    assert dict_.keys_sorted_by_value(GDCT5) == ('a', 'b', 'c', 'd', 'e')


def test__filter():
    """ test dict_.filter_keys
    """

    ref_filt_dct1 = {
        'key2': 'val2',
        'key3': 'val3'
    }
    ref_filt_dct2 = {
        'key1': {
            'subkey1': 'subval1-1',
            'subkey2': 'subval1-2'
        },
        'key2': {
            'subkey1': 'subval2-1',
            'subkey2': 'subval2-2'
        }
    }

    filt_dct1 = dict_.filter_keys(GDCT1, GDCT2)
    filt_dct2 = dict_.filter_keys(GDCT3, GDCT4)
    assert ref_filt_dct1 == filt_dct1
    assert ref_filt_dct2 == filt_dct2


def test__invert():
    """ test
    """

    ref_inv_dct = {
        'val1': 'key1',
        'val2': 'key2',
        'val3': 'key3'
    }

    inv_dct = dict_.invert(GDCT1)
    assert ref_inv_dct == inv_dct
