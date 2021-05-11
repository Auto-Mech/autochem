""" test automol.formula not utilized by other tests
"""

from automol import formula


FORMULA = {'C': 2, 'H': 6, 'O': 1}


def test__formula():
    """ test automol.formula.
    """

    # Test counters
    ref_n_elec = 26
    n_elec = formula.electron_count(FORMULA)

    assert n_elec == ref_n_elec

    ref_n_carbon = 2
    ref_n_hydrogen = 6
    n_carbon = formula.element_count(FORMULA, 'C')
    n_hydrogen = formula.element_count(FORMULA, 'H')

    assert n_carbon == ref_n_carbon
    assert n_hydrogen == ref_n_hydrogen

    # Test string
    fml_str1 = formula.string(FORMULA)
    fml_str2 = formula.string2(FORMULA)
    print(fml_str1)
    print(fml_str2)
