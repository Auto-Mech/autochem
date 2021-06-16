""" test automol.formula not utilized by other tests
"""

from automol import formula


FORMULA = {'C': 2, 'H': 6, 'O': 1}


def test__formula():
    """ test automol.formula.
    """

    assert formula.electron_count(FORMULA) == 26

    assert formula.element_count(FORMULA, 'C') == 2
    assert  formula.element_count(FORMULA, 'H') == 6

    assert formula.string(FORMULA) == 'C2H6O'
    assert formula.string2(FORMULA) == 'C2H6O1'
