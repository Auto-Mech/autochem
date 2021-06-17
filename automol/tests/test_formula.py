""" test automol.formula not utilized by other tests
"""

from automol import formula


FORMULA = {'C': 2, 'H': 6, 'O': 1}


def test__formula():
    """ test automol.formula.
    """

    assert formula.electron_count(FORMULA) == 26

    assert formula.element_count(FORMULA, 'C') == 2
    assert formula.element_count(FORMULA, 'H') == 6

    assert formula.string(FORMULA) == 'C2H6O'
    assert formula.string2(FORMULA) == 'C2H6O1'

    assert formula.from_string('C2H6O') == {'C': 2, 'H': 6, 'O': 1}
    assert formula.from_string('C2H6O1') == {'C': 2, 'H': 6, 'O': 1}
    assert formula.from_string('C2O1H6') == {'C': 2, 'O': 1, 'H': 6}
    assert formula.from_string('C2OH6') == {'C': 2, 'O': 1, 'H': 6}
