""" test automol.formula not utilized by other tests
"""

from automol import form


FORMULA = {'C': 2, 'H': 6, 'O': 1}


def test__formula():
    """ test automol.formula.electron_count
        test automol.formula.element_count
        test automol.formula.string
        test automol.formula.string2
        test automol.formula.from_string
    """

    assert form.electron_count(FORMULA) == 26

    assert form.element_count(FORMULA, 'C') == 2
    assert form.element_count(FORMULA, 'H') == 6

    assert form.string(FORMULA) == 'C2H6O'
    assert form.string2(FORMULA) == 'C2H6O1'

    assert form.from_string('C2H6O') == {'C': 2, 'H': 6, 'O': 1}
    assert form.from_string('C2H6O1') == {'C': 2, 'H': 6, 'O': 1}
    assert form.from_string('C2O1H6') == {'C': 2, 'O': 1, 'H': 6}
    assert form.from_string('C2OH6') == {'C': 2, 'O': 1, 'H': 6}


def test__sort():
    """ test automol.formula.sorted_symbols
    """

    ref_sort_symbs = ('C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'O', 'N')

    symbs = ('C', 'N', 'H', 'H', 'C', 'O', 'H', 'H', 'H', 'H', 'H')
    sort_symbs = form.sorted_symbols(
        symbs, symbs_first=('C', 'H'), symbs_last=('N',))

    assert ref_sort_symbs == sort_symbs
