""" tests for random utility functions
"""

import numpy
from automol import util


def test__separate_negatives():
    """test util.separate_negatives"""

    alst1 = (1.0, 2.0, 3.0)
    alst2 = (-1.0, -2.0, -3.0)
    alst3 = (1.0, 2.0, -3.0)

    nlst1, plst1 = util.separate_negatives(alst1)
    nlst2, plst2 = util.separate_negatives(alst2)
    nlst3, plst3 = util.separate_negatives(alst3)
    assert nlst1 == () and plst1 == (1.0, 2.0, 3.0)
    assert nlst2 == (-1.0, -2.0, -3.0) and plst2 == ()
    assert nlst3 == (-3.0,) and plst3 == (1.0, 2.0)


def test__separate_value_similar_to():
    """test util.value_similar_to"""

    vals = (1.111, 2.222, 3.333)
    testval1 = 1.121
    testval2 = 2.232
    testval3 = 3.533

    assert util.value_similar_to(testval1, vals, 0.1)
    assert util.value_similar_to(testval2, vals, 0.1)
    assert not util.value_similar_to(testval3, vals, 0.1)


def test__scale_iterable():
    """test util.scale_iterable"""

    ref_scale_iter = (5.0, 10.0, 15.0)

    scale = 5.0
    iter1 = [1.0, 2.0, 3.0]
    iter2 = (1.0, 2.0, 3.0)

    scale_iter1 = util.scale_iterable(iter1, scale)
    scale_iter2 = util.scale_iterable(iter2, scale)

    assert numpy.allclose(ref_scale_iter, scale_iter1)
    assert numpy.allclose(ref_scale_iter, scale_iter2)
    assert isinstance(scale_iter1, list)
    assert isinstance(scale_iter2, tuple)


def test__sort_by_lst():
    """test util.sort_by_lst()"""

    lst = ("q", "a", "e", "x", "f", "t")
    ref_lst = ("x", "a", "q", "e")

    assert ("x", "a", "q", "e", "f", "t") == util.sort_by_list(
        lst, ref_lst, include_missing=True
    )
    assert ("x", "a", "q", "e") == util.sort_by_list(
        lst, ref_lst, include_missing=False
    )
