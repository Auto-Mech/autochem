""" test ring functionality in automol.graph
"""
from automol import graph


def test__rings():
    """ test graph.rings
    """
    c5h5n5o_cgr = (
        {0: ('C', 1, None), 1: ('C', 0, None), 2: ('C', 0, None),
         3: ('C', 0, None), 4: ('C', 0, None), 5: ('N', 2, None),
         6: ('N', 0, None), 7: ('N', 0, None), 8: ('N', 0, None),
         9: ('N', 1, None), 10: ('O', 1, None)},
        {frozenset({10, 4}): (1, None), frozenset({8, 2}): (1, None),
         frozenset({0, 6}): (1, None), frozenset({9, 3}): (1, None),
         frozenset({1, 2}): (1, None), frozenset({3, 7}): (1, None),
         frozenset({2, 5}): (1, None), frozenset({1, 6}): (1, None),
         frozenset({0, 7}): (1, None), frozenset({9, 4}): (1, None),
         frozenset({1, 3}): (1, None), frozenset({8, 4}): (1, None)})

    assert graph.rings(c5h5n5o_cgr) == (
        ({0: ('C', 1, None), 1: ('C', 0, None), 3: ('C', 0, None),
          6: ('N', 0, None), 7: ('N', 0, None)},
         {frozenset({0, 6}): (1, None), frozenset({3, 7}): (1, None),
          frozenset({0, 7}): (1, None), frozenset({1, 6}): (1, None),
          frozenset({1, 3}): (1, None)}),
        ({1: ('C', 0, None), 2: ('C', 0, None), 3: ('C', 0, None),
          4: ('C', 0, None), 8: ('N', 0, None), 9: ('N', 1, None)},
         {frozenset({8, 2}): (1, None), frozenset({9, 3}): (1, None),
          frozenset({1, 2}): (1, None), frozenset({9, 4}): (1, None),
          frozenset({1, 3}): (1, None), frozenset({8, 4}): (1, None)})
    )


def test__ring_systems():
    """ test automol.graph.ring_systems
    """
    # molecule:
    # InChI=1S/C13H20/c1-2-9-5-8(1)6-10(9)7-13-11-3-4-12(11)13/h8-13H,1-7H2
    gra = ({0: ('C', 1, None), 1: ('C', 2, None), 2: ('C', 1, None),
            5: ('C', 2, None), 8: ('C', 1, None), 11: ('C', 2, None),
            13: ('C', 2, None), 15: ('C', 2, None), 19: ('C', 1, None),
            20: ('C', 1, None), 22: ('C', 1, None), 27: ('C', 2, None),
            29: ('C', 2, None)},
           {frozenset({13, 15}): (1, None), frozenset({0, 1}): (1, None),
            frozenset({0, 2}): (1, None), frozenset({2, 5}): (1, None),
            frozenset({2, 15}): (1, None), frozenset({19, 22}): (1, None),
            frozenset({8, 1}): (1, None), frozenset({19, 29}): (1, None),
            frozenset({0, 11}): (1, None), frozenset({11, 22}): (1, None),
            frozenset({8, 13}): (1, None), frozenset({8, 5}): (1, None),
            frozenset({19, 20}): (1, None), frozenset({20, 22}): (1, None),
            frozenset({27, 29}): (1, None), frozenset({27, 20}): (1, None)})
    rsys = graph.ring_systems(gra)
    assert len(rsys) == 2

    rsy_rngs = list(map(graph.rings, rsys))
    assert tuple(map(len, rsy_rngs)) == (2, 2)


if __name__ == '__main__':
    test__rings()
    test__ring_systems()
