"""test automol.graph"""

import itertools

import automol
import numpy
import pytest
from automol import graph

# Vinyl radical with E/Z stereo
C3H5_SGR = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, None),
        3: ("H", 0, None),
        4: ("H", 0, None),
        5: ("H", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
    },
    {
        frozenset({0, 2}): (1, False),
        frozenset({0, 3}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({1, 4}): (1, None),
        frozenset({1, 5}): (1, None),
        frozenset({1, 6}): (1, None),
        frozenset({2, 7}): (1, None),
    },
)


C8H13O_CGR = (
    {
        0: ("C", 3, None),
        1: ("C", 2, None),
        2: ("C", 3, None),
        3: ("C", 1, None),
        4: ("C", 1, None),
        5: ("C", 1, None),
        6: ("C", 1, None),
        7: ("C", 1, None),
        8: ("O", 0, None),
    },
    {
        frozenset({1, 4}): (1, None),
        frozenset({4, 6}): (1, None),
        frozenset({0, 3}): (1, None),
        frozenset({2, 6}): (1, None),
        frozenset({6, 7}): (1, None),
        frozenset({8, 7}): (1, None),
        frozenset({3, 5}): (1, None),
        frozenset({5, 7}): (1, None),
    },
)
C8H13O_RGR = (
    {
        0: ("C", 3, None),
        1: ("C", 2, None),
        2: ("C", 3, None),
        3: ("C", 1, None),
        4: ("C", 1, None),
        5: ("C", 1, None),
        6: ("C", 1, None),
        7: ("C", 1, None),
        8: ("O", 0, None),
    },
    {
        frozenset({1, 4}): (2, None),
        frozenset({4, 6}): (1, None),
        frozenset({0, 3}): (1, None),
        frozenset({2, 6}): (1, None),
        frozenset({6, 7}): (1, None),
        frozenset({8, 7}): (1, None),
        frozenset({3, 5}): (2, None),
        frozenset({5, 7}): (1, None),
    },
)
C8H13O_SGR = (
    {
        0: ("C", 3, None),
        1: ("C", 2, None),
        2: ("C", 3, None),
        3: ("C", 1, None),
        4: ("C", 1, None),
        5: ("C", 1, None),
        6: ("C", 1, False),
        7: ("C", 1, False),
        8: ("O", 0, None),
    },
    {
        frozenset({1, 4}): (1, None),
        frozenset({4, 6}): (1, None),
        frozenset({0, 3}): (1, None),
        frozenset({2, 6}): (1, None),
        frozenset({6, 7}): (1, None),
        frozenset({8, 7}): (1, None),
        frozenset({3, 5}): (1, False),
        frozenset({5, 7}): (1, None),
    },
)


C3H3_CGR = (
    {0: ("C", 1, None), 1: ("C", 1, None), 2: ("C", 1, None)},
    {
        frozenset({0, 1}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({2, 0}): (1, None),
    },
)
C3H3_RGRS = (
    (
        {0: ("C", 1, None), 1: ("C", 1, None), 2: ("C", 1, None)},
        {
            frozenset({0, 1}): (1, None),
            frozenset({1, 2}): (2, None),
            frozenset({2, 0}): (1, None),
        },
    ),
    (
        {0: ("C", 1, None), 1: ("C", 1, None), 2: ("C", 1, None)},
        {
            frozenset({0, 1}): (1, None),
            frozenset({1, 2}): (1, None),
            frozenset({2, 0}): (2, None),
        },
    ),
    (
        {0: ("C", 1, None), 1: ("C", 1, None), 2: ("C", 1, None)},
        {
            frozenset({0, 1}): (2, None),
            frozenset({1, 2}): (1, None),
            frozenset({2, 0}): (1, None),
        },
    ),
)

C2_CGR = ({0: ("C", 0, None), 1: ("C", 0, None)}, {frozenset({0, 1}): (1, None)})
C2_RGRS = (
    ({0: ("C", 0, None), 1: ("C", 0, None)}, {frozenset({0, 1}): (1, None)}),
    ({0: ("C", 0, None), 1: ("C", 0, None)}, {frozenset({0, 1}): (2, None)}),
    ({0: ("C", 0, None), 1: ("C", 0, None)}, {frozenset({0, 1}): (3, None)}),
)

CH2FH2H_CGR_IMP = (
    {0: ("F", 0, None), 1: ("C", 2, None), 2: ("H", 1, None), 3: ("H", 0, None)},
    {frozenset({0, 1}): (1, None)},
)
CH2FH2H_CGR_EXP = (
    {
        0: ("F", 0, None),
        1: ("C", 0, None),
        2: ("H", 0, None),
        3: ("H", 0, None),
        4: ("H", 0, None),
        5: ("H", 0, None),
        6: ("H", 0, None),
    },
    {
        frozenset({0, 1}): (1, None),
        frozenset({1, 4}): (1, None),
        frozenset({1, 5}): (1, None),
        frozenset({2, 6}): (1, None),
    },
)

C2H2CL2F2_CGR = (
    {
        0: ("C", 1, None),
        1: ("C", 1, None),
        2: ("F", 0, None),
        3: ("Cl", 0, None),
        4: ("F", 0, None),
        5: ("Cl", 0, None),
    },
    {
        frozenset({0, 1}): (1, None),
        frozenset({0, 2}): (1, None),
        frozenset({0, 3}): (1, None),
        frozenset({1, 4}): (1, None),
        frozenset({1, 5}): (1, None),
    },
)
C2H2CL2F2_SGRS = tuple(
    sorted(
        [
            (
                {
                    0: ("C", 1, False),
                    1: ("C", 1, False),
                    2: ("F", 0, None),
                    3: ("Cl", 0, None),
                    4: ("F", 0, None),
                    5: ("Cl", 0, None),
                },
                {
                    frozenset({0, 1}): (1, None),
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({1, 4}): (1, None),
                    frozenset({1, 5}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 1, False),
                    1: ("C", 1, True),
                    2: ("F", 0, None),
                    3: ("Cl", 0, None),
                    4: ("F", 0, None),
                    5: ("Cl", 0, None),
                },
                {
                    frozenset({0, 1}): (1, None),
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({1, 4}): (1, None),
                    frozenset({1, 5}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 1, True),
                    1: ("C", 1, True),
                    2: ("F", 0, None),
                    3: ("Cl", 0, None),
                    4: ("F", 0, None),
                    5: ("Cl", 0, None),
                },
                {
                    frozenset({0, 1}): (1, None),
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({1, 4}): (1, None),
                    frozenset({1, 5}): (1, None),
                },
            ),
        ],
        key=graph.frozen,
    )
)

C3H3CL2F3_CGR = (
    {
        0: ("C", 1, None),
        1: ("C", 1, None),
        2: ("C", 1, None),
        3: ("Cl", 0, None),
        4: ("Cl", 0, None),
        5: ("F", 0, None),
        6: ("F", 0, None),
        7: ("F", 0, None),
    },
    {
        frozenset({0, 1}): (1, None),
        frozenset({0, 2}): (1, None),
        frozenset({0, 5}): (1, None),
        frozenset({2, 4}): (1, None),
        frozenset({1, 3}): (1, None),
        frozenset({1, 6}): (1, None),
        frozenset({2, 7}): (1, None),
    },
)
C3H3CL2F3_SGRS = tuple(
    sorted(
        [
            (
                {
                    0: ("C", 1, None),
                    1: ("C", 1, False),
                    2: ("C", 1, False),
                    3: ("Cl", 0, None),
                    4: ("Cl", 0, None),
                    5: ("F", 0, None),
                    6: ("F", 0, None),
                    7: ("F", 0, None),
                },
                {
                    frozenset({0, 1}): (1, None),
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 5}): (1, None),
                    frozenset({2, 4}): (1, None),
                    frozenset({1, 3}): (1, None),
                    frozenset({1, 6}): (1, None),
                    frozenset({2, 7}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 1, None),
                    1: ("C", 1, True),
                    2: ("C", 1, True),
                    3: ("Cl", 0, None),
                    4: ("Cl", 0, None),
                    5: ("F", 0, None),
                    6: ("F", 0, None),
                    7: ("F", 0, None),
                },
                {
                    frozenset({0, 1}): (1, None),
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 5}): (1, None),
                    frozenset({2, 4}): (1, None),
                    frozenset({1, 3}): (1, None),
                    frozenset({1, 6}): (1, None),
                    frozenset({2, 7}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 1, False),
                    1: ("C", 1, False),
                    2: ("C", 1, True),
                    3: ("Cl", 0, None),
                    4: ("Cl", 0, None),
                    5: ("F", 0, None),
                    6: ("F", 0, None),
                    7: ("F", 0, None),
                },
                {
                    frozenset({0, 1}): (1, None),
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 5}): (1, None),
                    frozenset({2, 4}): (1, None),
                    frozenset({1, 3}): (1, None),
                    frozenset({1, 6}): (1, None),
                    frozenset({2, 7}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 1, True),
                    1: ("C", 1, False),
                    2: ("C", 1, True),
                    3: ("Cl", 0, None),
                    4: ("Cl", 0, None),
                    5: ("F", 0, None),
                    6: ("F", 0, None),
                    7: ("F", 0, None),
                },
                {
                    frozenset({0, 1}): (1, None),
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 5}): (1, None),
                    frozenset({2, 4}): (1, None),
                    frozenset({1, 3}): (1, None),
                    frozenset({1, 6}): (1, None),
                    frozenset({2, 7}): (1, None),
                },
            ),
        ],
        key=graph.frozen,
    )
)

C3H5N3_CGR = (
    {
        0: ("C", 1, None),
        1: ("C", 1, None),
        2: ("C", 0, None),
        3: ("N", 1, None),
        4: ("N", 1, None),
        5: ("N", 1, None),
    },
    {
        frozenset({1, 4}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({0, 3}): (1, None),
        frozenset({0, 2}): (1, None),
        frozenset({2, 5}): (1, None),
    },
)
C3H5N3_SGRS = tuple(
    sorted(
        [
            (
                {
                    0: ("C", 1, None),
                    1: ("C", 1, None),
                    2: ("C", 0, None),
                    3: ("N", 1, None),
                    4: ("N", 1, None),
                    5: ("N", 1, None),
                },
                {
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 3}): (1, False),
                    frozenset({1, 2}): (1, None),
                    frozenset({1, 4}): (1, False),
                    frozenset({2, 5}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 1, None),
                    1: ("C", 1, None),
                    2: ("C", 0, None),
                    3: ("N", 1, None),
                    4: ("N", 1, None),
                    5: ("N", 1, None),
                },
                {
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 3}): (1, False),
                    frozenset({1, 2}): (1, None),
                    frozenset({1, 4}): (1, True),
                    frozenset({2, 5}): (1, False),
                },
            ),
            (
                {
                    0: ("C", 1, None),
                    1: ("C", 1, None),
                    2: ("C", 0, None),
                    3: ("N", 1, None),
                    4: ("N", 1, None),
                    5: ("N", 1, None),
                },
                {
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 3}): (1, False),
                    frozenset({1, 2}): (1, None),
                    frozenset({1, 4}): (1, True),
                    frozenset({2, 5}): (1, True),
                },
            ),
            (
                {
                    0: ("C", 1, None),
                    1: ("C", 1, None),
                    2: ("C", 0, None),
                    3: ("N", 1, None),
                    4: ("N", 1, None),
                    5: ("N", 1, None),
                },
                {
                    frozenset({0, 2}): (1, None),
                    frozenset({0, 3}): (1, True),
                    frozenset({1, 2}): (1, None),
                    frozenset({1, 4}): (1, True),
                    frozenset({2, 5}): (1, None),
                },
            ),
        ],
        key=graph.frozen,
    )
)

C8H13O_SGRS = tuple(
    sorted(
        [
            (
                {
                    0: ("C", 3, None),
                    1: ("C", 2, None),
                    2: ("C", 3, None),
                    3: ("C", 1, None),
                    4: ("C", 1, None),
                    5: ("C", 1, None),
                    6: ("C", 1, False),
                    7: ("C", 1, False),
                    8: ("O", 0, None),
                },
                {
                    frozenset({1, 4}): (1, None),
                    frozenset({4, 6}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({2, 6}): (1, None),
                    frozenset({6, 7}): (1, None),
                    frozenset({8, 7}): (1, None),
                    frozenset({3, 5}): (1, False),
                    frozenset({5, 7}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 3, None),
                    1: ("C", 2, None),
                    2: ("C", 3, None),
                    3: ("C", 1, None),
                    4: ("C", 1, None),
                    5: ("C", 1, None),
                    6: ("C", 1, False),
                    7: ("C", 1, False),
                    8: ("O", 0, None),
                },
                {
                    frozenset({1, 4}): (1, None),
                    frozenset({4, 6}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({2, 6}): (1, None),
                    frozenset({6, 7}): (1, None),
                    frozenset({8, 7}): (1, None),
                    frozenset({3, 5}): (1, True),
                    frozenset({5, 7}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 3, None),
                    1: ("C", 2, None),
                    2: ("C", 3, None),
                    3: ("C", 1, None),
                    4: ("C", 1, None),
                    5: ("C", 1, None),
                    6: ("C", 1, False),
                    7: ("C", 1, True),
                    8: ("O", 0, None),
                },
                {
                    frozenset({1, 4}): (1, None),
                    frozenset({4, 6}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({2, 6}): (1, None),
                    frozenset({6, 7}): (1, None),
                    frozenset({8, 7}): (1, None),
                    frozenset({3, 5}): (1, False),
                    frozenset({5, 7}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 3, None),
                    1: ("C", 2, None),
                    2: ("C", 3, None),
                    3: ("C", 1, None),
                    4: ("C", 1, None),
                    5: ("C", 1, None),
                    6: ("C", 1, False),
                    7: ("C", 1, True),
                    8: ("O", 0, None),
                },
                {
                    frozenset({1, 4}): (1, None),
                    frozenset({4, 6}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({2, 6}): (1, None),
                    frozenset({6, 7}): (1, None),
                    frozenset({8, 7}): (1, None),
                    frozenset({3, 5}): (1, True),
                    frozenset({5, 7}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 3, None),
                    1: ("C", 2, None),
                    2: ("C", 3, None),
                    3: ("C", 1, None),
                    4: ("C", 1, None),
                    5: ("C", 1, None),
                    6: ("C", 1, True),
                    7: ("C", 1, False),
                    8: ("O", 0, None),
                },
                {
                    frozenset({1, 4}): (1, None),
                    frozenset({4, 6}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({2, 6}): (1, None),
                    frozenset({6, 7}): (1, None),
                    frozenset({8, 7}): (1, None),
                    frozenset({3, 5}): (1, False),
                    frozenset({5, 7}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 3, None),
                    1: ("C", 2, None),
                    2: ("C", 3, None),
                    3: ("C", 1, None),
                    4: ("C", 1, None),
                    5: ("C", 1, None),
                    6: ("C", 1, True),
                    7: ("C", 1, False),
                    8: ("O", 0, None),
                },
                {
                    frozenset({1, 4}): (1, None),
                    frozenset({4, 6}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({2, 6}): (1, None),
                    frozenset({6, 7}): (1, None),
                    frozenset({8, 7}): (1, None),
                    frozenset({3, 5}): (1, True),
                    frozenset({5, 7}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 3, None),
                    1: ("C", 2, None),
                    2: ("C", 3, None),
                    3: ("C", 1, None),
                    4: ("C", 1, None),
                    5: ("C", 1, None),
                    6: ("C", 1, True),
                    7: ("C", 1, True),
                    8: ("O", 0, None),
                },
                {
                    frozenset({1, 4}): (1, None),
                    frozenset({4, 6}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({2, 6}): (1, None),
                    frozenset({6, 7}): (1, None),
                    frozenset({8, 7}): (1, None),
                    frozenset({3, 5}): (1, False),
                    frozenset({5, 7}): (1, None),
                },
            ),
            (
                {
                    0: ("C", 3, None),
                    1: ("C", 2, None),
                    2: ("C", 3, None),
                    3: ("C", 1, None),
                    4: ("C", 1, None),
                    5: ("C", 1, None),
                    6: ("C", 1, True),
                    7: ("C", 1, True),
                    8: ("O", 0, None),
                },
                {
                    frozenset({1, 4}): (1, None),
                    frozenset({4, 6}): (1, None),
                    frozenset({0, 3}): (1, None),
                    frozenset({2, 6}): (1, None),
                    frozenset({6, 7}): (1, None),
                    frozenset({8, 7}): (1, None),
                    frozenset({3, 5}): (1, True),
                    frozenset({5, 7}): (1, None),
                },
            ),
        ],
        key=graph.frozen,
    )
)

# FC=CC=CF + [OH] => FC=C[CH]C(O)F
C4H5F2O_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, None),
        3: ("C", 0, None),
        4: ("F", 0, None),
        5: ("F", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("O", 0, None),
        11: ("H", 0, None),
    },
    {
        frozenset({8, 2}): (1, None),
        frozenset({2, 10}): (0.1, None),
        frozenset({0, 6}): (1, None),
        frozenset({1, 7}): (1, None),
        frozenset({9, 3}): (1, None),
        frozenset({0, 1}): (1, None),
        frozenset({0, 2}): (1, True),
        frozenset({2, 4}): (1, None),
        frozenset({3, 5}): (1, None),
        frozenset({10, 11}): (1, None),
        frozenset({1, 3}): (1, False),
    },
)

# FC=C(C(O)F)C(O)F + [OH] => FC(O)[C](C(O)F)C(O)F
C4H5F3O2_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, False),
        3: ("C", 0, True),
        4: ("F", 0, None),
        5: ("F", 0, None),
        6: ("F", 0, None),
        7: ("O", 0, None),
        8: ("O", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
        13: ("H", 0, None),
        14: ("O", 0, None),
        15: ("H", 0, None),
    },
    {
        frozenset({12, 7}): (1, None),
        frozenset({2, 10}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({0, 1}): (1, True),
        frozenset({3, 6}): (1, None),
        frozenset({2, 7}): (1, None),
        frozenset({2, 5}): (1, None),
        frozenset({0, 4}): (1, None),
        frozenset({8, 3}): (1, None),
        frozenset({0, 14}): (0.1, None),
        frozenset({8, 13}): (1, None),
        frozenset({14, 15}): (1, None),
        frozenset({11, 3}): (1, None),
        frozenset({1, 3}): (1, None),
        frozenset({0, 9}): (1, None),
    },
)

# ISOBUTANE
C4H10_GRA = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, None),
        3: ("C", 0, None),
        4: ("H", 0, None),
        5: ("H", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
        13: ("H", 0, None),
    },
    {
        frozenset({0, 3}): (1, None),
        frozenset({0, 4}): (1, None),
        frozenset({0, 5}): (1, None),
        frozenset({0, 6}): (1, None),
        frozenset({1, 3}): (1, None),
        frozenset({1, 7}): (1, None),
        frozenset({8, 1}): (1, None),
        frozenset({1, 9}): (1, None),
        frozenset({2, 3}): (1, None),
        frozenset({2, 10}): (1, None),
        frozenset({2, 11}): (1, None),
        frozenset({2, 12}): (1, None),
        frozenset({3, 13}): (1, None),
    },
)

# C2H5 + H2
C2H7_CGR = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("H", 0, None),
        3: ("H", 0, None),
        4: ("H", 0, None),
        5: ("H", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
    },
    {
        frozenset({0, 1}): (1, None),
        frozenset({0, 2}): (1, None),
        frozenset({0, 3}): (1, None),
        frozenset({0, 4}): (1, None),
        frozenset({1, 6}): (1, None),
        frozenset({1, 7}): (1, None),
        frozenset({8, 5}): (1, None),
    },
)


def test__from_data():
    """test getters"""
    cgr = graph.from_data(
        atm_symb_dct=graph.atom_symbols(C8H13O_CGR),
        bnd_keys=graph.bond_keys(C8H13O_CGR),
        atm_imp_hyd_dct=graph.atom_implicit_hydrogens(C8H13O_CGR),
    )
    assert cgr == C8H13O_CGR

    rgr = graph.from_data(
        atm_symb_dct=graph.atom_symbols(C8H13O_RGR),
        bnd_keys=graph.bond_keys(C8H13O_RGR),
        atm_imp_hyd_dct=graph.atom_implicit_hydrogens(C8H13O_RGR),
        bnd_ord_dct=graph.bond_orders(C8H13O_RGR),
    )
    assert rgr == C8H13O_RGR

    sgr = graph.from_data(
        atm_symb_dct=graph.atom_symbols(C8H13O_SGR),
        bnd_keys=graph.bond_keys(C8H13O_SGR),
        atm_imp_hyd_dct=graph.atom_implicit_hydrogens(C8H13O_SGR),
        atm_ste_par_dct=graph.atom_stereo_parities(C8H13O_SGR),
        bnd_ste_par_dct=graph.bond_stereo_parities(C8H13O_SGR),
    )
    assert sgr == C8H13O_SGR


def test__setters():
    """test graph setters"""
    atm_symbs = numpy.array(list("CHON"))
    bnd_ords = numpy.arange(1, 4)
    atm_imp_hyds = numpy.arange(0, 4)
    pars = numpy.array([None, True, False])

    print("\nTesting setters for an ordinary molecular graph...")
    orig_gra = C8H13O_CGR

    atm_keys = graph.atom_keys(orig_gra)
    bnd_keys = graph.bond_keys(orig_gra)
    natms = len(atm_keys)
    nbnds = len(bnd_keys)

    # atom symbols
    orig_atm_symb_dct = graph.atom_symbols(orig_gra)
    atm_symb_dct = dict(zip(atm_keys, numpy.random.choice(atm_symbs, size=natms)))
    gra = graph.set_atom_symbols(orig_gra, atm_symb_dct)
    print(atm_symb_dct)
    assert atm_symb_dct == graph.atom_symbols(gra)
    print(graph.set_atom_symbols(gra, orig_atm_symb_dct))
    assert orig_gra == graph.set_atom_symbols(gra, orig_atm_symb_dct)

    # bond orders
    orig_bnd_ord_dct = graph.bond_orders(orig_gra)
    bnd_ord_dct = dict(zip(bnd_keys, numpy.random.choice(bnd_ords, size=nbnds)))
    gra = graph.set_bond_orders(orig_gra, bnd_ord_dct)
    print(bnd_ord_dct)
    assert bnd_ord_dct == graph.bond_orders(gra)
    assert orig_gra == graph.set_bond_orders(gra, orig_bnd_ord_dct)

    # atom implicit hydrogen valences
    orig_atm_imp_hyd_dct = graph.atom_implicit_hydrogens(orig_gra)
    atm_imp_hyd_dct = dict(zip(atm_keys, numpy.random.choice(atm_imp_hyds, size=natms)))
    gra = graph.set_atom_implicit_hydrogens(orig_gra, atm_imp_hyd_dct)
    print(atm_imp_hyd_dct)
    assert atm_imp_hyd_dct == graph.atom_implicit_hydrogens(gra)
    assert orig_gra == graph.set_atom_implicit_hydrogens(gra, orig_atm_imp_hyd_dct)

    # atom stereo parities
    orig_atm_par_dct = graph.atom_stereo_parities(orig_gra)
    atm_par_dct = dict(zip(atm_keys, numpy.random.choice(pars, size=natms)))
    gra = graph.set_atom_stereo_parities(orig_gra, atm_par_dct)
    print(atm_par_dct)
    assert atm_par_dct == graph.atom_stereo_parities(gra)
    assert orig_gra == graph.set_atom_stereo_parities(gra, orig_atm_par_dct)

    # bond stereo parities
    orig_bnd_par_dct = graph.bond_stereo_parities(orig_gra)
    bnd_par_dct = dict(zip(bnd_keys, numpy.random.choice(pars, size=nbnds)))
    gra = graph.set_bond_stereo_parities(orig_gra, bnd_par_dct)
    print(bnd_par_dct)
    assert bnd_par_dct == graph.bond_stereo_parities(gra)
    assert orig_gra == graph.set_bond_stereo_parities(gra, orig_bnd_par_dct)


def test__string():
    """test graph.string and graph.from_string"""
    for sgr in C8H13O_SGRS:
        assert sgr == graph.from_string(graph.string(sgr))


def test__without_bond_orders():
    """test graph.without_bond_orders"""
    assert C8H13O_CGR == graph.without_pi_bonds(C8H13O_RGR)


def test__without_stereo_parities():
    """test graph.without_stereo_parities"""
    assert C8H13O_CGR == graph.without_stereo(C8H13O_SGR)


def test__electron_count():
    """test graph.electron_count"""
    assert graph.electron_count(C8H13O_CGR) == 69


def test__atom_count():
    """test graph.atom_count"""
    assert graph.atom_count(C8H13O_CGR) == 22
    assert graph.atom_count(C8H13O_CGR, symb="C") == 8
    assert graph.atom_count(C8H13O_CGR, symb="O") == 1
    assert graph.atom_count(C8H13O_CGR, heavy_only=True) == 9
    assert graph.atom_count(C8H13O_CGR, symb="H") == 13
    assert graph.atom_count(C8H13O_CGR, symb="H", heavy_only=True) == 0
    assert graph.atom_count(C8H13O_CGR, with_implicit=False) == 9

    cgr = (
        {
            0: ("C", 0, None),
            1: ("C", 0, None),
            2: ("O", 0, None),
            3: ("H", 0, None),
            4: ("H", 0, None),
            5: ("H", 0, None),
            6: ("X", 0, None),
        },
        {
            frozenset({1, 4}): (1, None),
            frozenset({1, 2}): (1, None),
            frozenset({0, 3}): (1, None),
            frozenset({0, 1}): (1, None),
            frozenset({2, 5}): (1, None),
        },
    )
    assert graph.atom_count(cgr) == 6
    assert graph.atom_count(cgr, dummy=True) == 7


def test__atoms_neighbor_atom_keys():
    """test graph.atoms_neighbor_atom_keys"""
    assert graph.atoms_neighbor_atom_keys(C8H13O_CGR) == {
        0: frozenset({3}),
        1: frozenset({4}),
        2: frozenset({6}),
        3: frozenset({0, 5}),
        4: frozenset({1, 6}),
        5: frozenset({3, 7}),
        6: frozenset({2, 4, 7}),
        7: frozenset({8, 5, 6}),
        8: frozenset({7}),
    }


def test__atoms_bond_keys():
    """test graph.atoms_neighbor_atom_keys"""
    assert graph.atoms_bond_keys(C8H13O_CGR) == {
        0: frozenset({frozenset({0, 3})}),
        1: frozenset({frozenset({1, 4})}),
        2: frozenset({frozenset({2, 6})}),
        3: frozenset({frozenset({3, 5}), frozenset({0, 3})}),
        4: frozenset({frozenset({1, 4}), frozenset({4, 6})}),
        5: frozenset({frozenset({3, 5}), frozenset({5, 7})}),
        6: frozenset({frozenset({6, 7}), frozenset({4, 6}), frozenset({2, 6})}),
        7: frozenset({frozenset({6, 7}), frozenset({5, 7}), frozenset({8, 7})}),
        8: frozenset({frozenset({8, 7})}),
    }


# # bond properties
def test__bonds_neighbor_atom_keys():
    """test graph.bonds_neighbor_atom_keys"""

    assert graph.bonds_neighbor_atom_keys(C8H13O_CGR, group=False) == {
        frozenset({1, 4}): frozenset({6}),
        frozenset({4, 6}): frozenset({1, 2, 7}),
        frozenset({2, 6}): frozenset({4, 7}),
        frozenset({0, 3}): frozenset({5}),
        frozenset({6, 7}): frozenset({8, 2, 4, 5}),
        frozenset({8, 7}): frozenset({5, 6}),
        frozenset({3, 5}): frozenset({0, 7}),
        frozenset({5, 7}): frozenset({8, 3, 6}),
    }

    assert graph.bonds_neighbor_atom_keys(C8H13O_CGR, group=True) == {
        frozenset({1, 4}): (frozenset(), frozenset({6})),
        frozenset({4, 6}): (frozenset({1}), frozenset({2, 7})),
        frozenset({2, 6}): (frozenset(), frozenset({4, 7})),
        frozenset({0, 3}): (frozenset(), frozenset({5})),
        frozenset({6, 7}): (frozenset({2, 4}), frozenset({8, 5})),
        frozenset({7, 8}): (frozenset({5, 6}), frozenset()),
        frozenset({3, 5}): (frozenset({0}), frozenset({7})),
        frozenset({5, 7}): (frozenset({3}), frozenset({8, 6})),
    }


# # other properties
def test__branch():
    """test graph.branch"""
    # Using an atom key:
    assert graph.branch(C8H13O_CGR, 6, 4) == (
        {1: ("C", 2, None), 4: ("C", 1, None)},
        {frozenset({1, 4}): (1, None)},
    )
    assert graph.branch(C8H13O_CGR, 6, 4, keep_root=True) == (
        {1: ("C", 2, None), 4: ("C", 1, None), 6: ("C", 1, None)},
        {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None)},
    )
    # Using a bond key:
    bnd_key = frozenset({6, 4})
    assert graph.branch(C8H13O_CGR, 6, bnd_key) == (
        {1: ("C", 2, None), 4: ("C", 1, None)},
        {frozenset({1, 4}): (1, None)},
    )
    assert graph.branch(C8H13O_CGR, 6, bnd_key, keep_root=True) == (
        {1: ("C", 2, None), 4: ("C", 1, None), 6: ("C", 1, None)},
        {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None)},
    )


def test__connected_components():
    """test graph.connected_components"""
    (gra1, gra2), _ = graph.standard_keys_for_sequence([C3H3_CGR, C2_CGR])

    cmp_gras = graph.connected_components(graph.union(gra1, gra2))
    assert cmp_gras in [(gra1, gra2), (gra2, gra1)]


def test__subgraph():
    """test graph.subgraph"""
    assert graph.subgraph(C3H3_CGR, (1, 2)) == (
        {1: ("C", 1, None), 2: ("C", 1, None)},
        {frozenset({1, 2}): (1, None)},
    )


def test__bond_induced_subgraph():
    """test graph.bond_induced_subgraph"""
    assert graph.bond_induced_subgraph(
        C3H3_CGR, [frozenset({0, 1}), frozenset({1, 2})]
    ) == (
        {0: ("C", 1, None), 1: ("C", 1, None), 2: ("C", 1, None)},
        {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)},
    )


# # transformations
def test__relabel():
    """test graph.relabel"""
    assert graph.relabel(C3H3_CGR, {0: 10, 1: 11, 2: 12}) == (
        {10: ("C", 1, None), 11: ("C", 1, None), 12: ("C", 1, None)},
        {
            frozenset({10, 11}): (1, None),
            frozenset({11, 12}): (1, None),
            frozenset({12, 10}): (1, None),
        },
    )


def test__remove_atoms():
    """test graph.remove_atoms"""
    assert graph.remove_atoms(C3H3_CGR, (0,)) == (
        {1: ("C", 1, None), 2: ("C", 1, None)},
        {frozenset({1, 2}): (1, None)},
    )


def test__remove_bonds():
    """test graph.remove_bonds"""
    assert graph.remove_bonds(C3H3_CGR, [frozenset({1, 2})]) == (
        {0: ("C", 1, None), 1: ("C", 1, None), 2: ("C", 1, None)},
        {frozenset({0, 1}): (1, None), frozenset({2, 0}): (1, None)},
    )


# implicit/explicit hydrogen functions
# # atom properties
def test__atom_hydrogen_keys():
    """test graph.atom_hydrogen_keys"""
    assert graph.atom_nonbackbone_hydrogen_keys(CH2FH2H_CGR_EXP) == {
        0: frozenset(),
        1: frozenset({4, 5}),
        2: frozenset({6}),
        3: frozenset(),
        4: frozenset(),
        5: frozenset(),
        6: frozenset(),
    }


# # other properties
def test__backbone_keys():
    """test graph.backbone_keys"""
    assert graph.backbone_keys(CH2FH2H_CGR_EXP) == frozenset({0, 1, 2, 3})


def test__hydrogen_keys():
    """test graph.hydrogen_keys"""
    assert graph.nonbackbone_hydrogen_keys(CH2FH2H_CGR_EXP) == frozenset({4, 5, 6})


def test__explicit():
    """test graph.explicit"""
    assert CH2FH2H_CGR_EXP == graph.explicit(CH2FH2H_CGR_IMP)


def test__implicit():
    """test graph.implicit"""
    assert CH2FH2H_CGR_IMP == graph.implicit(graph.explicit(CH2FH2H_CGR_IMP))


# # comparisons
def test__isomorphic():
    """test graph.isomorphic"""
    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        assert graph.isomorphic(cgr, cgr_pmt)

    # Test backbone_only option, comparing implicit and explicit graphs
    assert graph.isomorphic(CH2FH2H_CGR_IMP, CH2FH2H_CGR_EXP, backbone_only=True)


def test__isomorphism():
    """test graph.isomorphism"""
    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        assert graph.isomorphism(cgr, cgr_pmt) == pmt_dct


def test__unique():
    """test graph.unique"""
    assert graph.unique(C3H3_RGRS) == C3H3_RGRS[:1]


# chemistry library
def test__atomic_valences():
    """test graph.atomic_valences"""
    assert graph.atomic_valences(C8H13O_CGR) == {
        0: 4,
        1: 4,
        2: 4,
        3: 4,
        4: 4,
        5: 4,
        6: 4,
        7: 4,
        8: 2,
    }


def test__atom_lone_pairs():
    """test graph.atom_lone_pairs"""
    assert graph.atom_lone_pairs(C8H13O_CGR) == {
        0: 0,
        1: 0,
        2: 0,
        3: 0,
        4: 0,
        5: 0,
        6: 0,
        7: 0,
        8: 2,
    }


def test__atom_bond_counts():
    """test graph.atom_bond_counts"""
    assert graph.atom_bond_counts(C8H13O_CGR) == {
        0: 4,
        1: 3,
        2: 4,
        3: 3,
        4: 3,
        5: 3,
        6: 4,
        7: 4,
        8: 1,
    }


def test__atom_unsaturations():
    """test graph.atom_unsaturations"""
    assert graph.atom_unpaired_electrons(C8H13O_CGR) == {
        0: 0,
        1: 1,
        2: 0,
        3: 1,
        4: 1,
        5: 1,
        6: 0,
        7: 0,
        8: 1,
    }


def test__unsaturated_atom_keys():
    """test graph.unsaturated_atom_keys"""
    assert graph.unsaturated_atom_keys(C8H13O_CGR) == frozenset({1, 3, 4, 5, 8})


def test__maximum_spin_multiplicity():
    """test graph.maximum_spin_multiplicity"""
    assert graph.maximum_spin_multiplicity(C2_CGR) == 7


def test__possible_spin_multiplicities():
    """test graph.possible_spin_multiplicities"""
    assert graph.possible_spin_multiplicities(C2_CGR) == (1, 3, 5, 7)


# resonance graph library
# # atom properties
def test__atom_hybridizations():
    """test graph.atom_hybridizations"""
    assert graph.atom_hybridizations(C3H3_CGR) == {0: 2, 1: 2, 2: 2}
    assert graph.atom_hybridizations(C8H13O_CGR) == {
        0: 3,
        1: 2,
        2: 3,
        3: 2,
        4: 2,
        5: 2,
        6: 3,
        7: 3,
        8: 3,
    }

    cgr = (
        {
            0: ("C", 0, None),
            1: ("C", 0, None),
            2: ("O", 0, None),
            3: ("H", 0, None),
            4: ("H", 0, None),
            5: ("H", 0, None),
            6: ("X", 0, None),
        },
        {
            frozenset({1, 4}): (1, None),
            frozenset({1, 2}): (1, None),
            frozenset({0, 3}): (1, None),
            frozenset({0, 1}): (1, None),
            frozenset({2, 5}): (1, None),
        },
    )
    assert graph.atom_hybridizations(cgr) == {0: 2, 1: 2, 2: 3, 3: 0, 4: 0, 5: 0, 6: -1}


def test__atom_centered_cumulene_keys():
    """test graph.atom_centered_cumulene_keys"""
    cgr = (
        {
            0: ("C", 1, None),
            1: ("C", 2, None),
            2: ("C", 0, None),
            3: ("C", 0, None),
            4: ("C", 1, None),
            5: ("C", 0, None),
            6: ("C", 0, None),
        },
        {
            frozenset({4, 6}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({2, 4}): (1, None),
            frozenset({5, 6}): (1, None),
            frozenset({3, 5}): (1, None),
            frozenset({1, 3}): (1, None),
        },
    )
    assert graph.atom_centered_cumulene_keys(cgr) == frozenset({(frozenset({1, 4}), 5)})


def test__bond_centered_cumulene_keys():
    """test graph.bond_centered_cumulene_keys"""
    cgr = (
        {
            0: ("C", 1, None),
            1: ("C", 2, None),
            2: ("C", 0, None),
            3: ("C", 0, None),
            4: ("C", 1, None),
            5: ("C", 0, None),
        },
        {
            frozenset({4, 5}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({2, 4}): (1, None),
            frozenset({3, 5}): (1, None),
            frozenset({1, 3}): (1, None),
        },
    )
    assert graph.bond_centered_cumulene_keys(cgr) == frozenset(
        {(frozenset({1, 4}), frozenset({3, 5}))}
    )


def test__radical_atom_keys():
    """test graph.radical_atom_keys"""
    assert graph.radical_atom_keys(C3H3_CGR) == frozenset({0, 1, 2})
    assert graph.radical_atom_keys(C8H13O_CGR) == frozenset({8})


def test__sigma_radical_atom_keys():
    """test graph.sigma_radical_atom_keys"""
    # CCC#[C]
    gra = (
        {0: ("C", 3, None), 1: ("C", 0, None), 2: ("C", 2, None), 3: ("C", 0, None)},
        {
            frozenset({0, 2}): (1, None),
            frozenset({1, 3}): (1, None),
            frozenset({2, 3}): (1, None),
        },
    )
    assert graph.sigma_radical_atom_keys(gra) == frozenset({1})

    # [C]#CC(CC)(CCC#[C])CC#[C]
    gra = (
        {
            0: ("C", 0, None),
            1: ("C", 0, None),
            2: ("C", 3, None),
            3: ("C", 0, None),
            4: ("C", 0, None),
            5: ("C", 0, None),
            6: ("C", 2, None),
            7: ("C", 0, None),
            8: ("C", 2, None),
            9: ("C", 2, None),
            10: ("C", 2, None),
            11: ("C", 0, None),
        },
        {
            frozenset({8, 4}): (1, None),
            frozenset({3, 7}): (1, None),
            frozenset({2, 6}): (1, None),
            frozenset({0, 4}): (1, None),
            frozenset({8, 10}): (1, None),
            frozenset({9, 11}): (1, None),
            frozenset({1, 5}): (1, None),
            frozenset({9, 5}): (1, None),
            frozenset({11, 7}): (1, None),
            frozenset({10, 11}): (1, None),
            frozenset({11, 6}): (1, None),
        },
    )
    assert graph.sigma_radical_atom_keys(gra) == frozenset({0, 1, 3})


# # bond properties
def test__kekules_bond_orders_collated():
    """test graph.kekules_bond_orders_collated"""
    print(graph.kekules_bond_orders_collated(C3H3_CGR))
    assert all(
        set(v) == {1, 2}
        for _, v in graph.kekules_bond_orders_collated(C3H3_CGR).items()
    )


# # transformations
def test__kekules():
    """test graph.kekules"""
    # C=C[CH2]
    gra = (
        {0: ("C", 2, None), 1: ("C", 1, None), 2: ("C", 2, None)},
        {frozenset({0, 1}): (2, None), frozenset({1, 2}): (1, None)},
    )

    gras = graph.kekules(gra)
    print(len(gras))
    assert len(gras) == 2 and gra in gras

    # C=CC=C[CH2]
    gra = (
        {
            0: ("C", 2, None),
            1: ("C", 1, None),
            2: ("C", 1, None),
            3: ("C", 1, None),
            4: ("C", 2, None),
        },
        {
            frozenset({3, 4}): (1, None),
            frozenset({0, 1}): (2, None),
            frozenset({2, 3}): (2, None),
            frozenset({1, 2}): (1, None),
        },
    )

    gras = graph.kekules(gra)
    print(len(gras))
    assert len(gras) == 3 and gra in gras

    # C=C=C=C
    gra = (
        {0: ("C", 2, None), 1: ("C", 0, None), 2: ("C", 0, None), 3: ("C", 2, None)},
        {
            frozenset({0, 1}): (2, None),
            frozenset({2, 3}): (2, None),
            frozenset({1, 2}): (2, None),
        },
    )

    gras = graph.kekules(gra)
    print(len(gras))
    assert len(gras) == 1 and gra in gras

    # C=CC=CC=C
    gra = (
        {
            0: ("C", 2, None),
            1: ("C", 1, None),
            2: ("C", 1, None),
            3: ("C", 1, None),
            4: ("C", 1, None),
            5: ("C", 2, None),
        },
        {
            frozenset({3, 4}): (1, None),
            frozenset({2, 3}): (2, None),
            frozenset({1, 2}): (1, None),
            frozenset({4, 5}): (2, None),
            frozenset({0, 1}): (2, None),
        },
    )

    gras = graph.kekules(gra)
    print(len(gras))
    assert len(gras) == 1 and gra in gras

    # C1=CC=CC=C1 (benzene)
    gra = (
        {
            0: ("C", 1, None),
            1: ("C", 1, None),
            2: ("C", 1, None),
            3: ("C", 1, None),
            4: ("C", 1, None),
            5: ("C", 1, None),
        },
        {
            frozenset({3, 4}): (1, None),
            frozenset({2, 3}): (2, None),
            frozenset({1, 2}): (1, None),
            frozenset({4, 5}): (2, None),
            frozenset({0, 1}): (2, None),
            frozenset({0, 5}): (1, None),
        },
    )

    gras = graph.kekules(gra)
    print(len(gras))
    assert len(gras) == 2 and gra in gras

    # # C12=CC=C1C=C2  _  _
    # #              ||_||_||
    # gra = ({0: ('C', 0, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({0, 3}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({0, 5}): (1, None)})

    # # C1=CC=C2C=CC=CC2=C1 (naphthalene)
    # gra = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({6, 7}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({8, 3}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({0, 9}): (1, None)})

    # # C1=CC=C2C=C3C=CC=CC3=CC2=C1 (anthracene)
    # gra = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 0, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 1, None),
    #         9: ('C', 1, None), 10: ('C', 0, None), 11: ('C', 1, None),
    #         12: ('C', 0, None), 13: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({0, 13}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({11, 12}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
    #         frozenset({8, 7}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({3, 12}): (1, None), frozenset({10, 5}): (1, None),
    #         frozenset({12, 13}): (1, None), frozenset({10, 11}): (1, None)})

    # # C1=CC=C2C(=C1)C=CC3=CC=CC=C32 (phenanthrene)
    # gra = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None), 10: ('C', 1, None), 11: ('C', 1, None),
    #         12: ('C', 1, None), 13: ('C', 0, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({4, 6}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({11, 12}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
    #         frozenset({3, 13}): (1, None), frozenset({0, 5}): (1, None),
    #         frozenset({8, 7}): (1, None), frozenset({8, 13}): (1, None),
    #         frozenset({12, 13}): (1, None), frozenset({10, 11}): (1, None)})

    # # C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2 (pyrene)
    # gra = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None), 10: ('C', 1, None), 11: ('C', 1, None),
    #         12: ('C', 0, None), 13: ('C', 0, None), 14: ('C', 1, None),
    #         15: ('C', 1, None)},
    #        {frozenset({4, 6}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({11, 12}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 5}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({14, 15}): (1, None), frozenset({10, 11}): (1, None),
    #         frozenset({3, 4}): (1, None), frozenset({3, 13}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({6, 7}): (1, None),
    #         frozenset({8, 13}): (1, None), frozenset({12, 14}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({2, 15}): (1, None),
    #         frozenset({12, 13}): (1, None)})

    # # C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7 (coronene)
    # gra = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 0, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 0, None), 10: ('C', 0, None), 11: ('C', 0, None),
    #         12: ('C', 1, None), 13: ('C', 1, None), 14: ('C', 1, None),
    #         15: ('C', 1, None), 16: ('C', 0, None), 17: ('C', 0, None),
    #         18: ('C', 0, None), 19: ('C', 0, None), 20: ('C', 1, None),
    #         21: ('C', 1, None), 22: ('C', 1, None), 23: ('C', 1, None)},
    #        {frozenset({17, 10}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({11, 12}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 5}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({9, 4}): (1, None), frozenset({14, 15}): (1, None),
    #         frozenset({10, 11}): (1, None), frozenset({3, 4}): (1, None),
    #         frozenset({16, 17}): (1, None), frozenset({22, 23}): (1, None),
    #         frozenset({16, 15}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({18, 19}): (1, None), frozenset({18, 3}): (1, None),
    #         frozenset({11, 14}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({2, 21}): (1, None),
    #         frozenset({17, 18}): (1, None), frozenset({8, 13}): (1, None),
    #         frozenset({20, 21}): (1, None), frozenset({19, 20}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({16, 23}): (1, None),
    #         frozenset({19, 22}): (1, None), frozenset({12, 13}): (1, None)})


def test__kekule():
    """test graph.kekule"""
    assert graph.kekule(C3H3_CGR) in C3H3_RGRS


@pytest.mark.parametrize(
    "smi,res",
    [
        ("CCC", frozenset({})),
        ("C=CC", frozenset({0, 1})),
        ("C=C[CH2]", frozenset({0, 2})),
        ("C=CCC[CH2]", frozenset({4})),
        ("O=O", frozenset({0, 1})),
    ],
)
def test__addition_atom_keys(smi, res):
    """test graph.addition_atom_keys"""
    gra = automol.smiles.graph(smi)
    assert automol.graph.addition_atom_keys(gra) == res


@pytest.mark.parametrize(
    "smi,res",
    [
        ("COC", frozenset({})),
        ("[CH2]OC", frozenset({frozenset({3, 4})})),
    ],
)
def test__beta_scission_bond_keys(smi, res):
    """test graph.beta_scission_bond_keys"""
    gra = automol.smiles.graph(smi)
    assert automol.graph.beta_scission_bond_keys(gra) == res


def test__rotational_bond_keys():
    """test graph.rotational_bond_keys"""
    cgr = (
        {0: ("C", 2, None), 1: ("C", 2, None), 2: ("C", 1, None), 3: ("C", 1, None)},
        {
            frozenset({0, 2}): (1, None),
            frozenset({1, 3}): (1, None),
            frozenset({2, 3}): (1, None),
        },
    )
    cgr = graph.explicit(cgr)
    assert graph.rotational_bond_keys(cgr) == frozenset({frozenset({2, 3})})

    cgr = (
        {0: ("C", 3, None), 1: ("C", 3, None), 2: ("C", 2, None), 3: ("C", 2, None)},
        {
            frozenset({0, 2}): (1, None),
            frozenset({1, 3}): (1, None),
            frozenset({2, 3}): (1, None),
        },
    )
    assert graph.rotational_bond_keys(cgr) == frozenset(
        {frozenset({0, 2}), frozenset({1, 3}), frozenset({2, 3})}
    )
    assert graph.rotational_bond_keys(cgr, with_h_rotors=False) == frozenset(
        {frozenset({2, 3})}
    )
    assert graph.rotational_bond_keys(cgr, with_ch_rotors=False) == frozenset(
        {frozenset({2, 3})}
    )

    # Check that we don't misidentify rotational bond keys
    cgr = automol.smiles.graph("C#CC=C")
    assert graph.rotational_bond_keys(cgr) == frozenset()


def test__rotational_segment_keys():
    """test graph.rotational_segment_keys"""
    gra = (
        {
            0: ("C", 3, None),
            1: ("C", 0, None),
            2: ("C", 0, None),
            3: ("C", 0, None),
            4: ("C", 0, None),
            5: ("C", 0, None),
            6: ("C", 0, None),
            7: ("C", 2, None),
            8: ("C", 1, None),
            9: ("C", 3, None),
            10: ("C", 2, None),
            11: ("C", 3, None),
        },
        {
            frozenset({3, 4}): (1, None),
            frozenset({2, 3}): (1, None),
            frozenset({1, 2}): (1, None),
            frozenset({4, 5}): (1, None),
            frozenset({0, 1}): (1, None),
            frozenset({6, 7}): (1, None),
            frozenset({8, 9}): (1, None),
            frozenset({7, 8}): (1, None),
            frozenset({8, 10}): (1, None),
            frozenset({5, 6}): (1, None),
            frozenset({10, 11}): (1, None),
        },
    )
    assert graph.rotational_segment_keys(gra) == frozenset(
        {(0, 1, 2, 3, 4, 5, 6, 7), (7, 8), (8, 9), (8, 10), (10, 11)}
    )


def test__rotational_coordinates():
    """test graph.rotational_coordinates"""
    # CC#CC#CC#CCC(C)CC (z-matrix)
    geo = automol.smiles.geometry("CC#CC#CC#CCC(C)CC")
    zma = automol.geom.zmatrix(geo)
    zgra = automol.zmat.graph(zma, stereo=False, dummy=True)

    tors_dct = automol.zmat.torsion_coordinates(zma, zgra)
    coo_key_lst0 = set(map(tuple, map(reversed, tors_dct.values())))

    coo_key_lst = graph.rotational_coordinates(zgra, segment=False)
    assert coo_key_lst == coo_key_lst0, f"{coo_key_lst} != {coo_key_lst0}"

    # Make sure the segment version runs
    coo_key_lst = graph.rotational_coordinates(zgra)
    print(coo_key_lst)


def test__align_with_geometry():
    """test graph.align_with_geometry"""
    gra = automol.smiles.graph("C#CC1CCC1C#C")
    geo = graph.geometry(gra)
    zma, zc_ = automol.geom.zmatrix_with_conversion_info(geo, gra=gra)
    zgra = graph.apply_zmatrix_conversion(gra, zc_)
    geo_idx_dct = automol.util.zmat_conv.relabel_dict(zc_, "zmat")
    args = sorted(graph.atom_keys(zgra))

    # Check the alignment
    gra1, geo1, args1, key_dct, idx_dct = graph.align_with_geometry(
        zgra, geo, args, geo_idx_dct
    )
    assert graph.geometry_matches(gra1, geo1)
    assert args1.count(None) == 4
    assert [k for k in args1 if k is not None] == list(range(16))
    gra0 = automol.graph.relabel(gra1, key_dct)
    assert automol.graph.without_dummy_atoms(zgra) == gra0
    geo0 = automol.geom.reorder(geo1, idx_dct)
    assert automol.geom.almost_equal(geo, geo0)

    # Check that the alignment works with dummy atoms included
    zgeo = automol.zmat.geometry(zma, dummy=True)
    zgra1, zgeo1, zargs1, zkey_dct, zidx_dct = graph.align_with_geometry(
        zgra, zgeo, args
    )
    assert graph.geometry_matches(zgra1, zgeo1)
    assert zargs1 == list(range(20))
    zgra0 = automol.graph.relabel(zgra1, zkey_dct)
    assert zgra == zgra0
    zgeo0 = automol.geom.reorder(zgeo1, zidx_dct)
    assert automol.geom.almost_equal(zgeo, zgeo0)


def test__species__graph_conversion():
    """test interchanging between graphs aligned by zma and geo"""
    geo = automol.smiles.geometry("CC#CC#CCCCC#CC")
    gra = automol.geom.graph(geo)

    zma, dc_ = automol.geom.zmatrix_with_conversion_info(geo)
    zgra = graph.apply_zmatrix_conversion(gra, dc_)
    assert zgra != gra

    geo_ = automol.zmat.geometry(zma, zc_=dc_)
    gra_ = automol.geom.graph(geo_, stereo=False)
    assert gra == gra_ == graph.undo_zmatrix_conversion(zgra, dc_)

    # extra test case from Luna
    ich = "InChI=1S/C11H8/c1-3-10(4-2)11-8-6-5-7-9-11/h1,5-9H,2H2"
    geo = automol.chi.geometry(ich)
    gra = automol.geom.graph(geo)

    zma, dc_ = automol.geom.zmatrix_with_conversion_info(geo)
    zgra = graph.apply_zmatrix_conversion(gra, dc_)
    assert zgra != gra

    geo_ = automol.zmat.geometry(zma, zc_=dc_)
    gra_ = automol.geom.graph(geo_, stereo=False)
    assert gra == gra_ == graph.undo_zmatrix_conversion(zgra, dc_)


# stereo graph library
def test__geometry_atom_parity():
    """test graph.geometry_atom_parity"""
    # '[C@H](Cl)(F)(O)'
    geo = (
        ("C", (-0.156548, 0.194561, -0.651837)),
        ("H", (0.352681, 0.459932, -2.635957)),
        ("Cl", (-2.670412, -1.956931, -0.346205)),
        ("F", (-0.861626, 2.505747, 0.34212)),
        ("O", (1.916674, -0.693354, 0.754657)),
        ("H", (1.419231, -0.509954, 2.537223)),
    )
    gra = automol.geom.graph(geo)
    assert graph.geometry_atom_parity(gra, geo, 0) is False

    # '[C@@H](Cl)(F)(O)'
    geo = (
        ("C", (-0.317722, 0.059267, -0.618703)),
        ("H", (-0.26878, 0.123481, -2.68267)),
        ("Cl", (-1.517891, 2.894393, 0.631708)),
        ("F", (-1.842322, -1.93151, 0.115519)),
        ("O", (2.108174, -0.367259, 0.380286)),
        ("H", (1.838541, -0.778373, 2.17386)),
    )
    gra = automol.geom.graph(geo)
    assert graph.geometry_atom_parity(gra, geo, 0) is True


def test__geometry_bond_parity():
    """test graph.geometry_bond_parity"""
    # r'F/C=N/[H]'
    geo = (
        ("F", (-2.483972, -1.744981, 0.024901)),
        ("C", (-0.860131, 0.558467, -0.008463)),
        ("N", (1.9591, 0.223915, -0.002683)),
        ("H", (3.348622, 1.709226, -0.024014)),
        ("H", (-1.577309, 2.507937, -0.03731)),
    )
    gra = automol.geom.graph(geo)
    assert graph.geometry_bond_parity(gra, geo, [1, 2]) is True

    # r'F/C=N\[H]'
    geo = (
        ("F", (-1.495767, 2.657838, -0.288402)),
        ("C", (-0.969504, -0.04963, -0.822887)),
        ("N", (1.620729, -1.139643, -0.788144)),
        ("H", (3.359017, -0.121361, -0.3958)),
        ("H", (-2.514475, -1.347203, -1.263959)),
    )
    gra = automol.geom.graph(geo)
    assert graph.geometry_bond_parity(gra, geo, [1, 2]) is False


def test__geometries_parity_mismatches():
    """test graph.geometries_parity_mismatches"""
    # F/C=N/[C@H](O)(F)
    geo1 = (
        ("F", (5.084539, -0.513665, -0.452995)),
        ("C", (2.63457, -0.807065, -0.671966)),
        ("N", (1.185305, 0.951476, 0.193626)),
        ("C", (-1.474258, 0.369908, -0.151588)),
        ("H", (-2.011785, 0.417195, -2.149368)),
        ("O", (-3.095413, 1.983267, 1.196529)),
        ("F", (-1.982076, -1.995096, 0.779896)),
        ("H", (1.893063, -2.574706, -1.552191)),
        ("H", (-2.233945, 2.168686, 2.808057)),
    )
    # F/C=N/[C@@H](O)(F)
    geo2 = (
        ("F", (3.769058, -2.906815, -2.104708)),
        ("C", (2.30642, -1.089166, -1.133306)),
        ("N", (0.446289, -1.783733, 0.261487)),
        ("C", (-1.006658, 0.368806, 1.181072)),
        ("H", (0.105919, 1.401487, 2.627576)),
        ("O", (-1.840221, 2.050204, -0.652784)),
        ("F", (-3.015519, -0.619458, 2.46581)),
        ("H", (2.796282, 0.888302, -1.585397)),
        ("H", (-3.561571, 1.690373, -1.05975)),
    )
    # F/C=N\[C@H](O)(F)
    geo3 = (
        ("F", (2.548759, -2.658852, 0.888608)),
        ("C", (2.755227, -0.684142, -0.677915)),
        ("N", (0.850395, 0.749212, -1.139181)),
        ("C", (-1.378559, -0.016047, 0.245876)),
        ("H", (-1.812511, -2.029738, 0.040561)),
        ("O", (-3.548738, 1.369169, -0.422483)),
        ("F", (-0.993379, 0.509014, 2.724287)),
        ("H", (4.599124, -0.330173, -1.574926)),
        ("H", (-3.020319, 3.091559, -0.084827)),
    )
    # F/C=N\[C@@H](O)(F)
    geo4 = (
        ("F", (0.178094, -1.31074, -2.345725)),
        ("C", (-1.723374, 0.199865, -1.629498)),
        ("N", (-1.704667, 1.254529, 0.554373)),
        ("C", (0.533485, 0.62288, 2.036335)),
        ("H", (0.438289, 1.705675, 3.82386)),
        ("O", (0.864003, -1.917842, 2.632388)),
        ("F", (2.594726, 1.541084, 0.791643)),
        ("H", (-3.28618, 0.525251, -2.961819)),
        ("H", (2.105626, -2.620703, 1.525167)),
    )
    gra = automol.geom.graph(geo1)
    keys = [3, (1, 2)]

    assert not graph.geometries_parity_mismatches(gra, geo1, geo1, keys)
    assert graph.geometries_parity_mismatches(gra, geo1, geo2, keys) == (3,)
    assert graph.geometries_parity_mismatches(gra, geo1, geo3, keys) == ((1, 2),)
    assert graph.geometries_parity_mismatches(gra, geo1, geo4, keys) == (3, (1, 2))


def test__stereogenic_keys():
    """test graph.stereogenic_atom_keys"""
    assert graph.unassigned_stereocenter_keys(C8H13O_CGR) == frozenset(
        {
            6,
            7,
            frozenset({3, 5}),
        }
    )
    assert graph.unassigned_stereocenter_keys(C3H3CL2F3_CGR) == frozenset({1, 2})
    assert graph.unassigned_stereocenter_keys(C3H5N3_CGR) == frozenset(
        {frozenset({1, 4}), frozenset({0, 3})}
    )

    # Atoms only
    assert graph.unassigned_stereocenter_keys(C8H13O_CGR, bond=False) == frozenset(
        {6, 7}
    )

    # Bonds only
    assert graph.unassigned_stereocenter_keys(C8H13O_CGR, atom=False) == frozenset(
        {
            frozenset({3, 5}),
        }
    )

    cgr = (
        {0: ("C", 2, None), 1: ("C", 3, None), 2: ("C", 1, None), 3: ("O", 1, None)},
        {
            frozenset({0, 2}): (1, None),
            frozenset({2, 3}): (1, None),
            frozenset({1, 2}): (1, None),
        },
    )
    print(graph.unassigned_stereocenter_keys(cgr))
    assert graph.unassigned_stereocenter_keys(cgr) == frozenset({2})

    # Bug fix:
    cgr = (
        {
            0: ("C", 3, None),
            1: ("C", 3, None),
            2: ("C", 2, None),
            3: ("C", 2, None),
            4: ("C", 2, None),
            5: ("C", 2, None),
            6: ("C", 2, None),
            7: ("C", 2, None),
            8: ("C", 2, None),
            9: ("C", 2, None),
            10: ("C", 1, None),
            11: ("C", 1, None),
            12: ("O", 0, None),
        },
        {
            frozenset({4, 6}): (1, None),
            frozenset({11, 12}): (1, None),
            frozenset({8, 9}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({10, 6}): (1, None),
            frozenset({8, 10}): (1, None),
            frozenset({2, 4}): (1, None),
            frozenset({9, 11}): (1, None),
            frozenset({10, 12}): (1, None),
            frozenset({3, 5}): (1, None),
            frozenset({11, 7}): (1, None),
            frozenset({1, 3}): (1, None),
            frozenset({5, 7}): (1, None),
        },
    )
    print(graph.unassigned_stereocenter_keys(cgr))
    assert graph.unassigned_stereocenter_keys(cgr) == frozenset({10, 11})


def test__expand_stereo():
    """test graph.expand_stereo"""
    print(graph.smiles(C2H2CL2F2_CGR))
    assert graph.expand_stereo(C2H2CL2F2_CGR) == C2H2CL2F2_SGRS
    assert graph.expand_stereo(C3H3CL2F3_CGR) == C3H3CL2F3_SGRS
    # When symmetry equivalents are filtered out, we can't guarantee that the
    # sequence will match exactly, but they will be isomorphic sequences.
    assert graph.sequence_isomorphism(
        graph.expand_stereo(C3H5N3_CGR), C3H5N3_SGRS, stereo=True
    )
    assert graph.expand_stereo(C8H13O_CGR) == C8H13O_SGRS

    # CC(OO)C(O[O])C(OO)C
    gra = (
        {
            0: ("C", 3, None),
            1: ("C", 3, None),
            2: ("C", 1, None),
            3: ("C", 1, None),
            4: ("C", 1, None),
            5: ("O", 1, None),
            6: ("O", 1, None),
            7: ("O", 0, None),
            8: ("O", 0, None),
            9: ("O", 0, None),
            10: ("O", 0, None),
        },
        {
            frozenset({10, 4}): (1, None),
            frozenset({8, 2}): (1, None),
            frozenset({3, 4}): (1, None),
            frozenset({9, 6}): (1, None),
            frozenset({9, 3}): (1, None),
            frozenset({10, 7}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({2, 4}): (1, None),
            frozenset({8, 5}): (1, None),
            frozenset({1, 3}): (1, None),
        },
    )
    assert len(graph.expand_stereo(gra, enant=True, symeq=True)) == 6
    assert len(graph.expand_stereo(gra, enant=True, symeq=False)) == 4
    assert len(graph.expand_stereo(gra, enant=False, symeq=True)) == 5
    assert len(graph.expand_stereo(gra, enant=False, symeq=False)) == 3

    # 'FC=CF.[CH2]C(O)C'
    gra = (
        {
            0: ("C", 0, None),
            1: ("C", 0, None),
            2: ("C", 0, None),
            3: ("O", 0, None),
            4: ("H", 0, None),
            5: ("H", 0, None),
            6: ("H", 0, None),
            7: ("H", 0, None),
            8: ("H", 0, None),
            9: ("H", 0, None),
            10: ("H", 0, None),
            11: ("C", 0, None),
            12: ("C", 0, None),
            13: ("F", 0, None),
            14: ("F", 0, None),
            15: ("H", 0, None),
            16: ("H", 0, None),
        },
        {
            frozenset({10, 3}): (1, None),
            frozenset({11, 12}): (1, None),
            frozenset({1, 2}): (1, None),
            frozenset({8, 1}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({16, 12}): (1, None),
            frozenset({0, 5}): (1, None),
            frozenset({1, 6}): (1, None),
            frozenset({1, 7}): (1, None),
            frozenset({2, 3}): (1, None),
            frozenset({0, 4}): (1, None),
            frozenset({11, 13}): (1, None),
            frozenset({11, 15}): (1, None),
            frozenset({9, 2}): (1, None),
            frozenset({12, 14}): (1, None),
        },
    )
    assert len(graph.expand_stereo(gra, enant=True, symeq=True)) == 4
    assert len(graph.expand_stereo(gra, enant=False, symeq=True)) == 2


@pytest.mark.parametrize(
    "smi,npars1,npars2,par_dct",
    [
        ("C1C2OC2CC1", 3, 1, {1: False, 3: True}),
        ("C1(O2)CC2CC1", 3, 1, {0: False, 3: True}),
        ("C=1C2CC(C=1)O2", 3, 1, {1: False, 3: True}),
    ],
)
def test__expand_stereo__strained(smi, npars1, npars2, par_dct):
    """test removal of strained stereoisomers from stereoexpansion"""
    print(f"smi = {smi}")
    gra = automol.smiles.graph(smi)

    # Without removing strained stereoisomers, there are three distinct possibilities
    sgras = graph.expand_stereo(gra, strained=True)
    assert len(sgras) == npars1, f"{len(sgras)} != {npars1}"

    # With removing strained stereoisomers, there is only one distinct possibility
    sgras = graph.expand_stereo(gra, strained=False)
    assert len(sgras) == npars2, f"{len(sgras)} != {npars2}"

    (sgra,) = sgras
    par_dct_ = automol.util.dict_.filter_by_value(
        automol.graph.stereo_parities(sgra), lambda x: x is not None
    )
    assert par_dct_ == par_dct, f"{par_dct_} != {par_dct}"


def test__ring_systems():
    """test graph.ring_systems"""
    chi = automol.smiles.chi("C12CC(C1)C2CC3C(C3)CCC4C5CCC(CC5)C4")
    gra = automol.chi.graph(chi)
    rsys = sorted(graph.ring_systems(gra), key=graph.atom_count)
    assert list(map(graph.atom_count, rsys)) == [7, 12, 21]


def test__equivalent_atoms():
    """test graph.equivalent_atoms"""
    # central carbon
    assert graph.equivalent_atoms(C4H10_GRA, 3) == {3}
    # central hydrogen
    assert graph.equivalent_atoms(C4H10_GRA, 13) == {13}
    # terminal carbons
    assert graph.equivalent_atoms(C4H10_GRA, 0) == {0, 1, 2}
    assert graph.equivalent_atoms(C4H10_GRA, 1) == {0, 1, 2}
    assert graph.equivalent_atoms(C4H10_GRA, 2) == {0, 1, 2}
    # terminal hydrogens
    assert graph.equivalent_atoms(C4H10_GRA, 4) == {4, 5, 6, 7, 8, 9, 10, 11, 12}
    assert graph.equivalent_atoms(C4H10_GRA, 5) == {4, 5, 6, 7, 8, 9, 10, 11, 12}
    assert graph.equivalent_atoms(C4H10_GRA, 6) == {4, 5, 6, 7, 8, 9, 10, 11, 12}
    assert graph.equivalent_atoms(C4H10_GRA, 11) == {4, 5, 6, 7, 8, 9, 10, 11, 12}
    assert graph.equivalent_atoms(C4H10_GRA, 12) == {4, 5, 6, 7, 8, 9, 10, 11, 12}


def test__equivalent_bonds():
    """test graph.equivalent_atoms"""
    assert graph.equivalent_bonds(C4H10_GRA, (2, 3)) == {(0, 3), (1, 3), (2, 3)}


def test__vmat__vmatrix():
    """test graph.vmat.vmatrix"""
    chi = automol.smiles.chi("C12CC(C1)C2CC3C(C3)CCC4C5CCC(CC5)C4")
    gra = automol.chi.graph(chi)
    _, zma_keys = graph.vmat.vmatrix(gra)
    assert set(zma_keys) == graph.atom_keys(gra)


def test__canonical():
    """test graph.canonical"""

    def _test_from_smiles(smi):
        chi = automol.smiles.chi(smi)
        print(chi)
        geo = automol.chi.geometry(chi)
        print(automol.geom.string(geo))
        gra = automol.geom.graph(geo)

        gra = graph.implicit(gra)

        can_gra = graph.canonical(gra)
        print(graph.string(can_gra, one_indexed=True))

        natms = len(graph.atom_keys(gra))

        for _ in range(10):
            pmt = list(map(int, numpy.random.permutation(natms)))
            print(pmt)
            pmt_gra = graph.relabel(gra, dict(enumerate(pmt)))
            can_pmt_gra = graph.canonical(pmt_gra)
            print(graph.string(can_pmt_gra, one_indexed=True))
            assert can_pmt_gra == can_gra

    # More tests can be added here
    smis = [
        "c1ccccc1CC",
        "F[C@@H]([C@@H](F)Cl)[C@H](F)Cl",
        r"[H]/N=C\C(\C=N\[H])=N\[H]",
    ]
    for smi in smis:
        _test_from_smiles(smi)


def test__calculate_priorities_and_assign_parities():
    """test graph.calculate_priorities_and_assign_parities"""

    def _test_from_smiles(smi, ref_atm_pars, ref_bnd_pars):
        print(smi)
        chi = automol.smiles.chi(smi)
        print(chi)
        geo = automol.chi.geometry(chi)
        gra = automol.geom.graph(geo)

        par_eval_ = graph.parity_evaluator_measure_from_geometry_(geo)

        gra, _, pri_dct, *_ = graph.calculate_stereo(gra, par_eval_=par_eval_)

        print(pri_dct)

        atm_par_dct = automol.util.dict_.filter_by_value(
            graph.atom_stereo_parities(gra), lambda x: x is not None
        )
        bnd_par_dct = automol.util.dict_.filter_by_value(
            graph.bond_stereo_parities(gra), lambda x: x is not None
        )
        # Frozensets don't sort properly, so use sorted tuples for the keys
        bnd_par_dct = automol.util.dict_.transform_keys(
            bnd_par_dct, lambda x: tuple(sorted(x))
        )

        atm_pars = [p for k, p in sorted(atm_par_dct.items()) if p is not None]
        bnd_pars = [p for k, p in sorted(bnd_par_dct.items()) if p is not None]

        print(ref_atm_pars, ref_bnd_pars)
        print(atm_pars, bnd_pars)
        assert atm_pars == ref_atm_pars
        assert bnd_pars == ref_bnd_pars

    # More tests can be added here
    args = [
        # Atom stereo tests
        ("C[C@@](F)(N)O", [True], []),  # R = clock  = '+' => True
        ("C[C@](F)(N)O", [False], []),  # S = aclock = '-' => False
        ("[C@@H](F)(N)O", [True], []),  # R = clock  = '+' => True
        ("[C@H](F)(N)O", [False], []),  # S = aclock = '-' => False
        # Bond stereo tests
        (r"F/C=C/F", [], [True]),  # trans = '+' => True
        (r"F/C=C\F", [], [False]),  # cis   = '-' => False
        (r"[H]/N=C(Cl)/F", [], [True]),  # trans = '+' => True
        (r"[H]/N=C(Cl)\F", [], [False]),  # cis   = '-' => False
        (r"[H]/N=C/F", [], [True]),  # trans = '+' => True
        (r"[H]/N=C\F", [], [False]),  # cis   = '-' => False
        (r"[H]/N=N/[H]", [], [True]),  # trans = '+' => True
        (r"[H]/N=N\[H]", [], [False]),  # cis   = '-' => False
        # Advanced tests
        ("F[C@@H]([C@@H](F)Cl)[C@H](F)Cl", [True, True, False], []),
        (r"[H]/N=C\C(\C=N\[H])=N\[H]", [], [False, True, False]),
    ]
    for smi, ref_atm_pars, ref_bnd_pars in args:
        _test_from_smiles(smi, ref_atm_pars, ref_bnd_pars)


def test__to_local_stereo():
    """test graph.to_local_stereo"""
    # Atom parity test:
    # Indices 3 and 4 are swapped so that local stereo will have opposite
    # parity
    can_gra = (
        {
            0: ("C", 3, None),
            1: ("C", 0, True),
            2: ("F", 0, None),
            4: ("N", 2, None),
            3: ("O", 1, None),
        },
        {
            frozenset({0, 1}): (1, None),
            frozenset({1, 4}): (1, None),
            frozenset({1, 3}): (1, None),
            frozenset({1, 2}): (1, None),
        },
    )
    can_par = graph.atom_stereo_parities(can_gra)[1]

    loc_gra = graph.to_local_stereo(can_gra)
    loc_par = graph.atom_stereo_parities(loc_gra)[1]
    print(can_par, loc_par)
    assert can_par is True
    assert loc_par is False

    # Bond parity test:
    # Indices 3 and 5 are swapped so that local stereo will have opposite
    # parity
    can_gra = (
        {
            0: ("C", 0, None),
            1: ("C", 0, None),
            2: ("Cl", 0, None),
            5: ("Cl", 0, None),
            4: ("F", 0, None),
            3: ("F", 0, None),
        },
        {
            frozenset({0, 1}): (1, True),
            frozenset({0, 2}): (1, None),
            frozenset({0, 4}): (1, None),
            frozenset({1, 3}): (1, None),
            frozenset({1, 5}): (1, None),
        },
    )
    can_par = graph.bond_stereo_parities(can_gra)[frozenset({0, 1})]

    loc_gra = graph.to_local_stereo(can_gra)
    loc_par = graph.bond_stereo_parities(loc_gra)[frozenset({0, 1})]
    print(can_par, loc_par)
    assert can_par is True
    assert loc_par is False


def test__from_local_stereo():
    """test graph.from_local_stereo"""
    gra = graph.explicit(C3H5_SGR)
    loc_gra = graph.to_local_stereo(gra)
    assert gra == graph.from_local_stereo(loc_gra)

    for gra in C2H2CL2F2_SGRS:
        gra = graph.explicit(gra)
        loc_gra = graph.to_local_stereo(gra)
        assert gra == graph.from_local_stereo(loc_gra)

    for gra in C3H3CL2F3_SGRS:
        gra = graph.explicit(gra)
        loc_gra = graph.to_local_stereo(gra)
        assert gra == graph.from_local_stereo(loc_gra)

    for gra in C3H5N3_SGRS:
        gra = graph.explicit(gra)
        loc_gra = graph.to_local_stereo(gra)
        assert gra == graph.from_local_stereo(loc_gra)

    for gra in C8H13O_SGRS:
        gra = graph.explicit(gra)
        loc_gra = graph.to_local_stereo(gra)
        assert gra == graph.from_local_stereo(loc_gra)


def test__has_resonance_bond_stereo():
    """test graph.has_resonance_bond_stereo"""
    gra = (
        {
            0: ("F", 0, None),
            1: ("C", 1, None),
            3: ("C", 1, None),
            4: ("C", 1, None),
            5: ("F", 0, None),
        },
        {
            frozenset({3, 4}): (1, True),
            frozenset({0, 1}): (1, None),
            frozenset({1, 3}): (1, True),
            frozenset({4, 5}): (1, None),
        },
    )
    assert graph.has_resonance_bond_stereo(gra)

    gra = (
        {
            0: ("F", 0, None),
            1: ("C", 2, None),
            4: ("C", 1, None),
            5: ("C", 1, None),
            6: ("F", 0, None),
        },
        {
            frozenset({4, 5}): (1, True),
            frozenset({0, 1}): (1, None),
            frozenset({1, 4}): (1, None),
            frozenset({5, 6}): (1, None),
        },
    )
    assert not graph.has_resonance_bond_stereo(gra)


def test__inchi_is_bad():
    """test graph.inchi_is_bad"""
    # This species is missing resonance bond stereo
    gra = (
        {
            0: ("F", 0, None),
            1: ("C", 1, None),
            3: ("C", 1, None),
            4: ("C", 1, None),
            5: ("F", 0, None),
        },
        {
            frozenset({3, 4}): (1, True),
            frozenset({0, 1}): (1, None),
            frozenset({1, 3}): (1, True),
            frozenset({4, 5}): (1, None),
        },
    )
    ich = graph.inchi(gra)
    print(ich)
    assert graph.inchi_is_bad(gra, ich)

    # This species is missing vinyl radical bond stereo
    gra = (
        {
            0: ("C", 1, None),
            1: ("C", 1, None),
            2: ("C", 2, None),
            3: ("C", 1, None),
            4: ("C", 0, None),
            5: ("C", 1, None),
            6: ("C", 1, None),
            7: ("C", 0, None),
            8: ("C", 1, None),
            9: ("C", 1, None),
            10: ("C", 0, None),
        },
        {
            frozenset({9, 6}): (1, None),
            frozenset({9, 10}): (1, None),
            frozenset({10, 7}): (1, None),
            frozenset({1, 2}): (1, None),
            frozenset({0, 1}): (1, True),
            frozenset({3, 6}): (1, None),
            frozenset({8, 10}): (1, None),
            frozenset({2, 4}): (1, None),
            frozenset({3, 5}): (1, None),
            frozenset({8, 5}): (1, None),
            frozenset({4, 7}): (1, None),
        },
    )
    ich = graph.inchi(gra)
    print(ich)
    assert graph.inchi_is_bad(gra, ich)

    # This species has mobile hydrogens
    gra = (
        {0: ("C", 0, None), 1: ("O", 1, None), 2: ("O", 0, None)},
        {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None)},
    )
    ich = graph.inchi(gra)
    print(ich)
    assert graph.inchi_is_bad(gra, ich)


def test__amchi():
    """test graph.amchi"""
    # bond stereo
    gra = (
        {
            0: ("C", 1, None),
            1: ("C", 1, None),
            2: ("C", 0, None),
            3: ("N", 1, None),
            4: ("N", 1, None),
            5: ("N", 1, None),
        },
        {
            frozenset({1, 4}): (1, True),
            frozenset({1, 2}): (1, None),
            frozenset({0, 3}): (1, False),
            frozenset({0, 2}): (1, None),
            frozenset({2, 5}): (1, False),
        },
    )
    chi = graph.amchi(gra)
    print(chi)

    assert chi == "AMChI=1/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3-"

    # atom stereo
    gra = (
        {
            0: ("C", 1, None),
            1: ("C", 1, True),
            2: ("C", 1, True),
            3: ("Cl", 0, None),
            4: ("Cl", 0, None),
            5: ("F", 0, None),
            6: ("F", 0, None),
            7: ("F", 0, None),
        },
        {
            frozenset({0, 1}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({0, 5}): (1, None),
            frozenset({2, 4}): (1, None),
            frozenset({1, 3}): (1, None),
            frozenset({1, 6}): (1, None),
            frozenset({2, 7}): (1, None),
        },
    )
    chi = graph.amchi(gra)
    print(chi)

    assert chi == "AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t2-,3-/m1/s1"


def test__amchi_with_numbers():
    """test graph.amchi"""
    # bond stereo
    gra1 = (
        {
            0: ("C", 1, None),
            1: ("C", 1, None),
            2: ("C", 0, None),
            3: ("N", 1, None),
            4: ("N", 1, None),
            5: ("N", 1, None),
        },
        {
            frozenset({1, 4}): (1, True),
            frozenset({1, 2}): (1, None),
            frozenset({0, 3}): (1, False),
            frozenset({0, 2}): (1, None),
            frozenset({2, 5}): (1, False),
        },
    )
    chi1, num_dcts1 = graph.amchi_with_numbers(gra1)
    print(chi1)
    print(num_dcts1)

    assert chi1 == "AMChI=1/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3-"

    # atom stereo
    gra2 = (
        {
            0: ("C", 1, None),
            1: ("C", 1, True),
            2: ("C", 1, True),
            3: ("Cl", 0, None),
            4: ("Cl", 0, None),
            5: ("F", 0, None),
            6: ("F", 0, None),
            7: ("F", 0, None),
        },
        {
            frozenset({0, 1}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({0, 5}): (1, None),
            frozenset({2, 4}): (1, None),
            frozenset({1, 3}): (1, None),
            frozenset({1, 6}): (1, None),
            frozenset({2, 7}): (1, None),
        },
    )
    chi2, num_dcts2 = graph.amchi_with_numbers(gra2)
    print(chi2)
    print(num_dcts2)

    assert chi2 == "AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t2-,3-/m1/s1"

    gra = graph.union_from_sequence([gra1, gra2], shift_keys=True)
    chi, num_dcts = graph.amchi_with_numbers(gra)
    print(chi)
    print(num_dcts)

    assert chi == (
        "AMChI=1/C3H3Cl2F3.C3H5N3/c4-2(7)1(6)3(5)8;4-1-3(6)2-5/"
        "h1-3H;1-2,4-6H/b;4-1-,5-2+,6-3-/t2-,3-;/m1./s1"
    )
    assert num_dcts == (
        {6: 0, 7: 1, 8: 2, 9: 3, 10: 4, 11: 5, 12: 6, 13: 7},
        {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5},
    )


def test__smiles():
    """test graph.smiles"""
    smis = [
        # Rings:
        r"[CH2]/C=C1\[N@H]C(=C1)O[O]",
        "C1CC1",
        # Stereo atoms:
        "N[C@](C)(F)C(=O)O",
        "N[C@H](C)C(=O)O",
        "C[C@H]1CCCCO1",
        # Stereo bonds:
        r"F/C=C/F",
        r"F/C=C\F",
        r"CC/C=C/C=C/CF",
        r"C1CCCCCCCCCC/N=N/1",
        r"[H]/N=N/[H]",
        r"[H]/N=N/N=N\[H]",
        "c1ccc2c(c1)Cc1c2cccc1",
    ]

    for smi in smis:
        print()
        print("STARTING SMILES:", smi)
        chi = automol.smiles.chi(smi)
        geo = automol.chi.geometry(chi)
        gra = automol.geom.graph(geo)
        print(graph.string(gra))
        print(automol.geom.string(geo))
        smi = graph.smiles(gra)
        print("smiles from code:", smi)

        chi_smi = automol.chi.smiles(chi)
        print("chi:", chi)
        print("smiles from chi:", chi_smi)

        schi = automol.smiles.chi(smi)
        print("chi from smiles:", schi)
        # print(graph.string(GRA, one_indexed=True))
        print(graph.rings_atom_keys(gra))
        assert schi == chi


def test__smiles__with_resonance():
    """test graph.smiles"""

    smis = [
        r"FC=C-C=C-[CH]O",
        r"F[CH]C=CF",
    ]
    for smi in smis:
        print()
        print("STARTING SMILES:", smi)
        chi = automol.smiles.chi(smi)
        geo = automol.chi.geometry(chi)
        gra = automol.geom.graph(geo)
        print(graph.string(gra))
        print(automol.geom.string(geo))

        ste_keys = graph.bond_stereo_keys(gra)
        nste = len(ste_keys)
        ste_pars_lst = list(itertools.product(map(bool, range(2)), repeat=nste))
        for ste_pars in ste_pars_lst:
            bnd_par_dct = dict(zip(ste_keys, ste_pars))
            print(bnd_par_dct)
            gra = graph.set_bond_stereo_parities(gra, bnd_par_dct)
            ach = graph.amchi(gra)
            print("amchi:", ach)
            smi = graph.smiles(gra)
            print("smiles from code:", smi)
            print()


def test__perturb_geometry_planar_dihedrals():
    """test graph.perturb_geometry_planar_dihedrals()"""
    geo = (
        ("F", (3.45091938160, -1.42093905806, -0.36575882605)),
        ("C", (1.32805118466, 0.33015256925, -0.48245655954)),
        ("C", (-1.32797704904, -0.32991985288, -0.37225479287)),
        ("F", (-3.44911016186, 1.42664994742, -0.41428201235)),
        ("H", (2.12703663434, 2.19198115968, -0.67990730009)),
        ("H", (-2.12891998970, -2.19792476541, -0.25903256241)),
    )
    gra = automol.geom.graph(geo)
    dih = automol.geom.dihedral_angle(geo, 0, 1, 2, 3, degree=True)
    assert not numpy.isclose(dih, 170.0)

    geo = graph.perturb_geometry_planar_dihedrals(gra, geo, ang=10.0, degree=True)

    dih = automol.geom.dihedral_angle(geo, 0, 1, 2, 3, degree=True)
    assert numpy.isclose(dih, 170.0)


def test__stereo_corrected_geometry():
    """test graph.stereo_corrected_geometry"""
    geo = automol.smiles.geometry("C1C(F)C12C(F)C2")
    gra = automol.geom.graph(geo)
    sgras = graph.expand_stereo(gra)
    for sgra in sgras:
        print(graph.smiles(sgra))
        sgeo = graph.stereo_corrected_geometry(sgra, geo)
        sgra_from_geo = automol.geom.graph(sgeo)
        assert sgra_from_geo == sgra


def test__embed__clean_geometry():
    """test graph.embed.clean_geometry"""
    # Make sure good geometries don't get messed up by cleaning
    # F/C=N\C#C
    geo1 = (
        ("C", (-3.321773, 0.223278, -1.798954)),
        ("C", (-1.526479, 0.653524, -0.482726)),
        ("N", (0.453152, 1.154608, 1.00976)),
        ("C", (2.562887, 0.064839, 0.451482)),
        ("F", (2.694212, -1.636094, -1.420516)),
        ("H", (-4.917597, -0.155311, -2.962995)),
        ("H", (4.300424, 0.512127, 1.504603)),
    )
    gra = automol.geom.graph(geo1)
    geo2 = graph.clean_geometry(gra, geo1)
    assert automol.geom.almost_equal(geo1, geo2)


def test__atom_hypervalencies():
    """also tests geometry conversion for this case"""
    geo = (
        ("C", (0.1604, -0.2027, 0.0202)),
        ("H", (-1.3399, 1.6931, -0.169)),
        ("H", (-0.7776, -1.5753, -1.1172)),
        ("H", (0.1415, -0.5361, 2.0067)),
        ("H", (1.8671, 0.556, -0.7342)),
        ("O", (-2.8677, 3.6235, -0.3617)),
        ("H", (-4.0545, 2.5878, -1.3276)),
    )
    gra1 = automol.geom.graph(geo, fix_hyper=False)
    nhyp_dct1 = automol.graph.atom_hypervalencies(gra1)
    assert automol.util.dict_.by_value(nhyp_dct1) == {1: 1}

    gra2 = automol.geom.graph(geo, fix_hyper=True)
    nhyp_dct2 = automol.graph.atom_hypervalencies(gra2)
    assert automol.util.dict_.by_value(nhyp_dct2) == {}


@pytest.mark.parametrize(
    "gra, lin_keys, extend, lin_cap_dct0",
    [
        (C2H7_CGR, (5,), False, {(5,): (8, 5)}),
    ],
)
def test__linear_segments_cap_keys(gra, lin_keys, extend, lin_cap_dct0):
    """Test graph.linear_segments_cap_keys."""
    lin_cap_dct = graph.linear_segment_cap_keys(
        gra=gra, lin_keys=lin_keys, extend=extend
    )
    assert lin_cap_dct == lin_cap_dct0, f"{lin_cap_dct} != {lin_cap_dct0}"


if __name__ == "__main__":
    # test__to_local_stereo()

    # test__has_resonance_bond_stereo()
    # test__amchi_with_indices()
    # test__stereogenic_keys()
    # test__kekules_bond_orders_collated()
    # test__inchi_is_bad()
    # test__expand_stereo()
    # test__species__graph_conversion()
    # test__canonical()
    # test__calculate_priorities_and_assign_parities()
    # test__smiles()
    # test__kekules()
    # test__geometry_atom_parity()
    # test__geometry_bond_parity()
    # test__geometries_parity_mismatches()
    # test__unique()
    # test__isomorphic()
    # test__branch()
    # test__perturb_geometry_planar_dihedrals()
    # test__from_data()
    # test__setters()
    # test__atom_count()
    # test__atom_hybridizations()
    # test__kekules()
    # test__stereo_corrected_geometry()
    # test__embed__clean_geometry()
    # test__rotational_coordinates()
    # test__stereo_corrected_geometry()
    # test__atom_hypervalencies()
    # test__align_with_geometry()
    # test__embed__clean_geometry()
    # test__expand_stereo__strained("C=1C2CC(C=1)O2", 3, 1, {1: False, 3: True})
    # test__expand_stereo()
    test__linear_segments_cap_keys(C2H7_CGR, (5,), False, {(5,): (8, 5)})
