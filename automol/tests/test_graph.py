""" test automol.graph
"""

import itertools
import numpy
import automol
from automol import graph


# Vinyl radical with E/Z stereo
C3H5_SGR = (
    {0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
     3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
     6: ('H', 0, None), 7: ('H', 0, None)},
    {frozenset({0, 2}): (1, False), frozenset({0, 3}): (1, None),
     frozenset({1, 2}): (1, None), frozenset({1, 4}): (1, None),
     frozenset({1, 5}): (1, None), frozenset({1, 6}): (1, None),
     frozenset({2, 7}): (1, None)})


C8H13O_CGR = (
    {0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 1, None), 7: ('C', 1, None), 8: ('O', 0, None)},
    {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
     frozenset({3, 5}): (1, None), frozenset({5, 7}): (1, None)})
C8H13O_RGR = (
    {0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 1, None), 7: ('C', 1, None), 8: ('O', 0, None)},
    {frozenset({1, 4}): (2, None), frozenset({4, 6}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
     frozenset({3, 5}): (2, None), frozenset({5, 7}): (1, None)})
C8H13O_SGR = (
    {0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 1, False), 7: ('C', 1, False), 8: ('O', 0, None)},
    {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
     frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)})


C3H3_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
    {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
     frozenset({2, 0}): (1, None)})
C3H3_RGRS = (
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
     {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
      frozenset({2, 0}): (1, None)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
     {frozenset({0, 1}): (1, None), frozenset({1, 2}): (2, None),
      frozenset({2, 0}): (1, None)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
     {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
      frozenset({2, 0}): (2, None)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
     {frozenset({0, 1}): (2, None), frozenset({1, 2}): (1, None),
      frozenset({2, 0}): (1, None)}),
)

C2_CGR = ({0: ('C', 0, None), 1: ('C', 0, None)},
          {frozenset({0, 1}): (1, None)})
C2_RGRS = (
    ({0: ('C', 0, None), 1: ('C', 0, None)},
     {frozenset({0, 1}): (1, None)}),
    ({0: ('C', 0, None), 1: ('C', 0, None)},
     {frozenset({0, 1}): (2, None)}),
    ({0: ('C', 0, None), 1: ('C', 0, None)},
     {frozenset({0, 1}): (3, None)}),
)

CH2FH2H_CGR_IMP = (
    {0: ('F', 0, None), 1: ('C', 2, None), 2: ('H', 1, None),
     3: ('H', 0, None)},
    {frozenset({0, 1}): (1, None)})
CH2FH2H_CGR_EXP = (
    {0: ('F', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
     3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
     6: ('H', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({1, 4}): (1, None),
     frozenset({1, 5}): (1, None), frozenset({2, 6}): (1, None)})

C2H2CL2F2_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('F', 0, None),
     3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
     frozenset({1, 5}): (1, None)})
C2H2CL2F2_SGRS = (
    ({0: ('C', 1, False), 1: ('C', 1, False), 2: ('F', 0, None),
      3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
      frozenset({1, 5}): (1, None)}),
    ({0: ('C', 1, False), 1: ('C', 1, True), 2: ('F', 0, None),
      3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
      frozenset({1, 5}): (1, None)}),
    ({0: ('C', 1, True), 1: ('C', 1, False), 2: ('F', 0, None),
      3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
      frozenset({1, 5}): (1, None)}),
    ({0: ('C', 1, True), 1: ('C', 1, True), 2: ('F', 0, None),
      3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
      frozenset({1, 5}): (1, None)})
)

C3H3CL2F3_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
     3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
     6: ('F', 0, None), 7: ('F', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
     frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
     frozenset({2, 7}): (1, None)})
C3H3CL2F3_SGRS = (
    ({0: ('C', 1, None), 1: ('C', 1, False), 2: ('C', 1, False),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
    ({0: ('C', 1, None), 1: ('C', 1, True), 2: ('C', 1, True),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
    ({0: ('C', 1, False), 1: ('C', 1, False), 2: ('C', 1, True),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
    ({0: ('C', 1, False), 1: ('C', 1, True), 2: ('C', 1, False),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
    ({0: ('C', 1, True), 1: ('C', 1, False), 2: ('C', 1, True),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
    ({0: ('C', 1, True), 1: ('C', 1, True), 2: ('C', 1, False),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
)

C3H5N3_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
     3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
    {frozenset({1, 4}): (1, None), frozenset({1, 2}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({2, 5}): (1, None)})
C3H5N3_SGRS = (
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, False), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, None)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, False)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, True)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, False), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, True), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, False)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, False), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, True), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, True)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, True), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, None)}),
)

C8H13O_SGRS = (
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, False), 7: ('C', 1, False), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, False), 7: ('C', 1, False), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, True), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, False), 7: ('C', 1, True), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, False), 7: ('C', 1, True), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, True), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, True), 7: ('C', 1, False), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, True), 7: ('C', 1, False), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, True), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, True), 7: ('C', 1, True), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, True), 7: ('C', 1, True), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, True), frozenset({5, 7}): (1, None)}),
)

# FC=CC=CF + [OH] => FC=C[CH]C(O)F
C4H5F2O_TSG = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
                3: ('C', 0, None), 4: ('F', 0, None), 5: ('F', 0, None),
                6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None),
                9: ('H', 0, None), 10: ('O', 0, None), 11: ('H', 0, None)},
               {frozenset({8, 2}): (1, None), frozenset({2, 10}): (0.1, None),
                frozenset({0, 6}): (1, None), frozenset({1, 7}): (1, None),
                frozenset({9, 3}): (1, None), frozenset({0, 1}): (1, None),
                frozenset({0, 2}): (1, True), frozenset({2, 4}): (1, None),
                frozenset({3, 5}): (1, None), frozenset({10, 11}): (1, None),
                frozenset({1, 3}): (1, False)})

# FC=C(C(O)F)C(O)F + [OH] => FC(O)[C](C(O)F)C(O)F
C4H5F3O2_TSG = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, False),
                 3: ('C', 0, True), 4: ('F', 0, None), 5: ('F', 0, None),
                 6: ('F', 0, None), 7: ('O', 0, None), 8: ('O', 0, None),
                 9: ('H', 0, None), 10: ('H', 0, None), 11: ('H', 0, None),
                 12: ('H', 0, None), 13: ('H', 0, None), 14: ('O', 0, None),
                 15: ('H', 0, None)},
                {frozenset({12, 7}): (1, None), frozenset({2, 10}): (1, None),
                 frozenset({1, 2}): (1, None), frozenset({0, 1}): (1, True),
                 frozenset({3, 6}): (1, None), frozenset({2, 7}): (1, None),
                 frozenset({2, 5}): (1, None), frozenset({0, 4}): (1, None),
                 frozenset({8, 3}): (1, None), frozenset({0, 14}): (0.1, None),
                 frozenset({8, 13}): (1, None), frozenset({14, 15}): (1, None),
                 frozenset({11, 3}): (1, None), frozenset({1, 3}): (1, None),
                 frozenset({0, 9}): (1, None)})

# ISOBUTANE
C4H10_GRA = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
              3: ('C', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
              6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None),
              9: ('H', 0, None), 10: ('H', 0, None), 11: ('H', 0, None),
              12: ('H', 0, None), 13: ('H', 0, None)},
             {frozenset({0, 3}): (1, None), frozenset({0, 4}): (1, None),
              frozenset({0, 5}): (1, None), frozenset({0, 6}): (1, None),
              frozenset({1, 3}): (1, None), frozenset({1, 7}): (1, None),
              frozenset({8, 1}): (1, None), frozenset({1, 9}): (1, None),
              frozenset({2, 3}): (1, None), frozenset({2, 10}): (1, None),
              frozenset({2, 11}): (1, None), frozenset({2, 12}): (1, None),
              frozenset({3, 13}): (1, None)})


def test__from_data():
    """ test getters
    """
    cgr = automol.graph.from_data(
        atm_symb_dct=graph.atom_symbols(C8H13O_CGR),
        bnd_keys=graph.bond_keys(C8H13O_CGR),
        atm_imp_hyd_vlc_dct=(
            graph.atom_implicit_hydrogen_valences(C8H13O_CGR)),
    )
    assert cgr == C8H13O_CGR

    rgr = automol.graph.from_data(
        atm_symb_dct=graph.atom_symbols(C8H13O_RGR),
        bnd_keys=graph.bond_keys(C8H13O_RGR),
        atm_imp_hyd_vlc_dct=(
            graph.atom_implicit_hydrogen_valences(C8H13O_RGR)),
        bnd_ord_dct=graph.bond_orders(C8H13O_RGR),
    )
    assert rgr == C8H13O_RGR

    sgr = automol.graph.from_data(
        atm_symb_dct=graph.atom_symbols(C8H13O_SGR),
        bnd_keys=graph.bond_keys(C8H13O_SGR),
        atm_imp_hyd_vlc_dct=(
            graph.atom_implicit_hydrogen_valences(C8H13O_SGR)),
        atm_ste_par_dct=graph.atom_stereo_parities(C8H13O_SGR),
        bnd_ste_par_dct=graph.bond_stereo_parities(C8H13O_SGR)
    )
    assert sgr == C8H13O_SGR


def test__set_atom_implicit_hydrogen_valences():
    """ test graph.set_atom_implicit_hydrogen_valences
    """
    atm_keys = graph.atom_keys(C8H13O_CGR)
    cgr = graph.set_atom_implicit_hydrogen_valences(
        C8H13O_CGR, {atm_key: 0 for atm_key in atm_keys})

    assert cgr == automol.graph.from_data(
        graph.atom_symbols(C8H13O_CGR), graph.bond_keys(C8H13O_CGR))


def test__string():
    """ test graph.string and graph.from_string
    """
    for sgr in C8H13O_SGRS:
        assert sgr == automol.graph.from_string(automol.graph.string(sgr))


def test__without_bond_orders():
    """ test graph.without_bond_orders
    """
    assert C8H13O_CGR == graph.without_bond_orders(C8H13O_RGR)


def test__without_stereo_parities():
    """ test graph.without_stereo_parities
    """
    assert C8H13O_CGR == graph.without_stereo_parities(C8H13O_SGR)


def test__electron_count():
    """ test graph.electron_count
    """
    assert graph.electron_count(C8H13O_CGR) == 69


def test__atom_count():
    """ test graph.electron_count
    """
    assert graph.atom_count(C8H13O_CGR) == 22
    assert graph.atom_count(C8H13O_CGR, with_implicit=False) == 9


def test__heavy_atom_count():
    """ test graph.explicit_hydrogen_count
    """
    cgr = graph.explicit(C8H13O_CGR)
    assert graph.heavy_atom_count(cgr) == 9


def test__atoms_neighbor_atom_keys():
    """ test graph.atoms_neighbor_atom_keys
    """
    assert graph.atoms_neighbor_atom_keys(C8H13O_CGR) == {
        0: frozenset({3}),
        1: frozenset({4}),
        2: frozenset({6}),
        3: frozenset({0, 5}),
        4: frozenset({1, 6}),
        5: frozenset({3, 7}),
        6: frozenset({2, 4, 7}),
        7: frozenset({8, 5, 6}),
        8: frozenset({7})
    }


def test__atoms_second_degree_neighbor_atom_keys():
    """ test graph.atoms_neighbor_atom_keys
    """
    assert graph.atoms_second_degree_neighbor_atom_keys(C8H13O_CGR) == {
        0: frozenset({5}),
        1: frozenset({6}),
        2: frozenset({4, 7}),
        3: frozenset({7}),
        4: frozenset({2, 7}),
        5: frozenset({0, 8, 6}),
        6: frozenset({8, 1, 5}),
        7: frozenset({2, 3, 4}),
        8: frozenset({5, 6}),
    }


def test__atoms_bond_keys():
    """ test graph.atoms_neighbor_atom_keys
    """
    assert graph.atoms_bond_keys(C8H13O_CGR) == {
        0: frozenset({frozenset({0, 3})}),
        1: frozenset({frozenset({1, 4})}),
        2: frozenset({frozenset({2, 6})}),
        3: frozenset({frozenset({3, 5}), frozenset({0, 3})}),
        4: frozenset({frozenset({1, 4}), frozenset({4, 6})}),
        5: frozenset({frozenset({3, 5}), frozenset({5, 7})}),
        6: frozenset({frozenset({6, 7}), frozenset({4, 6}),
                      frozenset({2, 6})}),
        7: frozenset({frozenset({6, 7}), frozenset({5, 7}),
                      frozenset({8, 7})}),
        8: frozenset({frozenset({8, 7})})
    }


# # bond properties
def test__bonds_neighbor_atom_keys():
    """ test graph.bonds_neighbor_atom_keys
    """

    assert graph.bonds_neighbor_atom_keys(C8H13O_CGR) == {
        frozenset({1, 4}): frozenset({6}),
        frozenset({4, 6}): frozenset({1, 2, 7}),
        frozenset({2, 6}): frozenset({4, 7}),
        frozenset({0, 3}): frozenset({5}),
        frozenset({6, 7}): frozenset({8, 2, 4, 5}),
        frozenset({8, 7}): frozenset({5, 6}),
        frozenset({3, 5}): frozenset({0, 7}),
        frozenset({5, 7}): frozenset({8, 3, 6})
    }


# # other properties
def test__branch():
    """ test graph.branch
    """
    assert graph.branch(C8H13O_CGR, 6, frozenset({6, 4})) == (
        {1: ('C', 2, None), 4: ('C', 1, None), 6: ('C', 1, None)},
        {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None)}
    )


def test__connected_components():
    """ test graph.connected_components
    """
    gra1 = C3H3_CGR
    gra2 = C2_CGR
    gra1_natms = automol.formula.atom_count(graph.formula(C3H3_CGR))
    gra2 = graph.transform_keys(gra2, lambda x: x + gra1_natms)

    gra = graph.union(gra1, gra2)
    cmp_gras = graph.connected_components(gra)
    assert cmp_gras in [(gra1, gra2), (gra2, gra1)]


def test__subgraph():
    """ test graph.subgraph
    """
    assert graph.subgraph(C3H3_CGR, (1, 2)) == (
        {1: ('C', 1, None), 2: ('C', 1, None)},
        {frozenset({1, 2}): (1, None)})


def test__bond_induced_subgraph():
    """ test graph.bond_induced_subgraph
    """
    assert graph.bond_induced_subgraph(
        C3H3_CGR, [frozenset({0, 1}), frozenset({1, 2})]) == (
            {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
            {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})


# # transformations
def test__relabel():
    """ test graph.relabel
    """
    assert graph.relabel(C3H3_CGR, {0: 10, 1: 11, 2: 12}) == (
        {10: ('C', 1, None), 11: ('C', 1, None), 12: ('C', 1, None)},
        {frozenset({10, 11}): (1, None), frozenset({11, 12}): (1, None),
         frozenset({12, 10}): (1, None)})


def test__remove_atoms():
    """ test graph.remove_atoms
    """
    assert graph.remove_atoms(C3H3_CGR, (0,)) == (
        {1: ('C', 1, None), 2: ('C', 1, None)},
        {frozenset({1, 2}): (1, None)})


def test__remove_bonds():
    """ test graph.remove_bonds
    """
    assert graph.remove_bonds(C3H3_CGR, [frozenset({1, 2})]) == (
        {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
        {frozenset({0, 1}): (1, None), frozenset({2, 0}): (1, None)})


# implicit/explicit hydrogen functions
# # atom properties
def test__atom_explicit_hydrogen_valences():
    """ test graph.atom_explicit_hydrogen_valences
    """
    assert graph.atom_explicit_hydrogen_valences(CH2FH2H_CGR_EXP) == {
        0: 0, 1: 2, 2: 1, 3: 0, 4: 0, 5: 0, 6: 0
    }


def test__atom_explicit_hydrogen_keys():
    """ test graph.atom_explicit_hydrogen_keys
    """
    assert graph.atom_explicit_hydrogen_keys(CH2FH2H_CGR_EXP) == {
        0: frozenset(),
        1: frozenset({4, 5}),
        2: frozenset({6}),
        3: frozenset(),
        4: frozenset(),
        5: frozenset(),
        6: frozenset()
    }


# # other properties
def test__backbone_keys():
    """ test graph.backbone_keys
    """
    assert graph.backbone_keys(CH2FH2H_CGR_EXP) == frozenset({0, 1, 2, 3})


def test__explicit_hydrogen_keys():
    """ test graph.explicit_hydrogen_keys
    """
    assert graph.explicit_hydrogen_keys(CH2FH2H_CGR_EXP) == frozenset(
        {4, 5, 6})


def test__explicit():
    """ test graph.explicit
    """
    assert CH2FH2H_CGR_EXP == graph.explicit(CH2FH2H_CGR_IMP)


def test__implicit():
    """ test graph.implicit
    """
    assert CH2FH2H_CGR_IMP == graph.implicit(graph.explicit(CH2FH2H_CGR_IMP))


# # comparisons
def test__backbone_isomorphic():
    """ test graph.backbone_isomorphic
    """
    assert graph.backbone_isomorphic(CH2FH2H_CGR_IMP, CH2FH2H_CGR_EXP)

    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        assert graph.backbone_isomorphic(cgr, cgr_pmt)


def test__backbone_isomorphism():
    """ test graph.backbone_isomorphism
    """
    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        assert graph.backbone_isomorphism(cgr, cgr_pmt) == pmt_dct


def test__backbone_unique():
    """ test graph.backbone_unique
    """
    assert graph.backbone_unique(C3H3_RGRS) == C3H3_RGRS[:2]


# chemistry library
def test__atom_element_valences():
    """ test graph.atom_element_valences
    """
    assert graph.atom_element_valences(C8H13O_CGR) == {
        0: 4, 1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4, 7: 4, 8: 2}


def test__atom_lone_pair_counts():
    """ test graph.atom_lone_pair_counts
    """
    assert graph.atom_lone_pair_counts(C8H13O_CGR) == {
        0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 2}


def test__atom_bond_valences():
    """ test graph.atom_bond_valences
    """
    assert graph.atom_bond_valences(C8H13O_CGR) == {
        0: 4, 1: 3, 2: 4, 3: 3, 4: 3, 5: 3, 6: 4, 7: 4, 8: 1}


def test__atom_unsaturated_valences():
    """ test graph.atom_unsaturated_valences
    """
    assert graph.atom_unsaturated_valences(C8H13O_CGR) == {
        0: 0, 1: 1, 2: 0, 3: 1, 4: 1, 5: 1, 6: 0, 7: 0, 8: 1}


def test__unsaturated_atom_keys():
    """ test graph.unsaturated_atom_keys
    """
    assert graph.unsaturated_atom_keys(C8H13O_CGR) == frozenset(
        {1, 3, 4, 5, 8})


def test__maximum_spin_multiplicity():
    """ test graph.maximum_spin_multiplicity
    """
    assert graph.maximum_spin_multiplicity(C2_CGR) == 7


def test__possible_spin_multiplicities():
    """ test graph.possible_spin_multiplicities
    """
    assert graph.possible_spin_multiplicities(C2_CGR) == (1, 3, 5, 7)


# miscellaneous
def test__bond_symmetry_numbers():
    """ test graph.bond_symmetry_numbers
    """
    assert graph.bond_symmetry_numbers(C8H13O_CGR) == {
        frozenset({1, 4}): 1, frozenset({4, 6}): 1, frozenset({2, 6}): 3,
        frozenset({0, 3}): 3, frozenset({6, 7}): 1, frozenset({8, 7}): 1,
        frozenset({3, 5}): 1, frozenset({5, 7}): 1}


# resonance graph library
# # atom properties
def test__resonance_dominant_atom_hybridizations():
    """ test graph.resonance_dominant_atom_hybridizations
    """
    assert graph.resonance_dominant_atom_hybridizations(C3H3_CGR) == {
        0: 2, 1: 2, 2: 2}
    assert graph.resonance_dominant_atom_hybridizations(C8H13O_CGR) == {
        0: 3, 1: 2, 2: 3, 3: 2, 4: 2, 5: 2, 6: 3, 7: 3, 8: 3}

    cgr = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('O', 0, None),
            3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
            6: ('X', 0, None)},
           {frozenset({1, 4}): (1, None), frozenset({1, 2}): (1, None),
            frozenset({0, 3}): (1, None), frozenset({0, 1}): (1, None),
            frozenset({2, 5}): (1, None)})
    assert graph.resonance_dominant_atom_hybridizations(cgr) == {
        0: 2, 1: 2, 2: 3, 3: 0, 4: 0, 5: 0, 6: -1}


def test__resonance_dominant_atom_centered_cumulene_keys():
    """ test graph.resonance_dominant_atom_centered_cumulene_keys
    """
    cgr = ({0: ('C', 1, None), 1: ('C', 2, None), 2: ('C', 0, None),
            3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 0, None),
            6: ('C', 0, None)},
           {frozenset({4, 6}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({2, 4}): (1, None), frozenset({5, 6}): (1, None),
            frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})
    assert (graph.resonance_dominant_atom_centered_cumulene_keys(cgr) ==
            frozenset({(frozenset({1, 4}), 5)}))


def test__resonance_dominant_bond_centered_cumulene_keys():
    """ test graph.resonance_dominant_bond_centered_cumulene_keys
    """
    cgr = ({0: ('C', 1, None), 1: ('C', 2, None), 2: ('C', 0, None),
            3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 0, None)},
           {frozenset({4, 5}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
            frozenset({1, 3}): (1, None)})
    assert (graph.resonance_dominant_bond_centered_cumulene_keys(cgr) ==
            frozenset({(frozenset({1, 4}), frozenset({3, 5}))}))


def test__resonance_dominant_radical_atom_keys():
    """ test graph.resonance_dominant_radical_atom_keys
    """
    assert graph.resonance_dominant_radical_atom_keys(C3H3_CGR) == frozenset(
        {0, 1, 2})
    assert graph.resonance_dominant_radical_atom_keys(C8H13O_CGR) == frozenset(
        {8})


def test__sigma_radical_atom_keys():
    """ test graph.sigma_radical_atom_keys
    """
    # CCC#[C]
    gra = ({0: ('C', 3, None), 1: ('C', 0, None), 2: ('C', 2, None),
            3: ('C', 0, None)},
           {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
            frozenset({2, 3}): (1, None)})
    assert graph.sigma_radical_atom_keys(gra) == frozenset({1})

    # [C]#CC(CC)(CCC#[C])CC#[C]
    gra = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 3, None),
            3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 0, None),
            6: ('C', 2, None), 7: ('C', 0, None), 8: ('C', 2, None),
            9: ('C', 2, None), 10: ('C', 2, None), 11: ('C', 0, None)},
           {frozenset({8, 4}): (1, None), frozenset({3, 7}): (1, None),
            frozenset({2, 6}): (1, None), frozenset({0, 4}): (1, None),
            frozenset({8, 10}): (1, None), frozenset({9, 11}): (1, None),
            frozenset({1, 5}): (1, None), frozenset({9, 5}): (1, None),
            frozenset({11, 7}): (1, None), frozenset({10, 11}): (1, None),
            frozenset({11, 6}): (1, None)})
    assert graph.sigma_radical_atom_keys(gra) == frozenset({0, 1, 3})


# # bond properties
def test__resonance_dominant_bond_orders():
    """ test graph.resonance_dominant_bond_orders
    """
    assert graph.resonance_dominant_bond_orders(C3H3_CGR) == {
        frozenset({0, 1}): frozenset({1, 2}),
        frozenset({0, 2}): frozenset({1, 2}),
        frozenset({1, 2}): frozenset({1, 2})
    }


# # transformations
def test__resonances():
    """ test graph.resonances
    """
    assert graph.resonances(C3H3_CGR) == C3H3_RGRS


def test__subresonances():
    """ test graph.subresonances
    """
    assert graph.subresonances(C2_RGRS[1]) == C2_RGRS[1:]


def test__dominant_resonances():
    """ test graph.dominant_resonances
    """
    assert graph.dominant_resonances(C3H3_CGR) == C3H3_RGRS[1:]


def test__dominant_resonance():
    """ test graph.dominant_resonance
    """
    assert graph.dominant_resonance(C3H3_CGR) == C3H3_RGRS[1]


def test__rotational_bond_keys():
    """ test graph.rotational_bond_keys
    """
    cgr = ({0: ('C', 2, None), 1: ('C', 2, None), 2: ('C', 1, None),
            3: ('C', 1, None)},
           {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
            frozenset({2, 3}): (1, None)})
    cgr = automol.graph.explicit(cgr)
    assert (automol.graph.rotational_bond_keys(cgr) ==
            frozenset({frozenset({2, 3})}))

    cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 2, None),
            3: ('C', 2, None)},
           {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
            frozenset({2, 3}): (1, None)})
    assert (automol.graph.rotational_bond_keys(cgr) ==
            frozenset({frozenset({0, 2}), frozenset({1, 3}),
                       frozenset({2, 3})}))
    assert (automol.graph.rotational_bond_keys(cgr, with_h_rotors=False) ==
            frozenset({frozenset({2, 3})}))
    assert (automol.graph.rotational_bond_keys(cgr, with_chx_rotors=False) ==
            frozenset({frozenset({2, 3})}))


def test__species__graph_conversion():
    """ test interchanging between graphs aligned by zma and geo
    """

    ich = automol.smiles.inchi('CC#CC#CCCCC#CC')
    geo = automol.inchi.geometry(ich)
    gra = automol.geom.graph(geo)

    zma, zma_keys, dummy_key_dct = (
        automol.geom.zmatrix_with_conversion_info(geo))
    zgra = automol.graph.relabel_for_zmatrix(gra, zma_keys, dummy_key_dct)

    geo, gdummy_key_dct = automol.zmat.geometry_with_conversion_info(zma)
    ggra = automol.graph.relabel_for_geometry(zgra)

    old_zgra = zgra
    zgra = automol.graph.insert_dummy_atoms(ggra, gdummy_key_dct)
    assert zgra == old_zgra


# stereo graph library
def test__stereogenic_atom_keys():
    """ test graph.stereogenic_atom_keys
    """
    assert graph.stereogenic_atom_keys(C8H13O_CGR) == frozenset({6, 7})
    assert graph.stereogenic_atom_keys(C3H3CL2F3_CGR) == frozenset({1, 2})

    cgr = ({0: ('C', 2, None), 1: ('C', 3, None), 2: ('C', 1, None),
            3: ('O', 1, None)},
           {frozenset({0, 2}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({1, 2}): (1, None)})
    print(graph.stereogenic_atom_keys(cgr))
    assert graph.stereogenic_atom_keys(cgr) == frozenset({2})

    # Bug fix:
    cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 2, None),
            3: ('C', 2, None), 4: ('C', 2, None), 5: ('C', 2, None),
            6: ('C', 2, None), 7: ('C', 2, None), 8: ('C', 2, None),
            9: ('C', 2, None), 10: ('C', 1, None), 11: ('C', 1, None),
            12: ('O', 0, None)},
           {frozenset({4, 6}): (1, None), frozenset({11, 12}): (1, None),
            frozenset({8, 9}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({10, 6}): (1, None), frozenset({8, 10}): (1, None),
            frozenset({2, 4}): (1, None), frozenset({9, 11}): (1, None),
            frozenset({10, 12}): (1, None), frozenset({3, 5}): (1, None),
            frozenset({11, 7}): (1, None), frozenset({1, 3}): (1, None),
            frozenset({5, 7}): (1, None)})
    print(graph.stereogenic_atom_keys(cgr))
    assert graph.stereogenic_atom_keys(cgr) == frozenset({10, 11})


def test__stereogenic_bond_keys():
    """ test graph.stereogenic_bond_keys
    """
    assert graph.stereogenic_bond_keys(C8H13O_CGR) == frozenset(
        {frozenset({3, 5})})
    assert graph.stereogenic_bond_keys(C3H5N3_CGR) == frozenset(
        {frozenset({1, 4}), frozenset({0, 3})})


def test__stereomers():
    """ test graph.stereomers
    """
    assert graph.stereomers(C2H2CL2F2_CGR) == C2H2CL2F2_SGRS
    assert graph.stereomers(C3H3CL2F3_CGR) == C3H3CL2F3_SGRS
    assert graph.stereomers(C3H5N3_CGR) == C3H5N3_SGRS
    assert graph.stereomers(C8H13O_CGR) == C8H13O_SGRS


def test__ring_systems():
    """ test graph.ring_systems
    """
    ich = automol.smiles.inchi('C12CC(C1)C2CC3C(C3)CCC4C5CCC(CC5)C4')
    gra = automol.inchi.graph(ich)
    rsys = sorted(graph.ring_systems(gra), key=graph.atom_count)
    assert list(map(graph.atom_count, rsys)) == [7, 12, 21]


def test__equivalent_atoms():
    """ test graph.equivalent_atoms
    """
    # central carbon
    assert graph.equivalent_atoms(C4H10_GRA, 3) == {3}
    # central hydrogen
    assert graph.equivalent_atoms(C4H10_GRA, 13) == {13}
    # terminal carbons
    assert graph.equivalent_atoms(C4H10_GRA, 0) == {0, 1, 2}
    assert graph.equivalent_atoms(C4H10_GRA, 1) == {0, 1, 2}
    assert graph.equivalent_atoms(C4H10_GRA, 2) == {0, 1, 2}
    # terminal hydrogens
    assert graph.equivalent_atoms(C4H10_GRA, 4) == {4, 5, 6, 7, 8, 9, 10,
                                                    11, 12}
    assert graph.equivalent_atoms(C4H10_GRA, 5) == {4, 5, 6, 7, 8, 9, 10,
                                                    11, 12}
    assert graph.equivalent_atoms(C4H10_GRA, 6) == {4, 5, 6, 7, 8, 9, 10,
                                                    11, 12}
    assert graph.equivalent_atoms(C4H10_GRA, 11) == {4, 5, 6, 7, 8, 9, 10,
                                                     11, 12}
    assert graph.equivalent_atoms(C4H10_GRA, 12) == {4, 5, 6, 7, 8, 9, 10,
                                                     11, 12}


def test__equivalent_bonds():
    """ test graph.equivalent_atoms
    """
    assert graph.equivalent_bonds(C4H10_GRA, (2, 3)) == {
        (0, 3), (1, 3), (2, 3)}


def test__vmat__vmatrix():
    """ test graph.vmat.vmatrix
    """
    ich = automol.smiles.inchi('C12CC(C1)C2CC3C(C3)CCC4C5CCC(CC5)C4')
    gra = automol.inchi.graph(ich)
    _, zma_keys = graph.vmat.vmatrix(gra)
    assert set(zma_keys) == graph.atom_keys(gra)


def test__ts__nonconserved_atom_stereo_keys():
    """ test graph.ts.nonconserved_atom_stereo_keys
    """
    assert graph.ts.nonconserved_atom_stereo_keys(C4H5F2O_TSG) == (
        (frozenset({2}), frozenset()))
    assert graph.ts.nonconserved_atom_stereo_keys(C4H5F3O2_TSG) == (
        (frozenset({0}), frozenset()))


def test__ts__nonconserved_bond_stereo_keys():
    """ test graph.ts.nonconserved_bond_stereo_keys
    """
    assert graph.ts.nonconserved_bond_stereo_keys(C4H5F2O_TSG) == (
        (frozenset({frozenset({0, 1})}), frozenset({frozenset({0, 2})})))
    assert graph.ts.nonconserved_bond_stereo_keys(C4H5F3O2_TSG) == (
        (frozenset(), frozenset({frozenset({0, 1})})))


def test__ts__compatible_reverse_stereomers():
    """ test graph.ts.stereo_expand_reverse_graphs
    """
    for ste_tsg in graph.ts.stereomers(C4H5F2O_TSG):
        ste_tsgs = [
            s
            for r in graph.ts.compatible_reverse_stereomers(ste_tsg)
            for s in graph.ts.compatible_reverse_stereomers(r)]
        assert any(s == ste_tsg for s in ste_tsgs)

    for ste_tsg in graph.ts.stereomers(C4H5F3O2_TSG):
        ste_tsgs = [
            s
            for r in graph.ts.compatible_reverse_stereomers(ste_tsg)
            for s in graph.ts.compatible_reverse_stereomers(r)]
        assert any(s == ste_tsg for s in ste_tsgs)


def test__canonical():
    """ test graph.canonical
    """

    def _test_from_smiles(smi):
        ich = automol.smiles.inchi(smi)
        print(ich)
        geo = automol.inchi.geometry(ich)
        gra = automol.geom.graph(geo)

        gra = automol.graph.implicit(gra)

        can_gra = automol.graph.canonical(gra)
        print(automol.graph.string(can_gra, one_indexed=True))

        natms = len(automol.graph.atom_keys(gra))

        for _ in range(10):
            pmt = list(map(int, numpy.random.permutation(natms)))
            print(pmt)
            pmt_gra = automol.graph.relabel(gra, dict(enumerate(pmt)))
            can_pmt_gra = automol.graph.canonical(pmt_gra)
            print(automol.graph.string(can_pmt_gra, one_indexed=True))
            assert can_pmt_gra == can_gra

    # More tests can be added here
    smis = [
        'c1ccccc1CC',
        'F[C@@H]([C@@H](F)Cl)[C@H](F)Cl',
        r'[H]/N=C\C(\C=N\[H])=N\[H]',
    ]
    for smi in smis:
        _test_from_smiles(smi)


def test__canonical_priorities_and_stereo_parities():
    """ test graph.canonical_priorities_and_stereo_parities
    """

    def _test_from_smiles(smi, ref_atm_pars, ref_bnd_pars):
        ich = automol.smiles.inchi(smi)
        print(ich)
        geo = automol.inchi.geometry(ich)
        gra = automol.geom.graph(geo)

        par_eval_ = graph.parity_evaluator_from_geometry_(gra, geo)

        can_key_dct, atm_par_dct, bnd_par_dct = (
            graph.canonical_priorities_and_stereo_parities(
                gra, par_eval_=par_eval_))
        print(can_key_dct)

        atm_pars = [p for k, p in sorted(atm_par_dct.items()) if p is not None]
        bnd_pars = [p for k, p in sorted(bnd_par_dct.items()) if p is not None]

        print(atm_pars, bnd_pars)
        assert atm_pars == ref_atm_pars
        assert bnd_pars == ref_bnd_pars

    # More tests can be added here
    args = [
        # Atom stereo tests
        ('C[C@@](F)(N)O', [True], []),  # R = clock  = '+' => True
        ('C[C@](F)(N)O', [False], []),  # S = aclock = '-' => False
        ('[C@@H](F)(N)O', [True], []),  # R = clock  = '+' => True
        ('[C@H](F)(N)O', [False], []),  # S = aclock = '-' => False
        # Bond stereo tests
        (r'F/C=C/F', [], [True]),           # trans = '+' => True
        (r'F/C=C\F', [], [False]),          # cis   = '-' => False
        (r'[H]/N=C(Cl)/F', [], [True]),     # trans = '+' => True
        (r'[H]/N=C(Cl)\F', [], [False]),    # cis   = '-' => False
        (r'[H]/N=C/F', [], [True]),         # trans = '+' => True
        (r'[H]/N=C\F', [], [False]),        # cis   = '-' => False
        (r'[H]/N=N/[H]', [], [True]),       # trans = '+' => True
        (r'[H]/N=N\[H]', [], [False]),      # cis   = '-' => False
        # Advanced tests
        ('F[C@@H]([C@@H](F)Cl)[C@H](F)Cl', [True, True, False], []),
        (r'[H]/N=C\C(\C=N\[H])=N\[H]', [], [True, False, False]),
    ]
    for smi, ref_atm_pars, ref_bnd_pars in args:
        _test_from_smiles(smi, ref_atm_pars, ref_bnd_pars)


def test__to_local_stereo():
    """ test graph.to_local_stereo
    """
    # Atom parity test:
    # Indices 3 and 4 are swapped so that local stereo will have opposite
    # parity
    can_gra = ({0: ('C', 3, None), 1: ('C', 0, True), 2: ('F', 0, None),
                4: ('N', 2, None), 3: ('O', 1, None)},
               {frozenset({0, 1}): (1, None), frozenset({1, 4}): (1, None),
                frozenset({1, 3}): (1, None), frozenset({1, 2}): (1, None)})
    can_par = graph.atom_stereo_parities(can_gra)[1]

    loc_gra = graph.to_local_stereo(can_gra)
    loc_par = graph.atom_stereo_parities(loc_gra)[1]
    print(can_par, loc_par)
    assert can_par is True
    assert loc_par is False

    # Bond parity test:
    # Indices 3 and 5 are swapped so that local stereo will have opposite
    # parity
    can_gra = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('Cl', 0, None),
                5: ('Cl', 0, None), 4: ('F', 0, None), 3: ('F', 0, None)},
               {frozenset({0, 1}): (1, True), frozenset({0, 2}): (1, None),
                frozenset({0, 4}): (1, None), frozenset({1, 3}): (1, None),
                frozenset({1, 5}): (1, None)})
    can_par = graph.bond_stereo_parities(can_gra)[frozenset({0, 1})]

    loc_gra = graph.to_local_stereo(can_gra)
    loc_par = graph.bond_stereo_parities(loc_gra)[frozenset({0, 1})]
    print(can_par, loc_par)
    assert can_par is True
    assert loc_par is False


def test__from_local_stereo():
    """ test graph.from_local_stereo
    """
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
    """ test graph.has_resonance_bond_stereo
    """
    gra = ({0: ('F', 0, None), 1: ('C', 1, None), 3: ('C', 1, None),
            4: ('C', 1, None), 5: ('F', 0, None)},
           {frozenset({3, 4}): (1, True), frozenset({0, 1}): (1, None),
            frozenset({1, 3}): (1, True), frozenset({4, 5}): (1, None)})
    assert graph.has_resonance_bond_stereo(gra)

    gra = ({0: ('F', 0, None), 1: ('C', 2, None), 4: ('C', 1, None),
            5: ('C', 1, None), 6: ('F', 0, None)},
           {frozenset({4, 5}): (1, True), frozenset({0, 1}): (1, None),
            frozenset({1, 4}): (1, None), frozenset({5, 6}): (1, None)})
    assert not graph.has_resonance_bond_stereo(gra)


def test__amchi():
    """ test graph.amchi
    """
    # bond stereo
    gra = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
            3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
           {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
            frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
            frozenset({2, 5}): (1, False)})
    chi = graph.amchi(gra)
    print(chi)

    assert chi == 'AMChI=1/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3-'

    # atom stereo
    gra = ({0: ('C', 1, None), 1: ('C', 1, True), 2: ('C', 1, True),
            3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
            6: ('F', 0, None), 7: ('F', 0, None)},
           {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
            frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
            frozenset({2, 7}): (1, None)})
    chi = graph.amchi(gra)
    print(chi)

    assert chi == 'AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t2-,3-/m1/s1'


def test__amchi_with_indices():
    """ test graph.amchi
    """
    # bond stereo
    gra1 = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
             3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
            {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
             frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
             frozenset({2, 5}): (1, False)})
    chi1, chi1_idx_dcts = graph.amchi_with_indices(gra1)
    print(chi1)
    print(chi1_idx_dcts)

    assert chi1 == 'AMChI=1/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3-'

    # atom stereo
    gra2 = ({0: ('C', 1, None), 1: ('C', 1, True), 2: ('C', 1, True),
             3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
             6: ('F', 0, None), 7: ('F', 0, None)},
            {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
             frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
             frozenset({2, 7}): (1, None)})
    chi2, chi2_idx_dcts = graph.amchi_with_indices(gra2)
    print(chi2)
    print(chi2_idx_dcts)

    assert chi2 == 'AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t2-,3-/m1/s1'

    gra = automol.graph.union_from_sequence([gra1, gra2], shift_keys=True)
    chi, chi_idx_dcts = graph.amchi_with_indices(gra)
    print(chi)
    print(chi_idx_dcts)

    assert chi == ('AMChI=1/C3H3Cl2F3.C3H5N3/c4-2(7)1(6)3(5)8;4-1-3(6)2-5/'
                   'h1-3H;1-2,4-6H/b;4-1-,5-2+,6-3-/t2-,3-;/m1./s1')
    assert chi_idx_dcts == (
        {6: 0, 7: 1, 8: 2, 9: 3, 10: 4, 11: 5, 12: 6, 13: 7},
        {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5})


def test__smiles():
    """ test graph.smiles
    """
    smis = [
        # Rings:
        'C1[C@@]2(C3)O[C@@]2(C[N@H]3)OC1',
        r'[CH2]/C=C1\[N@H]C(=C1)O[O]',
        'C1CC1',
        # Stereo atoms:
        'N[C@](C)(F)C(=O)O',
        'N[C@H](C)C(=O)O',
        'C[C@H]1CCCCO1',
        # Stereo bonds:
        r'F/C=C/F',
        r'F/C=C\F',
        r'CC/C=C/C=C/CF',
        r'C1CCCCCCCCCC/N=N/1',
        r'[H]/N=N/[H]',
        r'[H]/N=N/N=N\[H]',
    ]

    for smi in smis:
        print()
        print("STARTING SMILES:", smi)
        ich = automol.smiles.inchi(smi)
        geo = automol.inchi.geometry(ich)
        gra = automol.geom.graph(geo)
        print(automol.graph.string(gra))
        print(automol.geom.string(geo))
        smi = automol.graph.smiles(gra)
        print('smiles from code:', smi)

        ich_smi = automol.inchi.smiles(ich)
        print('inchi:', ich)
        print('smiles from inchi:', ich_smi)

        sich = automol.smiles.inchi(smi)
        print('inchi from smiles:', sich)
        # print(automol.graph.string(GRA, one_indexed=True))
        print(automol.graph.rings_atom_keys(gra))
        assert sich == ich


def test__smiles__with_resonance():
    """ test graph.smiles
    """

    smis = [
        r'FC=C-C=C-[CH]O',
        r'F[CH]C=CF',
    ]
    for smi in smis:
        print()
        print("STARTING SMILES:", smi)
        ich = automol.smiles.inchi(smi)
        geo = automol.inchi.geometry(ich)
        gra = automol.geom.graph(geo)
        print(automol.graph.string(gra))
        print(automol.geom.string(geo))

        ste_keys = automol.graph.bond_stereo_keys(gra)
        nste = len(ste_keys)
        ste_pars_lst = list(
            itertools.product(map(bool, range(2)), repeat=nste))
        for ste_pars in ste_pars_lst:
            bnd_par_dct = dict(zip(ste_keys, ste_pars))
            print(bnd_par_dct)
            gra = automol.graph.set_bond_stereo_parities(gra, bnd_par_dct)
            ach = automol.graph.amchi(gra)
            print('amchi:', ach)
            smi = automol.graph.smiles(gra)
            print('smiles from code:', smi)
            print()


if __name__ == '__main__':
    # test__smiles()
    # test__smiles()
    # test__canonical()
    test__canonical_priorities_and_stereo_parities()
    # test__to_local_stereo()

    # test__has_resonance_bond_stereo()
    # test__amchi_with_indices()
    # test__stereogenic_atom_keys()
    # test__ts__nonconserved_atom_stereo_keys()
    # test__ts__compatible_reverse_stereomers()
    # test__stereogenic_atom_keys()
