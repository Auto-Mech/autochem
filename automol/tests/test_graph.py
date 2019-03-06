""" test the automechanc.mol.graph module
"""
import numpy
from automol import graph

C_ICH = 'InChI=1S/C'
C_CGR = (
    {0: ('C', 0, None)}, {})

C2_ICH = 'InChI=1S/C2/c1-2'
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

C3H3_ICH = 'InChI=1S/C3H3/c1-2-3-1/h1-3H'
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

CH2FH2H_ICH = 'InChI=1S/CH2F.H2.H/c1-2;;/h1H2;1H;'
CH2FH2H_CGR = (
    {0: ('F', 0, None), 1: ('C', 2, None), 2: ('H', 1, None),
     3: ('H', 0, None)},
    {frozenset({0, 1}): (1, None)})

CH2FH2H_CGR_EXP = (
    {0: ('F', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
     3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
     6: ('H', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({1, 4}): (1, None),
     frozenset({1, 5}): (1, None), frozenset({2, 6}): (1, None)})

C5H5N5O_ICH = (
    'InChI=1S/C5H5N5O/c6-3-2-4(8-1-7-2)10-5(11)9-3/h1H,(H4,6,7,8,9,10,11)')
C5H5N5O_CGR = (
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

C2H2CL2F2_ICH = 'InChI=1S/C2H2Cl2F2/c3-1(5)2(4)6/h1-2H'
C2H2CL2F2_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('F', 0, None),
     3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
     frozenset({1, 5}): (1, None)})
C2H2CL2F2_STE_ICHS = (
    'InChI=1S/C2H2Cl2F2/c3-1(5)2(4)6/h1-2H/t1-,2-/m0/s1',
    'InChI=1S/C2H2Cl2F2/c3-1(5)2(4)6/h1-2H/t1-,2+',
    'InChI=1S/C2H2Cl2F2/c3-1(5)2(4)6/h1-2H/t1-,2+',
    'InChI=1S/C2H2Cl2F2/c3-1(5)2(4)6/h1-2H/t1-,2-/m1/s1'
)
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
      frozenset({1, 5}): (1, None)}),
)

C3H3CL2F3_ICH = 'InChI=1S/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H'
C3H3CL2F3_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
     3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
     6: ('F', 0, None), 7: ('F', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
     frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
     frozenset({2, 7}): (1, None)})
# these are incorrect -- currently we can't handle InChI generation for
# higher-order stereo
# C3H3CL2F3_STE_ICHS = (
#     'InChI=1S/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t2-,3-/m0/s1',
#     'InChI=1S/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t2-,3-/m1/s1',
#     'InChI=1S/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t1-,2-,3+',
#     'InChI=1S/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t1-,2+,3-',
#     'InChI=1S/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t1-,2-,3+',
#     'InChI=1S/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t1-,2+,3-',
# )
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

C3H5N3_ICH = 'InChI=1S/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H'
C3H5N3_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
     3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
    {frozenset({1, 4}): (1, None), frozenset({1, 2}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({2, 5}): (1, None)})
# these are incorrect -- currently we can't handle InChI generation for
# higher-order stereo
# C3H5N3_STE_ICHS = (
#     'InChI=1S/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1+,5-2+',
#     'InChI=1S/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3+',
#     'InChI=1S/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3-',
#     'InChI=1S/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3+',
#     'InChI=1S/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3-',
#     'InChI=1S/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2-',
# )
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

C8H13O_ICH = 'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3'
C8H13O_CGR = (
    {0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 1, None), 7: ('C', 1, None), 8: ('O', 0, None)},
    {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
     frozenset({3, 5}): (1, None), frozenset({5, 7}): (1, None)})
C8H13O_RGRS = (
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, None), 7: ('C', 1, None), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, None), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, None), 7: ('C', 1, None), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (2, None), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, None), 7: ('C', 1, None), 8: ('O', 0, None)},
     {frozenset({1, 4}): (2, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, None), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, None), 7: ('C', 1, None), 8: ('O', 0, None)},
     {frozenset({1, 4}): (2, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (2, None), frozenset({5, 7}): (1, None)}),
)
C8H13O_STE_ICHS = (
    'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3/b6-4-/t7-,8-/m0/s1',
    'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3/b6-4+/t7-,8-/m0/s1',
    'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3/b6-4-/t7-,8+/m0/s1',
    'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3/b6-4+/t7-,8+/m0/s1',
    'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3/b6-4-/t7-,8+/m1/s1',
    'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3/b6-4+/t7-,8+/m1/s1',
    'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3/b6-4-/t7-,8-/m1/s1',
    'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3/b6-4+/t7-,8-/m1/s1'
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


# core library
# # constructors and value getters
def test__from_atoms_and_bonds():
    """ test graph.from_atoms_and_bonds
    """
    assert C8H13O_CGR == graph.from_atoms_and_bonds(
        atm_dct=graph.atoms(C8H13O_CGR), bnd_dct=graph.bonds(C8H13O_CGR))


def test__from_dictionaries():
    """ test graph.from_dictionaries
    """
    assert graph.from_dictionaries(
        graph.atom_symbols(CH2FH2H_CGR_EXP),
        graph.bond_keys(CH2FH2H_CGR_EXP)
    ) == CH2FH2H_CGR_EXP

    assert graph.from_dictionaries(
        graph.atom_symbols(C8H13O_RGRS[0]),
        graph.bond_keys(C8H13O_RGRS[0]),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(
            C8H13O_RGRS[0]),
        bnd_ord_dct=graph.bond_orders(C8H13O_RGRS[0])
    ) == C8H13O_RGRS[0]

    assert graph.from_dictionaries(
        graph.atom_symbols(C8H13O_SGRS[0]),
        graph.bond_keys(C8H13O_SGRS[0]),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(
            C8H13O_SGRS[0]),
        atm_ste_par_dct=graph.atom_stereo_parities(C8H13O_SGRS[0]),
        bnd_ste_par_dct=graph.bond_stereo_parities(C8H13O_SGRS[0])
    ) == C8H13O_SGRS[0]

    # a litte ridiculous, but make sure we get the keys right
    sgr_ref = C8H13O_SGRS[0]
    natms = len(graph.atoms(sgr_ref))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        sgr = graph.relabel(sgr_ref, pmt_dct)
        assert graph.from_dictionaries(
            graph.atom_symbols(sgr),
            graph.bond_keys(sgr),
            atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(sgr),
            atm_ste_par_dct=graph.atom_stereo_parities(sgr),
            bnd_ste_par_dct=graph.bond_stereo_parities(sgr)
        ) == sgr


# # transformations
def test__graph__without_bond_orders():
    """ test graph.without_bond_orders
    """
    assert graph.without_stereo_parities(C8H13O_RGRS[0]) == C8H13O_CGR


def test__graph__without_stereo_parities():
    """ test graph.without_stereo_parities
    """
    assert graph.without_stereo_parities(C8H13O_SGRS[0]) == C8H13O_CGR


# graph theory library
# # atom properties
def test__atom_neighbor_keys():
    """ test graph.atom_neighbor_keys
    """
    assert graph.atom_neighbor_keys(C8H13O_CGR) == {
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


def test__atom_bond_keys():
    """ test graph.atom_neighbor_keys
    """
    assert graph.atom_bond_keys(C8H13O_CGR) == {
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
def test__bond_neighbor_keys():
    """ test graph.bond_neighbor_keys
    """
    assert graph.bond_neighbor_keys(C8H13O_CGR) == {
        frozenset({1, 4}): frozenset({frozenset({4, 6})}),
        frozenset({4, 6}): frozenset({frozenset({6, 7}), frozenset({1, 4}),
                                      frozenset({2, 6})}),
        frozenset({2, 6}): frozenset({frozenset({6, 7}), frozenset({4, 6})}),
        frozenset({0, 3}): frozenset({frozenset({3, 5})}),
        frozenset({6, 7}): frozenset({frozenset({4, 6}), frozenset({8, 7}),
                                      frozenset({5, 7}), frozenset({2, 6})}),
        frozenset({8, 7}): frozenset({frozenset({6, 7}), frozenset({5, 7})}),
        frozenset({3, 5}): frozenset({frozenset({5, 7}), frozenset({0, 3})}),
        frozenset({5, 7}): frozenset({frozenset({6, 7}), frozenset({3, 5}),
                                      frozenset({8, 7})})
    }


# # other properties
def test__branch():
    """ test graph.branch
    """
    assert graph.branch(C8H13O_CGR, 6, frozenset({6, 4})) == (
        {1: ('C', 2, None), 4: ('C', 1, None), 6: ('C', 1, None)},
        {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None)}
    )


def test__rings():
    """ test graph.rings
    """
    assert graph.rings(C5H5N5O_CGR) == (
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


def test__delete_atoms():
    """ test graph.delete_atoms
    """
    assert graph.delete_atoms(C3H3_CGR, (0,)) == (
        {1: ('C', 1, None), 2: ('C', 1, None)},
        {frozenset({1, 2}): (1, None)})


# connectivity graph library
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


# # transformations
def test__add_explicit_hydrogens():
    """ test graph.add_explicit_hydrogens
    """
    assert graph.add_explicit_hydrogens(C_CGR, {0: 1}) == (
        {0: ('C', 0, None), 1: ('H', 0, None)},
        {frozenset({0, 1}): (1, None)})


def test__implicit():
    """ test graph.implicit
    """
    assert graph.implicit(CH2FH2H_CGR_EXP) == CH2FH2H_CGR


def test__explicit():
    """ test graph.explicit
    """
    assert graph.explicit(CH2FH2H_CGR) == CH2FH2H_CGR_EXP


# # comparisons
def test__backbone_isomorphic():
    """ test graph.backbone_isomorphic
    """
    assert graph.backbone_isomorphic(CH2FH2H_CGR, CH2FH2H_CGR_EXP)

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
    assert graph.backbone_unique(C2H2CL2F2_SGRS) == (
        C2H2CL2F2_SGRS[0], C2H2CL2F2_SGRS[1], C2H2CL2F2_SGRS[3]
    )


# inchi conversion library
def test__atom_inchi_numbers():
    """ test graph.atom_inchi_numbers
    """
    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        inv_pmt_dct = dict(map(reversed, pmt_dct.items()))
        assert graph.atom_inchi_numbers(cgr_pmt) == inv_pmt_dct


def test__inchi():
    """ test graph.inchi
    """
    assert graph.inchi(C_CGR) == C_ICH
    assert graph.inchi(C2_CGR) == C2_ICH
    assert graph.inchi(C3H3_CGR) == C3H3_ICH
    assert graph.inchi(CH2FH2H_CGR) == CH2FH2H_ICH
    assert graph.inchi(CH2FH2H_CGR_EXP) == CH2FH2H_ICH
    assert graph.inchi(C5H5N5O_CGR) == C5H5N5O_ICH
    assert graph.inchi(C8H13O_CGR) == C8H13O_ICH
    assert graph.inchi(C3H3CL2F3_CGR) == C3H3CL2F3_ICH
    assert graph.inchi(C3H5N3_CGR) == C3H5N3_ICH
    assert graph.inchi(C2H2CL2F2_CGR) == C2H2CL2F2_ICH


# resonance library
# # atom properties
def test__atom_bond_valences():
    """ test graph.atom_bond_valences
    """
    assert graph.atom_bond_valences(C8H13O_CGR) == {
        0: 4, 1: 3, 2: 4, 3: 3, 4: 3, 5: 3, 6: 4, 7: 4, 8: 1}


def test__atom_radical_valences():
    """ test graph.atom_radical_valences
    """
    assert graph.atom_radical_valences(C8H13O_CGR) == {
        0: 0, 1: 1, 2: 0, 3: 1, 4: 1, 5: 1, 6: 0, 7: 0, 8: 1}


# # bond properties
def test__resonance_dominant_bond_orders():
    """ test grpah.resonance_dominant_bond_orders
    """
    assert graph.resonance_dominant_bond_orders(C3H3_CGR) == {
        frozenset({0, 1}): frozenset({1, 2}),
        frozenset({0, 2}): frozenset({1, 2}),
        frozenset({1, 2}): frozenset({1, 2})
    }


# # other properties
def test__maximum_spin_multiplicity():
    """ test graph.maximum_spin_multiplicity
    """
    assert graph.maximum_spin_multiplicity(C_CGR) == 5
    assert graph.maximum_spin_multiplicity(C2_CGR) == 7


def test__possible_spin_multiplicities():
    """ test graph.possible_spin_multiplicities
    """
    assert graph.possible_spin_multiplicities(C_CGR) == (1, 3, 5)
    assert graph.possible_spin_multiplicities(C2_CGR) == (1, 3, 5, 7)


# # transformations
def test__resonances():
    """ test graph.resonances
    """
    assert graph.resonances(C2_CGR) == C2_RGRS
    assert graph.resonances(C3H3_CGR) == C3H3_RGRS
    assert graph.resonances(C8H13O_CGR) == C8H13O_RGRS


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


# stereo library
# # properties
def test__stereo_inchi():
    """ test graph.stereo_inchi
    """
    assert tuple(map(graph.stereo_inchi, C2H2CL2F2_SGRS)) == C2H2CL2F2_STE_ICHS
    assert tuple(map(graph.stereo_inchi, C8H13O_SGRS)) == C8H13O_STE_ICHS


def test__is_chiral():
    """ test graph.is_chiral
    """
    assert graph.is_chiral(C2H2CL2F2_SGRS[0]) is True
    assert graph.is_chiral(C2H2CL2F2_SGRS[1]) is False
    assert graph.is_chiral(C3H3CL2F3_SGRS[0]) is True
    assert graph.is_chiral(C3H3CL2F3_SGRS[2]) is False


def test__stereogenic_atom_keys():
    """ test graph.stereogenic_atom_keys
    """
    assert graph.stereogenic_atom_keys(C8H13O_CGR) == frozenset({6, 7})
    assert graph.stereogenic_atom_keys(C3H3CL2F3_CGR) == frozenset({1, 2})


def test__stereogenic_bond_keys():
    """ test graph.stereogenic_bond_keys
    """
    assert graph.stereogenic_bond_keys(C8H13O_CGR) == frozenset(
        {frozenset({3, 5})})
    assert graph.stereogenic_bond_keys(C3H5N3_CGR) == frozenset(
        {frozenset({1, 4}), frozenset({0, 3})})


# # transformations
def test__stereomers():
    """ test graph.stereomers
    """
    assert graph.stereomers(C2H2CL2F2_CGR) == C2H2CL2F2_SGRS
    assert graph.stereomers(C3H3CL2F3_CGR) == C3H3CL2F3_SGRS
    assert graph.stereomers(C3H5N3_CGR) == C3H5N3_SGRS
    assert graph.stereomers(C8H13O_CGR) == C8H13O_SGRS


def test__substereomers():
    """ test graph.substereomers
    """
    partial_sgr = graph.set_atom_stereo_parities(C8H13O_CGR, {6: True})
    assert graph.substereomers(partial_sgr) == C8H13O_SGRS[4:]

    partial_sgr = graph.set_bond_stereo_parities(
        C8H13O_CGR, {frozenset({3, 5}): False})
    assert graph.substereomers(partial_sgr) == C8H13O_SGRS[0::2]


# # comparisons
def test__enantiomerically_unique():
    """ test graph.enantiomerically_unique
    """
    assert graph.enantiomerically_unique(C3H3CL2F3_SGRS) == (
        C3H3CL2F3_SGRS[0], C3H3CL2F3_SGRS[2], C3H3CL2F3_SGRS[4]
    )


# misc
def test__bond_symmetry_numbers():
    """ test graph.bond_symmetry_numbers
    """
    assert graph.bond_symmetry_numbers(C8H13O_CGR) == {
        frozenset({1, 4}): 1, frozenset({4, 6}): 1, frozenset({2, 6}): 3,
        frozenset({0, 3}): 3, frozenset({6, 7}): 1, frozenset({8, 7}): 1,
        frozenset({3, 5}): 1, frozenset({5, 7}): 1}


if __name__ == '__main__':
    # test__from_atoms_and_bonds()
    # test__from_dictionaries()
    # test__atom_neighbor_keys()
    # test__atom_bond_keys()
    # test__bond_neighbor_keys()
    # test__branch()
    # test__rings()
    # test__bond_induced_subgraph()
    # test__relabel()
    # test__delete_atoms()
    # test__atom_explicit_hydrogen_valences()
    # test__atom_explicit_hydrogen_keys()
    # test__backbone_keys()
    # test__explicit_hydrogen_keys()
    # test__add_explicit_hydrogens()
    # test__implicit()
    # test__inchi()
    # test__atom_bond_valences()
    # test__atom_radical_valences()
    # test__maximum_spin_multiplicity()
    # test__possible_spin_multiplicities()
    # test__resonances()
    # test__subresonances()
    # test__dominant_resonances()
    # test__stereo_inchi()
    # test__stereomers()
    # test__is_chiral()
    # test__stereogenic_atom_keys()
    # test__stereogenic_bond_keys()
    # test__stereomers()
    # test__substereomers()
    test__stereomers()
    test__backbone_unique()
    test__enantiomerically_unique()
    test__stereo_inchi()
    test__bond_symmetry_numbers()
