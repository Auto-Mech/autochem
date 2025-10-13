"""test graph.ts"""

import automol
import numpy
import pytest
from automol import graph

# Sn2 Atom Stereo
# FC(Cl)N + [OH] => FC(O)N + [Cl]
#  *
# [* marks an Sn2 inversion site]
CH4CLFNO_TSG = (
    {
        0: ("C", 0, True),
        1: ("Cl", 0, None),
        2: ("F", 0, None),
        3: ("N", 0, None),
        4: ("H", 0, None),
        5: ("H", 0, None),
        6: ("H", 0, None),
        7: ("O", 0, None),
        8: ("H", 0, None),
    },
    {
        frozenset({0, 3}): (1, None),
        frozenset({0, 1}): (0.9, None),
        frozenset({0, 2}): (1, None),
        frozenset({3, 6}): (1, None),
        frozenset({0, 4}): (1, None),
        frozenset({7, 8}): (1, None),
        frozenset({3, 5}): (1, None),
        frozenset({0, 7}): (0.1, None),
    },
)
CH4CLFNO_GEO = (
    ("C", (-1.001524, -0.381178, 0.270439)),
    ("Cl", (-1.672496, 2.593418, -1.321151)),
    ("F", (-1.819841, -1.982093, -1.511983)),
    ("N", (1.692691, -0.383608, 0.480908)),
    ("H", (-2.016387, -0.255401, 2.022491)),
    ("H", (2.404653, -0.577165, -1.299886)),
    ("H", (2.266771, 1.372256, 1.055224)),
    ("O", (-0.80741, -3.604572, 2.106322)),
    ("H", (-1.123524, -4.811378, 0.495377)),
)

# Fleeting Atom Stereo
# CCOCC + [OH] => C[CH]OCC + O
#  *
# [* marks a fleeting TS stereosite]
C4H11O2_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, False),
        3: ("C", 0, None),
        4: ("O", 0, None),
        5: ("H", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
        13: ("H", 0, None),
        14: ("H", 0, None),
        15: ("O", 0, None),
        16: ("H", 0, None),
    },
    {
        frozenset({0, 5}): (1, None),
        frozenset({3, 14}): (1, None),
        frozenset({11, 15}): (0.1, None),
        frozenset({2, 11}): (0.9, None),
        frozenset({1, 10}): (1, None),
        frozenset({3, 4}): (1, None),
        frozenset({1, 9}): (1, None),
        frozenset({0, 6}): (1, None),
        frozenset({0, 2}): (1, None),
        frozenset({3, 13}): (1, None),
        frozenset({2, 12}): (1, None),
        frozenset({15, 16}): (1, None),
        frozenset({2, 4}): (1, None),
        frozenset({1, 8}): (1, None),
        frozenset({0, 7}): (1, None),
        frozenset({1, 3}): (1, None),
    },
)
C4H11O2_GEO = (
    ("C", (4.33819, 1.441087, 0.295484)),
    ("C", (-4.574912, 0.214792, -0.174836)),
    ("C", (2.316425, -0.587509, 0.260823)),
    ("C", (-2.077647, -1.187579, -0.205727)),
    ("O", (-0.0881, 0.60958, 0.056152)),
    ("H", (6.220916, 0.603391, 0.451236)),
    ("H", (4.047672, 2.732272, 1.886248)),
    ("H", (4.255128, 2.582155, -1.4288)),
    ("H", (-6.159774, -1.096745, -0.371638)),
    ("H", (-4.65463, 1.594987, -1.714572)),
    ("H", (-4.796227, 1.270025, 1.591399)),
    ("H", (2.682893, -1.715305, -1.436134)),
    ("H", (2.341862, -1.826667, 1.919052)),
    ("H", (-2.01011, -2.538954, 1.361152)),
    ("H", (-1.866407, -2.209146, -1.993927)),
    ("O", (3.7089, -3.401617, -3.726459)),
    ("H", (3.282065, -5.016675, -2.48937)),
)

# Fleeting Bond Stereo
# C=C(O[O])OO => [CH]=C(OO)OO
#  *
# [* marks a fleeting TS stereosite]
C2H3O4_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("O", 0, None),
        3: ("O", 0, None),
        4: ("O", 0, None),
        5: ("O", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
    },
    {
        frozenset({2, 8}): (1, None),
        frozenset({1, 4}): (1, None),
        frozenset({0, 6}): (0.9, None),
        frozenset({0, 1}): (1, None),
        frozenset({3, 6}): (0.1, None),
        frozenset({2, 4}): (1, None),
        frozenset({1, 5}): (1, None),
        frozenset({3, 5}): (1, None),
        frozenset({0, 7}): (1, None),
    },
)
C2H3O4_GEO = (
    ("C", (-2.812896, 0.535324, 0.864225)),
    ("C", (-0.773943, -0.715224, 0.054102)),
    ("O", (3.585469, -2.437992, -0.732028)),
    ("O", (0.896435, 2.730653, -1.397099)),
    ("O", (1.399056, -2.02279, 0.305438)),
    ("O", (-0.803513, 0.372439, -2.38697)),
    ("H", (-2.22919, 0.44838, 2.621235)),
    ("H", (-4.102653, 2.150414, 1.056667)),
    ("H", (4.841235, -1.061205, -0.38557)),
)

# Fleeting Atom Stereo + Conserved Atom Stereo
# CCO[C@H](O[O])C => C[CH]O[C@H](OO)C
#  *
# [* marks a fleeting stereosite]
C4H9O3_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, False),
        3: ("C", 0, True),
        4: ("O", 0, None),
        5: ("O", 0, None),
        6: ("O", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
        13: ("H", 0, None),
        14: ("H", 0, None),
        15: ("H", 0, None),
    },
    {
        frozenset({4, 6}): (1, None),
        frozenset({3, 6}): (1, None),
        frozenset({1, 12}): (1, None),
        frozenset({2, 14}): (1, None),
        frozenset({1, 10}): (1, None),
        frozenset({0, 8}): (1, None),
        frozenset({0, 9}): (1, None),
        frozenset({2, 13}): (0.9, None),
        frozenset({3, 15}): (1, None),
        frozenset({4, 13}): (0.1, None),
        frozenset({1, 11}): (1, None),
        frozenset({0, 2}): (1, None),
        frozenset({2, 5}): (1, None),
        frozenset({3, 5}): (1, None),
        frozenset({0, 7}): (1, None),
        frozenset({1, 3}): (1, None),
    },
)
C4H9O3_GEO = (
    ("C", (-3.858804, -0.017433, -0.816915)),
    ("C", (3.284911, 2.431032, 0.462603)),
    ("C", (-1.959761, -1.667044, 0.44462)),
    ("C", (2.02811, -0.049194, 0.203738)),
    ("O", (1.521606, -4.246381, 1.551043)),
    ("O", (-0.22096, 0.069541, 1.635187)),
    ("O", (3.101169, -2.349906, 0.925334)),
    ("H", (-3.764188, 2.014997, -0.161893)),
    ("H", (-5.805311, -0.707311, -0.433068)),
    ("H", (-3.483251, 0.040265, -2.87643)),
    ("H", (3.006888, 3.239081, 2.415149)),
    ("H", (2.436056, 3.760618, -0.93411)),
    ("H", (5.336483, 2.305681, 0.107293)),
    ("H", (-1.992606, -3.145298, 1.82341)),
    ("H", (-1.040456, -2.832608, -1.161927)),
    ("H", (1.546318, -0.226049, -1.899431)),
)

# (Reactant Bond Stereo => Product Atom Stereo)
# F/C=C/F + [OH] => F[CH][C@H](O)F
C2H3F2O_TSG = (
    {
        0: ("C", 0, False),
        1: ("C", 0, None),
        2: ("F", 0, None),
        3: ("F", 0, None),
        4: ("H", 0, None),
        5: ("H", 0, None),
        6: ("O", 0, None),
        7: ("H", 0, None),
    },
    {
        frozenset({6, 7}): (1, None),
        frozenset({1, 5}): (1, None),
        frozenset({0, 6}): (0.1, None),
        frozenset({0, 1}): (1, True),
        frozenset({0, 2}): (1, None),
        frozenset({1, 3}): (1, None),
        frozenset({0, 4}): (1, None),
    },
)
C2H3F2O_GEO = (
    ("C", (1.297288, 0.39764, 0.454963)),
    ("C", (-1.326186, -0.366377, 0.15782)),
    ("F", (3.477276, -1.308171, 0.493162)),
    ("F", (-3.438916, 1.322689, -0.431268)),
    ("H", (2.079235, 2.259128, 0.133371)),
    ("H", (-2.02146, -2.249426, -0.230421)),
    ("O", (1.234258, 0.222229, -4.125039)),
    ("H", (0.617002, 2.080766, -3.428028)),
)

# Conserved Atom Stereo + (Reactant Bond Stereo => Product Atom Stereo)
# F/C=C([C@@H](F)O)\[C@H](F)O + [OH] => F[C@H]([C]([C@@H](F)O)[C@H](F)O)O
C4H5F3O2_TSG = (
    {
        0: ("C", 0, False),
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
        frozenset({7, 12}): (1, None),
        frozenset({2, 10}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({0, 1}): (1, False),
        frozenset({3, 6}): (1, None),
        frozenset({2, 5}): (1, None),
        frozenset({0, 4}): (1, None),
        frozenset({3, 8}): (1, None),
        frozenset({0, 14}): (0.1, None),
        frozenset({1, 3}): (1, None),
        frozenset({8, 13}): (1, None),
        frozenset({14, 15}): (1, None),
        frozenset({3, 11}): (1, None),
        frozenset({2, 7}): (1, None),
        frozenset({0, 9}): (1, None),
    },
)
C4H5F3O2_GEO = (
    ("C", (-1.141471, 3.060908, -1.372282)),
    ("C", (-0.335935, 0.809777, -0.552249)),
    ("C", (-2.015423, -0.956788, 0.900211)),
    ("C", (2.245586, -0.127934, -1.109997)),
    ("F", (0.402603, 4.723624, -2.508427)),
    ("F", (-2.955248, -2.766816, -0.553588)),
    ("F", (2.089015, -2.665134, -1.889925)),
    ("O", (-0.608473, -1.927482, 3.005918)),
    ("O", (3.834623, -0.044323, 0.926263)),
    ("H", (-3.009159, 3.830615, -1.183382)),
    ("H", (-3.595336, 0.12939, 1.788313)),
    ("H", (3.088668, 0.820776, -2.776933)),
    ("H", (-1.019998, -3.680692, 3.265653)),
    ("H", (3.061224, -1.14277, 2.195312)),
    ("O", (-1.286205, 1.482841, -4.596348)),
    ("H", (-3.168531, 0.849827, -3.99654)),
)

# Non-strict Fleeting Atom Stereo + (Reactant Atom Stereo => Product Bond Stereo)
# CC[C@H](O[O])C => C/C=C/C + O[O]
#  *
# [* a non-strict fleeting stereocenter, which is implied by the reactant and product
# configurations]
C4H9O2_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, True),
        3: ("C", 0, True),
        4: ("O", 0, None),
        5: ("O", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
        13: ("H", 0, None),
        14: ("H", 0, None),
    },
    {
        frozenset({2, 13}): (1, None),
        frozenset({1, 9}): (1, None),
        frozenset({0, 6}): (1, None),
        frozenset({2, 3}): (1, None),
        frozenset({1, 11}): (1, None),
        frozenset({4, 5}): (1, None),
        frozenset({0, 2}): (1, None),
        frozenset({4, 12}): (0.1, None),
        frozenset({2, 12}): (0.9, None),
        frozenset({3, 14}): (1, None),
        frozenset({3, 5}): (0.9, None),
        frozenset({1, 3}): (1, None),
        frozenset({0, 7}): (1, None),
        frozenset({1, 10}): (1, None),
        frozenset({0, 8}): (1, None),
    },
)
C4H9O2_GEO = (
    ("C", (2.840842, -0.918093, 0.301416)),
    ("C", (-2.107066, -1.899196, -0.110792)),
    ("C", (0.987221, 1.052291, 1.070762)),
    ("C", (-1.279309, 0.781779, -0.668652)),
    ("O", (-1.659829, 4.677101, -0.256324)),
    ("O", (-2.952545, 2.503599, 0.312223)),
    ("H", (2.360173, -1.778041, -1.574524)),
    ("H", (4.773451, -0.246703, 0.187287)),
    ("H", (2.703926, -2.528999, 1.689233)),
    ("H", (-3.298375, -2.425681, -1.754879)),
    ("H", (-0.350923, -3.043722, 0.066643)),
    ("H", (-3.108548, -1.931614, 1.691764)),
    ("H", (1.635613, 2.957619, 0.475404)),
    ("H", (0.428192, 1.019683, 3.075372)),
    ("H", (-0.462632, 1.023199, -2.601261)),
)


# Strict Fleeting Atom Stereo + Product Bond Stereo
# CCCC(O[O])CCC => CC/C=C/CCC + [O]O
#    *
# [* similar to the case above, but not present in reactant due to symmetry]
C7H15O2_TSG = (
    {
        0: ("C", 3, None),
        1: ("C", 3, None),
        2: ("C", 1, False),
        3: ("C", 1, False),
        4: ("C", 2, None),
        5: ("C", 2, None),
        6: ("C", 2, None),
        7: ("O", 0, None),
        8: ("O", 0, None),
        9: ("H", 0, None),
    },
    {
        frozenset({2, 3}): (1, None),
        frozenset({3, 6}): (1, None),
        frozenset({0, 4}): (1, None),
        frozenset({7, 9}): (0.1, None),
        frozenset({7, 8}): (1, None),
        frozenset({2, 4}): (1, None),
        frozenset({3, 8}): (0.9, None),
        frozenset({5, 6}): (1, None),
        frozenset({1, 5}): (1, None),
        frozenset({2, 9}): (0.9, None),
    },
)


# Reactant Bond Stereo => Product Bond Stereo (flipped)
# [H]\[C]=C/O + [OH] => O/C=C/O
#        *                 *
# [* parity flips due to vinyl adddition]
C2H4O2_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("O", 0, None),
        3: ("H", 0, None),
        4: ("H", 0, None),
        5: ("H", 0, None),
        6: ("O", 0, None),
        7: ("H", 0, None),
    },
    {
        frozenset({1, 4}): (1, None),
        frozenset({0, 6}): (0.1, None),
        frozenset({1, 2}): (1, None),
        frozenset({0, 3}): (1, None),
        frozenset({0, 1}): (1, True),
        frozenset({6, 7}): (1, None),
        frozenset({2, 5}): (1, None),
    },
)
C2H4O2_GEO = (
    ("C", (-1.759994, 0.859003, 0.020027)),
    ("C", (-0.10802, -0.946706, -0.598034)),
    ("O", (2.406749, -0.523634, -0.917962)),
    ("H", (-0.841145, 2.68758, 0.235109)),
    ("H", (-0.675523, -2.920932, -0.895759)),
    ("H", (2.952195, 1.207736, -0.664496)),
    ("O", (-4.196811, -1.722699, -0.095209)),
    ("H", (-3.621029, -2.105771, 1.867762)),
)

# 2x Reactant Bond Stereo => Product Bond Stereo (flipped)
# F/C=[C]/[H] + F/C=[C]/[H] => F/C=C\C=C/F
#    *             *              *   *
# [* parity flips due to vinyl adddition]
C4H4F2_TSG = (
    {
        0: ("F", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, None),
        3: ("H", 0, None),
        4: ("H", 0, None),
        5: ("H", 0, None),
        6: ("C", 0, None),
        7: ("C", 0, None),
        8: ("F", 0, None),
        9: ("H", 0, None),
    },
    {
        frozenset({1, 4}): (1, None),
        frozenset({2, 3}): (1, None),
        frozenset({1, 2}): (1, False),
        frozenset({2, 6}): (0.1, None),
        frozenset({0, 1}): (1, None),
        frozenset({6, 7}): (1, False),
        frozenset({7, 9}): (1, None),
        frozenset({7, 8}): (1, None),
        frozenset({5, 6}): (1, None),
    },
)
C4H4F2_GEO = (
    ("F", (-2.632951, -1.443686, -0.650948)),
    ("C", (-0.79737, 0.282099, -0.296882)),
    ("C", (1.578817, -0.181158, -0.27474)),
    ("H", (2.842947, 1.380913, -0.127814)),
    ("H", (-1.712989, 2.122868, -0.319623)),
    ("H", (4.801748, -4.030195, 0.216704)),
    ("C", (2.82878, -3.660381, 0.0426)),
    ("C", (1.291958, -5.531084, 0.046779)),
    ("F", (-1.246023, -5.354531, 0.023789)),
    ("H", (1.758437, -7.533401, 0.067458)),
)

# Constrained TS Bond Stereo
# [CH2]/C=C\CC=C => C1C[CH]CCC=1
#        *
# [* marks a constrained bond stereocenter -- must be cis]
C6H9_TSG = (
    {
        0: ("C", 2, None),
        1: ("C", 2, None),
        2: ("C", 1, None),
        3: ("C", 1, None),
        4: ("C", 1, None),
        5: ("C", 2, None),
    },
    {
        frozenset({4, 5}): (1, None),
        frozenset({0, 1}): (0.1, None),
        frozenset({0, 2}): (1, None),
        frozenset({2, 4}): (1, None),
        frozenset({3, 5}): (1, None),
        frozenset({1, 3}): (1, None),
    },
)

# Constrained TS Bond Stereo
# [CH2]/C=C\[C@@H](CC)OO => CC[C@H]1OCC=C1 + [OH]
#        *
# [* marks a constrained bond stereocenter -- must be cis]
C6H11O2_TSG = (
    {
        0: ("C", 0, None),
        1: ("H", 0, None),
        2: ("H", 0, None),
        3: ("C", 0, None),
        4: ("C", 0, None),
        5: ("C", 0, False),
        6: ("H", 0, None),
        7: ("C", 0, None),
        8: ("C", 0, None),
        9: ("O", 0, None),
        10: ("O", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
        13: ("H", 0, None),
        14: ("H", 0, None),
        15: ("H", 0, None),
        16: ("H", 0, None),
        17: ("H", 0, None),
        18: ("H", 0, None),
    },
    {
        frozenset({4, 5}): (1, None),
        frozenset({0, 1}): (1, None),
        frozenset({7, 8}): (1, None),
        frozenset({7, 14}): (1, None),
        frozenset({10, 18}): (1, None),
        frozenset({8, 16}): (1, None),
        frozenset({5, 9}): (1, None),
        frozenset({5, 7}): (1, None),
        frozenset({0, 9}): (0.1, None),
        frozenset({3, 4}): (1, None),
        frozenset({8, 17}): (1, None),
        frozenset({9, 10}): (0.9, None),
        frozenset({7, 13}): (1, None),
        frozenset({0, 3}): (1, None),
        frozenset({0, 2}): (1, None),
        frozenset({4, 12}): (1, None),
        frozenset({5, 6}): (1, None),
        frozenset({3, 11}): (1, None),
        frozenset({8, 15}): (1, None),
    },
)
# C6H11O2_GEO = (
#     ('C', (4.922102, -0.449109, -0.140716)),
#     ('H', (6.745736, 0.606405, 0.265351)),
#     ('H', (4.118399, -0.660212, 1.822398)),
#     ('C', (3.212445, 0.892843, -1.952966)),
#     ('C', (0.697707, 0.795413, -1.846221)),
#     ('C', (-0.747979, -0.689044, 0.330813)),
#     ('H', (0.20965, 0.275312, 2.099735)),
#     ('C', (-3.326397, 0.151782, 0.65621)),
#     ('C', (-3.153725, 2.957264, 1.641096)),
#     ('O', (1.273343, -2.63089, 0.385426)),
#     ('O', (0.295383, -4.928983, -0.643783)),
#     ('H', (4.311136, 1.35876, -3.540684)),
#     ('H', (-0.750563, 1.001214, -3.184155)),
#     ('H', (-4.26182, 0.460045, -1.176939)),
#     ('H', (-4.495947, -0.850012, 1.98946)),
#     ('H', (-1.97707, 2.856142, 3.362061)),
#     ('H', (-4.938737, 3.791284, 2.051566)),
#     ('H', (-2.076922, 4.046719, 0.258913)),
#     ('H', (-0.398204, -5.654819, 0.896648)))

# Constrained TS Vinyl Bond Stereo + Atom Stereo
# C=C[C@@H](CCC)O[O] => [H]/[C]=C/[C@@H](CCC)OO
#                              *
# [* marks a constrained bond stereocenter -- must be cis]
C6H11O2b_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, False),
        3: ("C", 0, None),
        4: ("C", 0, None),
        5: ("C", 0, None),
        6: ("O", 0, None),
        7: ("O", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
        13: ("H", 0, None),
        14: ("H", 0, None),
        15: ("H", 0, None),
        16: ("H", 0, None),
        17: ("H", 0, None),
        18: ("H", 0, None),
    },
    {
        frozenset({2, 3}): (1, None),
        frozenset({2, 6}): (1, None),
        frozenset({4, 5}): (1, None),
        frozenset({5, 17}): (1, None),
        frozenset({0, 1}): (1, None),
        frozenset({6, 7}): (1, None),
        frozenset({7, 8}): (0.1, None),
        frozenset({5, 16}): (1, None),
        frozenset({5, 18}): (1, None),
        frozenset({4, 14}): (1, None),
        frozenset({2, 11}): (1, None),
        frozenset({1, 10}): (1, None),
        frozenset({0, 8}): (0.9, None),
        frozenset({0, 9}): (1, None),
        frozenset({3, 4}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({4, 15}): (1, None),
        frozenset({3, 13}): (1, None),
        frozenset({3, 12}): (1, None),
    },
)
# C6H11O2b_GEO = (
#     ('C', (5.787848, -0.151178, -1.541546)),
#     ('C', (3.461157, -0.812373, -1.686212)),
#     ('C', (1.549103, -0.250423, 0.355867)),
#     ('C', (-1.120547, -0.550261, -0.473615)),
#     ('C', (-2.692891, 0.085392, 1.874711)),
#     ('C', (-5.418684, 0.373202, 1.179352)),
#     ('O', (2.431312, -1.320658, 2.595525)),
#     ('O', (4.750424, -0.318698, 2.953194)),
#     ('H', (6.447453, 0.650469, 0.710288)),
#     ('H', (7.414901, -0.679888, -2.576486)),
#     ('H', (2.744884, -1.92163, -3.292046)),
#     ('H', (1.822851, 1.835276, 0.730881)),
#     ('H', (-1.514309, 0.862512, -1.988341)),
#     ('H', (-1.622245, -2.418038, -1.184909)),
#     ('H', (-2.370337, -1.473728, 3.210659)),
#     ('H', (-2.037551, 1.880908, 2.67804)),
#     ('H', (-5.842784, 2.354336, 0.684703)),
#     ('H', (-5.89428, -0.848404, -0.496429)),
#     ('H', (-6.648168, -0.229583, 2.735668)))

# Constrained TS Vinyl Bond Stereo + Atom Stereo
# C=C[C]=COC => C=CC=CO[CH2]
#       *           *
# [* marks a constrained bond stereocenter -- must be cis]
C5H7O_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, None),
        3: ("C", 0, None),
        4: ("O", 0, None),
        5: ("C", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
    },
    {
        frozenset({3, 4}): (1, None),
        frozenset({5, 11}): (1, None),
        frozenset({0, 6}): (1, None),
        frozenset({2, 3}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({3, 9}): (1, None),
        frozenset({4, 5}): (1, None),
        frozenset({0, 1}): (1, None),
        frozenset({5, 12}): (1, None),
        frozenset({2, 10}): (0.1, None),
        frozenset({1, 8}): (1, None),
        frozenset({5, 10}): (0.9, None),
        frozenset({0, 7}): (1, None),
    },
)

# Constrained TS Bond Stereo + (Reactant Bond Stereo => Product Atom Stereo)
# C/C=C/C + C=CC=C => C[C@H]1[C@H](C)CC=CC1
#    ^                  ^     ^        *
# [* marks a constrained bond stereocenter -- must be cis]
# [^ marks bond stereo which becomes atom stereo]
# This is a Diels-Alder reaction, which we don't yet have a classifier for, but we
# should still be able to handle the other routines
C8H14_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, False),
        2: ("C", 0, None),
        3: ("C", 0, None),
        4: ("C", 0, None),
        5: ("C", 0, None),
        6: ("C", 0, False),
        7: ("C", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
        13: ("H", 0, None),
        14: ("H", 0, None),
        15: ("H", 0, None),
        16: ("H", 0, None),
        17: ("H", 0, None),
        18: ("H", 0, None),
        19: ("H", 0, None),
        20: ("H", 0, None),
        21: ("H", 0, None),
    },
    {
        frozenset({0, 1}): (1, None),
        frozenset({1, 2}): (0.1, None),
        frozenset({2, 3}): (1, None),
        frozenset({3, 4}): (1, None),
        frozenset({4, 5}): (1, None),
        frozenset({5, 6}): (0.1, None),
        frozenset({1, 6}): (1, None),
        frozenset({6, 7}): (1, None),
        frozenset({0, 8}): (1, None),
        frozenset({0, 9}): (1, None),
        frozenset({0, 10}): (1, None),
        frozenset({1, 11}): (1, None),
        frozenset({2, 12}): (1, None),
        frozenset({2, 13}): (1, None),
        frozenset({3, 14}): (1, None),
        frozenset({4, 15}): (1, None),
        frozenset({5, 16}): (1, None),
        frozenset({5, 17}): (1, None),
        frozenset({6, 18}): (1, None),
        frozenset({7, 19}): (1, None),
        frozenset({7, 20}): (1, None),
        frozenset({7, 21}): (1, None),
    },
)
C8H14_GEO = (
    ("C", (3.232417, -0.26753, -0.403764)),
    ("C", (0.661159, 0.068665, -1.68409)),
    ("C", (1.593441, 3.648706, -2.746415)),
    ("C", (-0.325404, 4.334356, -4.220785)),
    ("C", (-2.922068, 4.459019, -3.245847)),
    ("C", (-4.124938, 2.264194, -2.888306)),
    ("C", (-1.64472, 0.410532, -0.542585)),
    ("C", (-3.610366, -1.684613, -0.12205)),
    ("H", (3.819298, -2.232313, -0.54354)),
    ("H", (4.589412, 0.922721, -1.404988)),
    ("H", (3.162991, 0.299921, 1.570763)),
    ("H", (0.599165, -0.990403, -3.491541)),
    ("H", (1.620005, 3.988812, -0.692902)),
    ("H", (3.526105, 3.442252, -3.486798)),
    ("H", (-0.24523, 4.28609, -6.280982)),
    ("H", (-3.468856, 6.230918, -2.384826)),
    ("H", (-4.499243, 1.058283, -4.548879)),
    ("H", (-5.552061, 2.017693, -1.411803)),
    ("H", (-1.600619, 1.673979, 1.134201)),
    ("H", (-5.367419, -0.788789, 0.472015)),
    ("H", (-2.926932, -2.933633, 1.355101)),
    ("H", (-3.926188, -2.734536, -1.858791)),
)

# Unusual linearity case (biradical)
#
# [CH2]/C=[C]/C#C => C=C=C=C=[CH] + [H]
#          -  - -      - - -  -
# [- marks a potentially linear atom]
C5H4_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("H", 0, None),
        3: ("H", 0, None),
        4: ("C", 0, None),
        5: ("C", 0, None),
        6: ("C", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
    },
    {
        frozenset({1, 4}): (1, True),
        frozenset({0, 3}): (1, None),
        frozenset({4, 5}): (1, None),
        frozenset({6, 7}): (1, None),
        frozenset({0, 1}): (1, None),
        frozenset({0, 2}): (1, None),
        frozenset({5, 6}): (1, None),
        frozenset({1, 8}): (0.9, None),
    },
)

# Tricky linearity case (triple bond insertion)
#
# CC=[C]O => C#CO + [CH3]
#            - -
# [- marks a potentially linear atom]
C3H5O_TSG1 = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("O", 0, None),
        3: ("H", 0, None),
        4: ("H", 0, None),
        5: ("C", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
    },
    {
        frozenset({1, 2}): (1, None),
        frozenset({0, 3}): (1, None),
        frozenset({0, 1}): (1, False),
        frozenset({0, 5}): (0.9, None),
        frozenset({2, 4}): (1, None),
        frozenset({5, 6}): (1, None),
        frozenset({5, 8}): (1, None),
        frozenset({5, 7}): (1, None),
    },
)

# Tricky linearity case (triple bond insertion)
#
# CC=[C]O => CC#CO + [H]
#             - -
# [- marks a potentially linear atom]
C3H5O_TSG2 = (
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
    },
    {
        frozenset({0, 6}): (1, None),
        frozenset({2, 3}): (1, None),
        frozenset({3, 7}): (1, None),
        frozenset({1, 2}): (1, False),
        frozenset({0, 1}): (1, None),
        frozenset({0, 4}): (1, None),
        frozenset({0, 5}): (1, None),
        frozenset({1, 8}): (0.9, None),
    },
)

# Tricky linearity case (triple bond sigma addition)
#
# C#[C] + [OH] => C#CO
# -  -            - -
# [- marks a potentially linear atom]
C4H6_TSG = (
    {
        0: ("C", 0, None),
        1: ("C", 0, None),
        2: ("C", 0, None),
        3: ("H", 0, None),
        4: ("H", 0, None),
        5: ("H", 0, None),
        6: ("C", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
    },
    {
        frozenset({6, 9}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({2, 6}): (0.1, None),
        frozenset({0, 3}): (1, None),
        frozenset({0, 1}): (1, None),
        frozenset({6, 7}): (1, None),
        frozenset({0, 4}): (1, None),
        frozenset({0, 5}): (1, None),
        frozenset({6, 8}): (1, None),
    },
)

# Bridghead Atom Stereo Pair
# OO[C@@H]1CC=C[CH]1 => C=1[C@H]2C[C@H](C=1)O2 + [OH]
#    *          *           *      *
# [* marks a pair of bridgehead stereo atoms]
C5H7O2_TSG = (
    {
        0: ("O", 0, None),
        1: ("O", 0, None),
        2: ("C", 0, None),
        3: ("H", 0, None),
        4: ("C", 0, None),
        5: ("C", 0, None),
        6: ("C", 0, None),
        7: ("C", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
        13: ("H", 0, None),
    },
    {
        frozenset({4, 10}): (1, None),
        frozenset({2, 3}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({4, 5}): (1, None),
        frozenset({0, 1}): (0.9, None),
        frozenset({6, 7}): (1, None),
        frozenset({5, 12}): (1, None),
        frozenset({7, 8}): (1, None),
        frozenset({2, 4}): (1, None),
        frozenset({5, 6}): (1, None),
        frozenset({4, 11}): (1, None),
        frozenset({6, 13}): (1, None),
        frozenset({1, 5}): (0.1, None),
        frozenset({2, 7}): (1, None),
        frozenset({0, 9}): (1, None),
    },
)

# Double Bridgehead Atom Pair
# [O]OC1CCC2C1O2 => OOC1[CH]CC2C1O2
#     ^ ^ * *         ^  ^   * *
# [* marks the first pair of bridgehead stereo atoms]
# [^ marks the first pair of bridgehead stereo atoms]
C5H7O3_TSG = (
    {
        0: ("O", 0, None),
        1: ("O", 0, None),
        2: ("C", 0, None),
        3: ("C", 0, None),
        4: ("C", 0, None),
        5: ("C", 0, None),
        6: ("C", 0, None),
        7: ("O", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
        13: ("H", 0, None),
        14: ("H", 0, None),
    },
    {
        frozenset({8, 2}): (1, None),
        frozenset({2, 3}): (1, None),
        frozenset({9, 3}): (0.9, None),
        frozenset({2, 6}): (1, None),
        frozenset({4, 5}): (1, None),
        frozenset({0, 1}): (1, None),
        frozenset({6, 7}): (1, None),
        frozenset({6, 14}): (1, None),
        frozenset({5, 13}): (1, None),
        frozenset({5, 7}): (1, None),
        frozenset({0, 9}): (0.1, None),
        frozenset({3, 4}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({4, 12}): (1, None),
        frozenset({10, 3}): (1, None),
        frozenset({5, 6}): (1, None),
        frozenset({11, 4}): (1, None),
    },
)


# Double Bridgehead Pair with Bridgehead Transfer
# C1CC2[CH]C1O2 => [CH]1CC2CC1O2
# ^  *  ^  *        ^    * ^*
# [* marks the first pair of bridgehead stereo atoms]
# [^ marks the first pair of bridgehead stereo atoms]
C5H7O_B_TSG = (
    {
        0: ("C", 0, False),
        1: ("C", 0, False),
        2: ("C", 0, None),
        3: ("C", 0, False),
        4: ("C", 0, False),
        5: ("O", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
        8: ("H", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
        12: ("H", 0, None),
    },
    {
        frozenset({0, 6}): (1, None),
        frozenset({1, 7}): (1, None),
        frozenset({1, 2}): (1, None),
        frozenset({0, 3}): (1, None),
        frozenset({4, 5}): (1, None),
        frozenset({2, 10}): (1, None),
        frozenset({4, 12}): (1, None),
        frozenset({0, 4}): (1, None),
        frozenset({2, 4}): (1, None),
        frozenset({8, 1}): (0.9, None),
        frozenset({3, 5}): (1, None),
        frozenset({11, 3}): (1, None),
        frozenset({1, 3}): (1, None),
        frozenset({9, 2}): (1, None),
        frozenset({0, 8}): (0.1, None),
    },
)


# Simple Bridgehead Transfer
C5H7_TSG = (
    {
        0: ("C", 0, True),
        1: ("C", 0, None),
        2: ("H", 0, None),
        3: ("C", 0, False),
        4: ("H", 0, None),
        5: ("C", 0, None),
        6: ("H", 0, None),
        7: ("H", 0, None),
        8: ("C", 0, None),
        9: ("H", 0, None),
        10: ("H", 0, None),
        11: ("H", 0, None),
    },
    {
        frozenset({1, 4}): (1, None),
        frozenset({8, 11}): (1, None),
        frozenset({0, 6}): (0.1, None),
        frozenset({3, 7}): (1, None),
        frozenset({0, 1}): (1, None),
        frozenset({0, 2}): (1, None),
        frozenset({3, 6}): (0.9, None),
        frozenset({3, 5}): (1, None),
        frozenset({10, 5}): (1, None),
        frozenset({8, 5}): (1, None),
        frozenset({1, 3}): (1, None),
        frozenset({9, 5}): (1, None),
        frozenset({0, 8}): (1, None),
    },
)
C5H7_TS_GEO = (
    ("C", (2.411517, 0.003628, 0.0)),
    ("C", (0.9272, -2.162191, 0.0)),
    ("H", (4.474794, 0.007151, -0.0)),
    ("C", (-1.830041, -1.468223, -0.0)),
    ("H", (1.626336, -4.101193, 0.0)),
    ("C", (-1.834561, 1.463071, 0.0)),
    ("H", (-2.817187, -2.252453, -1.662617)),
    ("H", (-2.81719, -2.252455, 1.662611)),
    ("C", (0.920338, 2.16457, -0.0)),
    ("H", (-2.82419, 2.244254, 1.66257)),
    ("H", (-2.824192, 2.244256, -1.662566)),
    ("H", (1.614892, 4.105309, -2e-06)),
)


@pytest.mark.parametrize(
    "formula,tsg,geo,npars1,npars2",
    [
        ("CH4CLFNO", CH4CLFNO_TSG, CH4CLFNO_GEO, 1, 1),
        ("C4H11O2", C4H11O2_TSG, C4H11O2_GEO, 1, 1),
        ("C4H9O3", C4H9O3_TSG, C4H9O3_GEO, 2, 2),
        ("C2H3F2O", C2H3F2O_TSG, C2H3F2O_GEO, 2, 2),
        ("C4H5F3O2", C4H5F3O2_TSG, C4H5F3O2_GEO, 3, 4),
        ("C4H9O2", C4H9O2_TSG, C4H9O2_GEO, 2, 2),
        ("C2H4O2", C2H4O2_TSG, C2H4O2_GEO, 1, 1),
        ("C4H4F2", C4H4F2_TSG, C4H4F2_GEO, 2, 2),
        ("C8H14", C8H14_TSG, C8H14_GEO, 2, 2),
    ],
)
def test__set_stereo_from_geometry(formula, tsg, geo, npars1, npars2):
    """test graph.set_stereo_from_geometry"""

    print(f"{formula}: testing set_stereo_from_geometry")
    ftsg = graph.without_stereo(tsg)
    rtsg = graph.ts_reverse(ftsg)
    assert ftsg != rtsg
    # Check that they have the same stereogenic keys
    fste_keys = graph.unassigned_stereocenter_keys(ftsg)
    rste_keys = graph.unassigned_stereocenter_keys(rtsg)
    print(fste_keys)
    print(rste_keys)
    assert len(fste_keys) == len(rste_keys) == npars1
    assert fste_keys == rste_keys
    # Check that they have the same parities
    ftsg = graph.set_stereo_from_geometry(ftsg, geo)
    rtsg = graph.set_stereo_from_geometry(rtsg, geo)
    fste_par_dct = automol.util.dict_.filter_by_value(
        graph.stereo_parities(ftsg), lambda x: x is not None
    )
    rste_par_dct = automol.util.dict_.filter_by_value(
        graph.stereo_parities(rtsg), lambda x: x is not None
    )
    print(fste_par_dct)
    print(rste_par_dct)
    assert len(fste_par_dct) == len(rste_par_dct) == npars2
    assert fste_par_dct == rste_par_dct


def test__to_local_stereo():
    """test graph.to_local_stereo"""
    print("C4H11O2: testing to_local_stereo")
    # Swap two neighbors to reverse local parity in one copy
    tsg1 = C4H11O2_TSG
    tsg2 = graph.relabel(C4H11O2_TSG, {0: 4, 4: 0})
    # Convert both to local stereo
    loc_tsg1 = graph.to_local_stereo(tsg1)
    loc_tsg2 = graph.to_local_stereo(tsg2)
    # Check that the local parities have opposite values
    loc_par1 = graph.stereo_parities(loc_tsg1)[2]
    loc_par2 = graph.stereo_parities(loc_tsg2)[2]
    assert loc_par1 is False
    assert loc_par2 is True

    print("CH4ClFNO: testing to_local_stereo")
    # Reverse direction of an Sn2 reaction
    tsg1 = CH4CLFNO_TSG
    tsg2 = graph.ts.reverse(tsg1)
    # Convert both to local stereo
    loc_tsg1 = graph.to_local_stereo(tsg1)
    loc_tsg2 = graph.to_local_stereo(tsg2)
    # Check that the local parities have opposite values
    loc_par1 = graph.stereo_parities(loc_tsg1)[0]
    loc_par2 = graph.stereo_parities(loc_tsg2)[0]
    assert loc_par1 is False
    assert loc_par2 is True


def test__from_local_stereo():
    """test graph.from_local_stereo"""

    def _test(formula, tsg):
        print(f"{formula}: testing from_local_stereo")
        loc_tsg = graph.to_local_stereo(tsg)
        assert tsg == graph.from_local_stereo(loc_tsg), f"{tsg}\n\n{loc_tsg}"

    _test("CH4CLFNO", CH4CLFNO_TSG)
    _test("CH4CLFNO(rev)", graph.ts.reverse(CH4CLFNO_TSG))
    _test("C4H11O2", C4H11O2_TSG)
    # _test("C2H3O4", C2H3O4_TSG)
    _test("C4H9O3", C4H9O3_TSG)
    _test("C2H3F2O", C2H3F2O_TSG)
    _test("C4H5F3O2", C4H5F3O2_TSG)
    _test("C4H9O2", C4H9O2_TSG)
    _test("C2H4O2", C2H4O2_TSG)
    _test("C4H4F2", C4H4F2_TSG)
    _test("C8H14", C8H14_TSG)


@pytest.mark.parametrize(
    "formula,ts_gra,ts_counts",
    [
        ("CH4CLFNO", CH4CLFNO_TSG, [1, 1]),
        ("C4H11O2", C4H11O2_TSG, [2]),
        ("C2H3O4", C2H3O4_TSG, [1]),
        ("C4H9O3", C4H9O3_TSG, [2, 2]),
        ("C2H3F2O", C2H3F2O_TSG, [1, 1, 1, 1]),
        ("C4H5F3O2", C4H5F3O2_TSG, [1] * 12),
        ("C4H9O2", C4H9O2_TSG, [1, 1, 1, 1]),
        ("C2H4O2", C2H4O2_TSG, [1, 1]),
        ("C4H4F2", C4H4F2_TSG, [1, 1, 1, 1]),
        ("C8H14", C8H14_TSG, [1, 1, 1, 1]),
        ("C5H7O2", C5H7O2_TSG, [1, 1]),
        ("C5H7O3", C5H7O3_TSG, [1, 1, 1, 1]),
        ("C5H7", C5H7_TSG, [2]),
    ],
)
def test__ts__expand_reaction_stereo(formula, ts_gra, ts_counts):
    """test graph.ts.expand_reaction_stereo"""
    print(f"{formula}: testing ts.expand_reaction_stereo")
    _, _, ts_sgras_lst = zip(*graph.ts.expand_reaction_stereo(ts_gra))
    print(list(map(len, ts_sgras_lst)))
    assert list(map(len, ts_sgras_lst)) == ts_counts


def test__ts__fleeting_stereocenter_keys():
    """test graph.ts.fleeting_stereocenter_keys"""

    def _test(formula, tsg, ref_keys, ref_ns_keys):
        print(f"{formula}: testing ts.fleeting_stereocenter_keys")
        keys = graph.ts.fleeting_stereocenter_keys(tsg)
        assert keys == ref_keys, f"{keys} != {ref_keys}"

        ns_keys = graph.ts.fleeting_stereocenter_keys(tsg, strict=False)
        assert ns_keys == ref_ns_keys, f"{ns_keys} != {ref_ns_keys}"

        # These examples all only have fleeting atom stereo
        if keys:
            assert graph.ts.has_fleeting_atom_or_bond_stereo(tsg) == (True, False)
        else:
            assert graph.ts.has_fleeting_atom_or_bond_stereo(tsg) == (False, False)

    _test("CH4CLFNO", CH4CLFNO_TSG, frozenset(), frozenset())
    _test("C4H11O2", C4H11O2_TSG, frozenset({2}), frozenset({2}))
    _test("C2H3O4", C2H3O4_TSG, frozenset(), frozenset())
    _test("C4H9O3", C4H9O3_TSG, frozenset({2}), frozenset({2}))
    _test("C2H3F2O", C2H3F2O_TSG, frozenset(), frozenset())
    _test("C4H5F3O2", C4H5F3O2_TSG, frozenset(), frozenset())
    _test("C4H9O2", C4H9O2_TSG, frozenset(), frozenset({2}))
    _test("C7H15O2", C7H15O2_TSG, frozenset({2, 3}), frozenset({2, 3}))
    _test("C2H4O2", C2H4O2_TSG, frozenset(), frozenset())
    _test("C4H4F2", C4H4F2_TSG, frozenset(), frozenset())
    _test("C6H11O2", C6H11O2_TSG, frozenset(), frozenset())
    _test("C6H11O2b", C6H11O2b_TSG, frozenset(), frozenset())
    _test("C8H14", C8H14_TSG, frozenset(), frozenset())


@pytest.mark.parametrize(
    "formula, tsg, rcts_par_dct_ref, prds_par_dct_ref",
    [
        ("CH4CLFNO", CH4CLFNO_TSG, {0: False}, {0: True}),
        ("C4H11O2", C4H11O2_TSG, {}, {}),
        ("C2H3O4", C2H3O4_TSG, {}, {}),
        ("C4H9O3", C4H9O3_TSG, {3: True}, {3: True}),
        ("C2H3F2O", C2H3F2O_TSG, {frozenset({0, 1}): True}, {0: False}),
        (
            "C4H5F3O2",
            C4H5F3O2_TSG,
            {2: False, 3: True, frozenset({0, 1}): False},
            {0: False, 2: False, 3: True},
        ),
        ("C4H9O2", C4H9O2_TSG, {3: True}, {frozenset({2, 3}): True}),
        ("C2H4O2", C2H4O2_TSG, {frozenset({0, 1}): False}, {frozenset({0, 1}): True}),
        (
            "C4H4F2",
            C4H4F2_TSG,
            {frozenset({1, 2}): True, frozenset({6, 7}): True},
            {frozenset({1, 2}): False, frozenset({6, 7}): False},
        ),
        ("C6H9", C6H9_TSG, {frozenset({2, 4}): False}, {}),
        ("C6H11O2", C6H11O2_TSG, {5: False, frozenset({3, 4}): False}, {5: False}),
        ("C6H11O2b", C6H11O2b_TSG, {2: False}, {2: False, frozenset({0, 1}): True}),
        ("C8H14", C8H14_TSG, {frozenset({1, 6}): True}, {1: True, 6: True}),
        ("C5H7O", C5H7O_TSG, {}, {frozenset({2, 3}): True}),
        ("C5H7O_B", C5H7O_B_TSG, {3: True, 4: False}, {3: False, 4: True}),
    ],
)
def test__ts__reagents_graph(formula, tsg, rcts_par_dct_ref, prds_par_dct_ref):
    """test graph.ts.reactants_graph and graph.ts.products_graph"""

    print(f"{formula}: testing reactants_graph")
    rcts_gra = automol.graph.ts.reactants_graph(tsg)
    rcts_par_dct = automol.util.dict_.filter_by_value(
        automol.graph.stereo_parities(rcts_gra), lambda x: x is not None
    )
    print(f"asserting {rcts_par_dct} == {rcts_par_dct_ref}")
    assert rcts_par_dct == rcts_par_dct_ref

    print(f"{formula}: testing products_graph")
    prds_gra = automol.graph.ts.products_graph(tsg)
    prds_par_dct = automol.util.dict_.filter_by_value(
        automol.graph.stereo_parities(prds_gra), lambda x: x is not None
    )
    print(f"asserting {prds_par_dct} == {prds_par_dct_ref}")
    assert prds_par_dct == prds_par_dct_ref

    print("Re-testing with implicit graphs...")
    tsg = automol.graph.implicit(tsg)

    print(f"{formula}: testing reactants_graph")
    rcts_gra = automol.graph.ts.reactants_graph(tsg)
    rcts_par_dct = automol.util.dict_.filter_by_value(
        automol.graph.stereo_parities(rcts_gra), lambda x: x is not None
    )
    print(f"asserting {rcts_par_dct} == {rcts_par_dct_ref}")
    assert rcts_par_dct == rcts_par_dct_ref

    print(f"{formula}: testing products_graph")
    prds_gra = automol.graph.ts.products_graph(tsg)
    prds_par_dct = automol.util.dict_.filter_by_value(
        automol.graph.stereo_parities(prds_gra), lambda x: x is not None
    )
    print(f"asserting {prds_par_dct} == {prds_par_dct_ref}")
    assert prds_par_dct == prds_par_dct_ref

    print("---")


@pytest.mark.parametrize(
    "formula,tsg",
    [
        ("CH4CLFNO", CH4CLFNO_TSG),
        ("C4H11O2", C4H11O2_TSG),
        ("C2H3O4", C2H3O4_TSG),
        ("C4H9O3", C4H9O3_TSG),
        ("C2H3F2O", C2H3F2O_TSG),
        ("C4H5F3O2", C4H5F3O2_TSG),
        ("C4H9O2", C4H9O2_TSG),
        ("C2H4O2", C2H4O2_TSG),
        ("C4H4F2", C4H4F2_TSG),
        ("C8H14", C8H14_TSG),
    ],
)
def test__amchi(formula, tsg):
    """test graph.amchi"""
    print(f"{formula}: testing amchi")
    ftsg = tsg
    rtsg = graph.ts.reverse(ftsg)

    fchi = graph.amchi(ftsg)
    rchi = graph.amchi(rtsg)

    # Check interconversion
    print(fchi)
    print(graph.amchi(automol.amchi.graph(fchi)))
    assert graph.amchi(automol.amchi.graph(fchi)) == fchi
    assert graph.amchi(automol.amchi.graph(rchi)) == rchi

    print("fchi:", fchi)
    print("rchi:", rchi)
    assert fchi[:-1] == rchi[:-1]

    orig_keys = sorted(graph.atom_keys(tsg))
    for _ in range(5):
        perm_keys = numpy.random.permutation(orig_keys)
        perm_dct = dict(zip(orig_keys, perm_keys))

        perm_ftsg = graph.relabel(ftsg, perm_dct)
        perm_rtsg = graph.relabel(rtsg, perm_dct)

        perm_fchi = graph.amchi(perm_ftsg)
        perm_rchi = graph.amchi(perm_rtsg)

        assert perm_fchi == fchi, f"\n{perm_fchi} !=\n{fchi}"
        assert perm_rchi == rchi, f"\n{perm_rchi} !=\n{rchi}"


def test__linear_atom_keys():
    """test graph.linear_atom_keys"""
    tsg = C5H4_TSG
    rgra = graph.ts.reactants_graph(tsg, stereo=False)
    pgra = graph.ts.products_graph(tsg, stereo=False)
    assert graph.linear_atom_keys(tsg) == frozenset({1, 4, 5, 6})
    assert graph.linear_atom_keys(rgra) == frozenset({4, 5, 6})
    assert graph.linear_atom_keys(pgra) == frozenset({1, 4, 5, 6})


@pytest.mark.parametrize(
    "formula, tsg, extend, lin_seg_keys0",
    [
        ("C4H6", C4H6_TSG, False, ((1, 2),)),
        ("C4H6", C4H6_TSG, True, ((0, 1, 2, 6),)),
        ("C5H4", C5H4_TSG, False, ((1, 4, 5, 6),)),
        ("C5H4", C5H4_TSG, True, ((0, 1, 4, 5, 6, 7),)),
        ("C3H5O-1", C3H5O_TSG1, False, ((0, 1),)),
        ("C3H5O-1", C3H5O_TSG1, True, ((3, 0, 1, 2),)),
        ("C3H5O-2", C3H5O_TSG2, False, ((1, 2),)),
        ("C3H5O-2", C3H5O_TSG2, True, ((0, 1, 2, 3),)),
        ("CH4CLFNO", CH4CLFNO_TSG, False, ((0,),)),
        ("CH4CLFNO", CH4CLFNO_TSG, True, ((1, 0, 7),)),
    ],
)
def test__linear_segments_atom_keys(formula, tsg, extend, lin_seg_keys0):
    """Test graph.rotational_symmetry_number."""
    print(f"{formula}: testing rotational_symmetry_number")
    lin_seg_keys = graph.linear_segments_atom_keys(tsg, extend=extend)
    assert lin_seg_keys == lin_seg_keys0, f"{lin_seg_keys} != {lin_seg_keys0}"


@pytest.mark.parametrize(
    "formula, tsg, bkey, sym_num0",
    [
        ("C4H6", C4H6_TSG, (2, 6), 3),
        ("C3H5O-1", C3H5O_TSG1, (0, 5), 3),
        ("C3H5O-2", C3H5O_TSG2, (2, 3), 1),
    ],
)
def test__rotational_symmetry_number(formula, tsg, bkey, sym_num0):
    """Test graph.rotational_symmetry_number."""
    print(f"{formula}: testing rotational_symmetry_number")
    sym_num = graph.rotational_symmetry_number(tsg, *bkey)
    assert sym_num == sym_num0, f"{sym_num} != {sym_num0}"


def test__rotational_bond_keys():
    """test graph.rotational_bond_keys"""
    # C=C(O[O])C=C => [CH]=C(OO)C=C
    print("Rotational bond keys for C=C(O[O])C=C => [CH]=C(OO)C=C")
    tsg = (
        {
            0: ("C", 0, None),
            1: ("C", 0, None),
            2: ("C", 0, None),
            3: ("C", 0, None),
            4: ("O", 0, None),
            5: ("O", 0, None),
            6: ("H", 0, None),
            7: ("H", 0, None),
            8: ("H", 0, None),
            9: ("H", 0, None),
            10: ("H", 0, None),
        },
        {
            frozenset({1, 9}): (1, None),
            frozenset({2, 10}): (1, None),
            frozenset({0, 6}): (1, None),
            frozenset({2, 3}): (1, None),
            frozenset({4, 5}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({1, 8}): (0.9, None),
            frozenset({3, 5}): (1, None),
            frozenset({0, 7}): (1, None),
            frozenset({1, 3}): (1, False),
            frozenset({4, 8}): (0.1, None),
        },
    )
    rtsg = graph.ts.reverse(tsg)
    assert graph.rotational_bond_keys(tsg) == frozenset({frozenset({2, 3})})
    assert graph.rotational_bond_keys(rtsg) == frozenset({frozenset({2, 3})})

    # C=C(O[O])C=C => C=C(OO)[C]=C
    # (rotational bond is locked by TS ring)
    print("Rotational bond keys for C=C(O[O])C=C => C=C(OO)[C]=C")
    tsg = (
        {
            0: ("C", 0, None),
            1: ("C", 0, None),
            2: ("C", 0, None),
            3: ("C", 0, None),
            4: ("O", 0, None),
            5: ("O", 0, None),
            6: ("H", 0, None),
            7: ("H", 0, None),
            8: ("H", 0, None),
            9: ("H", 0, None),
            10: ("H", 0, None),
        },
        {
            frozenset({4, 10}): (0.1, None),
            frozenset({1, 9}): (1, None),
            frozenset({0, 6}): (1, None),
            frozenset({2, 3}): (1, None),
            frozenset({2, 10}): (0.9, None),
            frozenset({4, 5}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({1, 8}): (1, None),
            frozenset({3, 5}): (1, None),
            frozenset({0, 7}): (1, None),
            frozenset({1, 3}): (1, None),
        },
    )
    rtsg = graph.ts.reverse(tsg)
    assert graph.rotational_bond_keys(tsg) == frozenset()
    assert graph.rotational_bond_keys(rtsg) == frozenset()


def test__radical_atom_keys():
    """test graph.radical_atom_keys"""
    # [CH2]CC=C + [OH] => [CH2]C[CH]CO
    tsg = (
        {
            0: ("C", 0, None),
            1: ("H", 0, None),
            2: ("H", 0, None),
            3: ("C", 0, None),
            4: ("C", 0, None),
            5: ("C", 0, None),
            6: ("H", 0, None),
            7: ("H", 0, None),
            8: ("H", 0, None),
            9: ("H", 0, None),
            10: ("H", 0, None),
            11: ("O", 0, None),
            12: ("H", 0, None),
        },
        {
            frozenset({3, 4}): (1, None),
            frozenset({5, 11}): (0.1, None),
            frozenset({11, 12}): (1, None),
            frozenset({3, 7}): (1, None),
            frozenset({0, 3}): (1, None),
            frozenset({4, 5}): (1, None),
            frozenset({0, 1}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({3, 6}): (1, None),
            frozenset({4, 8}): (1, None),
            frozenset({5, 10}): (1, None),
            frozenset({5, 9}): (1, None),
        },
    )
    print(automol.graph.radical_atom_keys(tsg))
    assert automol.graph.radical_atom_keys(tsg) == frozenset({0})


def test__geometry():
    """test graph.geometry"""
    # C#C + [CH2] => [CH]=C[CH2]
    ts_gra = (
        {
            0: ("C", 0, None),
            1: ("H", 0, None),
            2: ("H", 0, None),
            3: ("C", 0, None),
            4: ("C", 0, None),
            5: ("H", 0, None),
            6: ("H", 0, None),
        },
        {
            frozenset({3, 4}): (1, True),
            frozenset({4, 6}): (1, None),
            frozenset({3, 5}): (1, None),
            frozenset({0, 1}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({0, 4}): (0.1, None),
        },
    )
    rct_geos = (
        (
            ("C", (0.0, 0.0, 0.20515244763195276)),
            ("H", (0.0, 1.8754738110196252, -0.6154554531697327)),
            ("H", (-0.0, -1.8754738110196252, -0.6154554531697327)),
        ),
        (
            ("C", (0.0, 0.0, 1.1364019233530684)),
            ("C", (0.0, 0.0, -1.1364019233530684)),
            ("H", (0.0, 0.0, 3.152539388247273)),
            ("H", (0.0, 0.0, -3.152539388247273)),
        ),
    )
    geo_idx_dct = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
    ts_geo = graph.ts_geometry_from_reactants(
        ts_gra, rct_geos, geo_idx_dct=geo_idx_dct, check=True
    )
    print(automol.geom.round_(ts_geo))

    # [CH2]CC=C[CH2] => C=CC=C[CH2] + [H]
    ts_gra = (
        {
            0: ("C", 0, None),
            1: ("C", 0, False),
            2: ("H", 0, None),
            3: ("H", 0, None),
            4: ("C", 0, None),
            5: ("H", 0, None),
            6: ("C", 0, None),
            7: ("H", 0, None),
            8: ("C", 0, None),
            9: ("H", 0, None),
            10: ("H", 0, None),
            11: ("H", 0, None),
            12: ("H", 0, None),
        },
        {
            frozenset({9, 6}): (1, None),
            frozenset({1, 4}): (1, True),
            frozenset({4, 6}): (1, True),
            frozenset({8, 11}): (1, None),
            frozenset({0, 3}): (1, None),
            frozenset({0, 1}): (1, None),
            frozenset({0, 2}): (1, None),
            frozenset({8, 10}): (1, None),
            frozenset({1, 12}): (0.9, None),
            frozenset({1, 5}): (1, None),
            frozenset({8, 6}): (1, None),
            frozenset({4, 7}): (1, None),
        },
    )
    rct_geos = (
        (
            ("C", (-4.672291, -0.457455, -0.152981)),
            ("C", (-2.307244, 0.626363, -0.386708)),
            ("H", (-4.95983, -2.139761, 0.988667)),
            ("H", (-6.310429, 0.3332, -1.096682)),
            ("C", (-0.126763, -0.273254, 0.751996)),
            ("H", (-2.136558, 2.316753, -1.554504)),
            ("C", (2.427578, 0.917615, 0.483968)),
            ("H", (-0.245351, -1.969796, 1.913306)),
            ("C", (4.36715, -0.817943, -0.625343)),
            ("H", (2.260094, 2.642939, -0.678196)),
            ("H", (3.103461, 1.582897, 2.340003)),
            ("H", (6.355744, -0.318011, -0.564376)),
            ("H", (3.802295, -2.420165, -1.773812)),
        ),
    )
    geo_idx_dct = {
        0: 8,
        1: 6,
        2: 11,
        3: 12,
        4: 4,
        5: 10,
        6: 1,
        7: 7,
        8: 0,
        9: 5,
        10: 2,
        11: 3,
        12: 9,
    }
    ts_geo = graph.ts_geometry_from_reactants(
        ts_gra, rct_geos, geo_idx_dct=geo_idx_dct, check=True
    )
    print(automol.geom.round_(ts_geo))

    # Bridgehead Atom Stereo Pair
    # [O]OC1OC=CCC1 => OOC1OC=CC[CH]1
    #     *      *       *       *
    # [* marks the pair of bridgehead stereo atoms]
    # One of the stereoatom neighbors is undergoing H migration, which presents a tricky
    # case for stereo correction
    ts_gra = (
        {
            0: ("O", 0, None),
            1: ("O", 0, None),
            2: ("C", 0, False),
            3: ("O", 0, None),
            4: ("C", 0, None),
            5: ("C", 0, None),
            6: ("C", 0, None),
            7: ("C", 0, False),
            8: ("H", 0, None),
            9: ("H", 0, None),
            10: ("H", 0, None),
            11: ("H", 0, None),
            12: ("H", 0, None),
            13: ("H", 0, None),
            14: ("H", 0, None),
        },
        {
            frozenset({8, 2}): (1, None),
            frozenset({2, 3}): (1, None),
            frozenset({4, 5}): (1, None),
            frozenset({0, 1}): (1, None),
            frozenset({6, 7}): (1, None),
            frozenset({14, 7}): (1, None),
            frozenset({10, 5}): (1, None),
            frozenset({9, 4}): (1, None),
            frozenset({2, 7}): (1, None),
            frozenset({11, 6}): (1, None),
            frozenset({3, 4}): (1, None),
            frozenset({0, 13}): (0.1, None),
            frozenset({13, 7}): (0.9, None),
            frozenset({1, 2}): (1, None),
            frozenset({5, 6}): (1, None),
            frozenset({12, 6}): (1, None),
        },
    )
    rct_geos = [
        (
            ("O", (-5.522426, -3.187961, -0.390535)),
            ("O", (-3.762872, -2.550969, 1.204516)),
            ("C", (-1.711414, -1.434957, 0.016648)),
            ("O", (0.534461, -2.85378, 0.309636)),
            ("C", (2.861053, -1.636541, 0.08201)),
            ("C", (3.043655, 0.805284, -0.391908)),
            ("C", (0.751244, 2.339971, -0.699047)),
            ("C", (-1.259943, 1.206692, 1.044995)),
            ("H", (-2.154465, -1.19976, -2.033012)),
            ("H", (4.600758, -2.688934, 0.293275)),
            ("H", (4.868492, 1.666952, -0.552345)),
            ("H", (0.020412, 1.995037, -2.666963)),
            ("H", (0.9976, 4.339791, -0.264024)),
            ("H", (-2.96238, 2.361521, 1.108978)),
            ("H", (-0.374629, 0.994683, 2.963635)),
        )
    ]
    ts_geo = graph.ts_geometry_from_reactants(ts_gra, rct_geos, check=True)
    print(automol.geom.round_(ts_geo))


def test__zmatrix():
    """test z-matrix generation"""
    # C[C]1O[C@H]1COO >> C/C([O])=C\COO
    ts_gra = (
        {
            0: ("C", 0, None),
            1: ("C", 0, False),
            2: ("C", 0, None),
            3: ("C", 0, None),
            4: ("O", 0, None),
            5: ("O", 0, None),
            6: ("O", 0, None),
            7: ("H", 0, None),
            8: ("H", 0, None),
            9: ("H", 0, None),
            10: ("H", 0, None),
            11: ("H", 0, None),
            12: ("H", 0, None),
            13: ("H", 0, None),
        },
        {
            frozenset({3, 4}): (1, None),
            frozenset({1, 4}): (0.9, None),
            frozenset({1, 2}): (1, None),
            frozenset({0, 3}): (1, None),
            frozenset({2, 6}): (1, None),
            frozenset({2, 12}): (1, None),
            frozenset({5, 6}): (1, None),
            frozenset({1, 3}): (1, True),
            frozenset({13, 5}): (1, None),
            frozenset({2, 11}): (1, None),
            frozenset({0, 7}): (1, None),
            frozenset({1, 10}): (1, None),
            frozenset({0, 8}): (1, None),
            frozenset({0, 9}): (1, None),
        },
    )

    ts_geo = automol.graph.geometry(ts_gra)
    ts_zma, ts_zc = automol.geom.zmatrix_with_conversion_info(ts_geo, gra=ts_gra)
    ts_zgra = automol.graph.apply_zmatrix_conversion(ts_gra, ts_zc)

    automol.graph.zmatrix_matches(ts_zgra, ts_zma)


@pytest.mark.parametrize(
    "tsg, bh_dct0",
    [(C5H7O_B_TSG, {(1, 0): ((2, 3, 8), (4, 3, 8)), (3, 4): ((0, 1, 5), (0, 2, 5))})],
)
def test__stereoatom_bridgehead_pairs(tsg, bh_dct0):
    """Test bridgehead pair detection."""
    bh_dct = automol.graph.stereoatom_bridgehead_pairs(tsg)
    assert bh_dct == bh_dct0, f"{bh_dct} != {bh_dct0}"


@pytest.mark.parametrize("tsg0, ts_geo0", [(C5H7_TSG, C5H7_TS_GEO)])
def test__stereo_corrected_geometry(tsg0, ts_geo0):
    # Invert the stereochemistry so that correction of atoms is required
    ts_geo0 = automol.geom.reflect_coordinates(ts_geo0)

    ts_geo = automol.graph.stereo_corrected_geometry(tsg0, ts_geo0)
    tsg = automol.graph.set_stereo_from_geometry(tsg0, ts_geo)
    assert tsg == tsg0, f"   {tsg}\n!= {tsg0}"


if __name__ == "__main__":
    # test__set_stereo_from_geometry()
    # test__to_local_stereo()
    # test__from_local_stereo()
    # test__ts__reagents_graph()
    # test__rotational_bond_keys()
    # test__ts__expand_reaction_stereo()
    # test__amchi("CH4CLFNO", CH4CLFNO_TSG)
    # test__ts__fleeting_stereocenter_keys()
    # test__linear_atom_keys()
    # test__radical_atom_keys()
    # test__geometry()
    # test__zmatrix()
    # test__ts__expand_reaction_stereo("C5H6O", C5H6O_TSG, [1, 1])
    # test__ts__expand_reaction_stereo("C5H7O3", C5H7O3_TSG, [1, 1, 1, 1])
    # test__ts__reagents_graph(
    #     "C5H7O_B", C5H7O_B_TSG, {3: True, 4: False}, {3: False, 4: True}
    # )
    # test__stereoatom_bridgehead_pairs(
    #     C5H7O_B_TSG, {(1, 0): ((2, 3, 8), (4, 3, 8)), (3, 4): ((0, 1, 5), (0, 2, 5))}
    # )
    # test__stereo_corrected_geometry(C5H7_TSG, C5H7_TS_GEO)
    test__ts__expand_reaction_stereo("C5H7", C5H7_TSG, [2])
