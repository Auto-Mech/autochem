""" test automol.graph.ts
"""
# import numpy
# from automol import graph

# Fleeting Atom Stereo
# CCOCC + [OH] => C[CH]OCC + O
#  *
# [* marks a fleeting TS stereosite]
C4H11O2_TSG = (
 {0: ('C', 0, None),
  1: ('C', 0, None),
  2: ('C', 0, False),
  3: ('C', 0, None),
  4: ('O', 0, None),
  5: ('H', 0, None),
  6: ('H', 0, None),
  7: ('H', 0, None),
  8: ('H', 0, None),
  9: ('H', 0, None),
  10: ('H', 0, None),
  11: ('H', 0, None),
  12: ('H', 0, None),
  13: ('H', 0, None),
  14: ('H', 0, None),
  15: ('O', 0, None),
  16: ('H', 0, None)},
 {frozenset({0, 5}): (1, None),
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
  frozenset({1, 3}): (1, None)})
C4H11O2_GEO = (
 ('C', (4.33819, 1.441087, 0.295484)),
 ('C', (-4.574912, 0.214792, -0.174836)),
 ('C', (2.316425, -0.587509, 0.260823)),
 ('C', (-2.077647, -1.187579, -0.205727)),
 ('O', (-0.0881, 0.60958, 0.056152)),
 ('H', (6.220916, 0.603391, 0.451236)),
 ('H', (4.047672, 2.732272, 1.886248)),
 ('H', (4.255128, 2.582155, -1.4288)),
 ('H', (-6.159774, -1.096745, -0.371638)),
 ('H', (-4.65463, 1.594987, -1.714572)),
 ('H', (-4.796227, 1.270025, 1.591399)),
 ('H', (2.682893, -1.715305, -1.436134)),
 ('H', (2.341862, -1.826667, 1.919052)),
 ('H', (-2.01011, -2.538954, 1.361152)),
 ('H', (-1.866407, -2.209146, -1.993927)),
 ('O', (3.7089, -3.401617, -3.726459)),
 ('H', (3.282065, -5.016675, -2.48937)))

# Fleeting Bond Stereo
# C=C(O[O])OO => [CH]=C(OO)OO
#  *
# [* marks a fleeting TS stereosite]
C2H3O4_TSG = (
 {0: ('C', 0, None),
  1: ('C', 0, None),
  2: ('O', 0, None),
  3: ('O', 0, None),
  4: ('O', 0, None),
  5: ('O', 0, None),
  6: ('H', 0, None),
  7: ('H', 0, None),
  8: ('H', 0, None)},
 {frozenset({2, 8}): (1, None),
  frozenset({1, 4}): (1, None),
  frozenset({0, 6}): (0.9, None),
  frozenset({0, 1}): (1, False),
  frozenset({3, 6}): (0.1, None),
  frozenset({2, 4}): (1, None),
  frozenset({1, 5}): (1, None),
  frozenset({3, 5}): (1, None),
  frozenset({0, 7}): (1, None)})
C2H3O4_GEO = (
 ('C', (-2.812896, 0.535324, 0.864225)),
 ('C', (-0.773943, -0.715224, 0.054102)),
 ('O', (3.585469, -2.437992, -0.732028)),
 ('O', (0.896435, 2.730653, -1.397099)),
 ('O', (1.399056, -2.02279, 0.305438)),
 ('O', (-0.803513, 0.372439, -2.38697)),
 ('H', (-2.22919, 0.44838, 2.621235)),
 ('H', (-4.102653, 2.150414, 1.056667)),
 ('H', (4.841235, -1.061205, -0.38557)))

# Fleeting Atom Stereo + Conserved Atom Stereo
# CCO[C@H](O[O])C => C[CH]O[C@H](OO)C
#  *
# [* marks a fleeting stereosite]
C4H9O3_TSG = (
 {0: ('C', 0, None),
  1: ('C', 0, None),
  2: ('C', 0, False),
  3: ('C', 0, True),
  4: ('O', 0, None),
  5: ('O', 0, None),
  6: ('O', 0, None),
  7: ('H', 0, None),
  8: ('H', 0, None),
  9: ('H', 0, None),
  10: ('H', 0, None),
  11: ('H', 0, None),
  12: ('H', 0, None),
  13: ('H', 0, None),
  14: ('H', 0, None),
  15: ('H', 0, None)},
 {frozenset({4, 6}): (1, None),
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
  frozenset({1, 3}): (1, None)})
C4H9O3_GEO = (
 ('C', (-3.858804, -0.017433, -0.816915)),
 ('C', (3.284911, 2.431032, 0.462603)),
 ('C', (-1.959761, -1.667044, 0.44462)),
 ('C', (2.02811, -0.049194, 0.203738)),
 ('O', (1.521606, -4.246381, 1.551043)),
 ('O', (-0.22096, 0.069541, 1.635187)),
 ('O', (3.101169, -2.349906, 0.925334)),
 ('H', (-3.764188, 2.014997, -0.161893)),
 ('H', (-5.805311, -0.707311, -0.433068)),
 ('H', (-3.483251, 0.040265, -2.87643)),
 ('H', (3.006888, 3.239081, 2.415149)),
 ('H', (2.436056, 3.760618, -0.93411)),
 ('H', (5.336483, 2.305681, 0.107293)),
 ('H', (-1.992606, -3.145298, 1.82341)),
 ('H', (-1.040456, -2.832608, -1.161927)),
 ('H', (1.546318, -0.226049, -1.899431)))

# (Reactant Bond Stereo => Product Atom Stereo)
# F/C=C/F + [OH] => F[CH][C@H](O)F
C2H3F2O_TSG = (
 {0: ('C', 0, False),
  1: ('C', 0, None),
  2: ('F', 0, None),
  3: ('F', 0, None),
  4: ('H', 0, None),
  5: ('H', 0, None),
  6: ('O', 0, None),
  7: ('H', 0, None)},
 {frozenset({6, 7}): (1, None),
  frozenset({1, 5}): (1, None),
  frozenset({0, 6}): (0.1, None),
  frozenset({0, 1}): (1, True),
  frozenset({0, 2}): (1, None),
  frozenset({1, 3}): (1, None),
  frozenset({0, 4}): (1, None)})
C2H3F2O_GEO = (
 ('C', (1.297288, 0.39764, 0.454963)),
 ('C', (-1.326186, -0.366377, 0.15782)),
 ('F', (3.477276, -1.308171, 0.493162)),
 ('F', (-3.438916, 1.322689, -0.431268)),
 ('H', (2.079235, 2.259128, 0.133371)),
 ('H', (-2.02146, -2.249426, -0.230421)),
 ('O', (1.234258, 0.222229, -4.125039)),
 ('H', (0.617002, 2.080766, -3.428028)))

# Conserved Atom Stereo + (Reactant Bond Stereo => Product Atom Stereo)
# F/C=C([C@@H](F)O)\[C@H](F)O + [OH] => F[C@H]([C]([C@@H](F)O)[C@H](F)O)O
C4H5F3O2_TSG = (
 {0: ('C', 0, False),
  1: ('C', 0, None),
  2: ('C', 0, False),
  3: ('C', 0, True),
  4: ('F', 0, None),
  5: ('F', 0, None),
  6: ('F', 0, None),
  7: ('O', 0, None),
  8: ('O', 0, None),
  9: ('H', 0, None),
  10: ('H', 0, None),
  11: ('H', 0, None),
  12: ('H', 0, None),
  13: ('H', 0, None),
  14: ('O', 0, None),
  15: ('H', 0, None)},
 {frozenset({7, 12}): (1, None),
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
  frozenset({0, 9}): (1, None)})
C4H5F3O2_GEO = (
 ('C', (-1.141471, 3.060908, -1.372282)),
 ('C', (-0.335935, 0.809777, -0.552249)),
 ('C', (-2.015423, -0.956788, 0.900211)),
 ('C', (2.245586, -0.127934, -1.109997)),
 ('F', (0.402603, 4.723624, -2.508427)),
 ('F', (-2.955248, -2.766816, -0.553588)),
 ('F', (2.089015, -2.665134, -1.889925)),
 ('O', (-0.608473, -1.927482, 3.005918)),
 ('O', (3.834623, -0.044323, 0.926263)),
 ('H', (-3.009159, 3.830615, -1.183382)),
 ('H', (-3.595336, 0.12939, 1.788313)),
 ('H', (3.088668, 0.820776, -2.776933)),
 ('H', (-1.019998, -3.680692, 3.265653)),
 ('H', (3.061224, -1.14277, 2.195312)),
 ('O', (-1.286205, 1.482841, -4.596348)),
 ('H', (-3.168531, 0.849827, -3.99654)))

# Fleeting Atom Stereo + (Reactant Atom Stereo => Product Bond Stereo)
# CC[C@H](O[O])C => C/C=C/C + O[O]
#  *
# [* marks a fleeting stereo site]
# An interesting case -- fleeting TS atom stereochemistry becomes bond
# stereochemistry
C4H9O2_TSG = (
 {0: ('C', 0, None),
  1: ('C', 0, None),
  2: ('C', 0, True),
  3: ('C', 0, True),
  4: ('O', 0, None),
  5: ('O', 0, None),
  6: ('H', 0, None),
  7: ('H', 0, None),
  8: ('H', 0, None),
  9: ('H', 0, None),
  10: ('H', 0, None),
  11: ('H', 0, None),
  12: ('H', 0, None),
  13: ('H', 0, None),
  14: ('H', 0, None)},
 {frozenset({2, 13}): (1, None),
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
  frozenset({0, 8}): (1, None)})
C4H9O2_GEO = (
 ('C', (2.840842, -0.918093, 0.301416)),
 ('C', (-2.107066, -1.899196, -0.110792)),
 ('C', (0.987221, 1.052291, 1.070762)),
 ('C', (-1.279309, 0.781779, -0.668652)),
 ('O', (-1.659829, 4.677101, -0.256324)),
 ('O', (-2.952545, 2.503599, 0.312223)),
 ('H', (2.360173, -1.778041, -1.574524)),
 ('H', (4.773451, -0.246703, 0.187287)),
 ('H', (2.703926, -2.528999, 1.689233)),
 ('H', (-3.298375, -2.425681, -1.754879)),
 ('H', (-0.350923, -3.043722, 0.066643)),
 ('H', (-3.108548, -1.931614, 1.691764)),
 ('H', (1.635613, 2.957619, 0.475404)),
 ('H', (0.428192, 1.019683, 3.075372)),
 ('H', (-0.462632, 1.023199, -2.601261)))


def test__set_stereo_from_geometry():
    """ test graph.set_stereo_from_geometry
    """
    pass