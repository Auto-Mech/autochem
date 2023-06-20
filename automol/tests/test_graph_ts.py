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

# Fleeting Atom Stereo (Symmetric Reaction)
# CC(O)C + [CH3] => CC(O)C + [CH3]
#  *
# [* marks a fleeting TS stereosite]
C4H11O_TSG = (
    {0: ('C', 0, None),
     1: ('C', 0, None),
     2: ('C', 0, True),
     3: ('O', 0, None),
     4: ('H', 0, None),
     5: ('H', 0, None),
     6: ('H', 0, None),
     7: ('H', 0, None),
     8: ('H', 0, None),
     9: ('H', 0, None),
     10: ('H', 0, None),
     11: ('H', 0, None),
     12: ('C', 0, None),
     13: ('H', 0, None),
     14: ('H', 0, None),
     15: ('H', 0, None)},
    {frozenset({1, 9}): (1, None),
     frozenset({1, 7}): (1, None),
     frozenset({0, 6}): (1, None),
     frozenset({2, 3}): (1, None),
     frozenset({1, 2}): (0.9, None),
     frozenset({2, 10}): (1, None),
     frozenset({12, 14}): (1, None),
     frozenset({0, 2}): (1, None),
     frozenset({12, 15}): (1, None),
     frozenset({0, 4}): (1, None),
     frozenset({0, 5}): (1, None),
     frozenset({2, 12}): (0.1, None),
     frozenset({1, 8}): (1, None),
     frozenset({12, 13}): (1, None),
     frozenset({3, 11}): (1, None)})

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

# Fleeting Atom Stereo + Conserved Atom Stereo
# CCO[C@H](O[O])C => C[CH]O[C@H](OO)C
#  *
# [* marks a fleeting stereosite]
C4H9O3_TSG = (
    {0: ('C', 0, None),
     1: ('C', 0, None),
     2: ('C', 0, True),
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

# Conserved Atom Stereo + (Reactant Bond Stereo => Product Atom Stereo) 
# F/C=C([C@@H](F)O)\[C@H](F)O + [OH] => F[C@H]([C]([C@@H](F)O)[C@H](F)O)O
C4H5F3O2_TSG = (
    {0: ('C', 0, None, False, None), 1: ('C', 0, None, None, None),
     2: ('C', 0, False, False, None), 3: ('C', 0, True, True, None),
     4: ('F', 0, None, None, None), 5: ('F', 0, None, None, None),
     6: ('F', 0, None, None, None), 7: ('O', 0, None, None, None),
     8: ('O', 0, None, None, None), 9: ('H', 0, None, None, None),
     10: ('H', 0, None, None, None), 11: ('H', 0, None, None, None),
     12: ('H', 0, None, None, None), 13: ('H', 0, None, None, None),
     14: ('O', 0, None, None, None), 15: ('H', 0, None, None, None)},
    {frozenset({7, 12}): (1, None, None, None),
     frozenset({2, 10}): (1, None, None, None),
     frozenset({1, 2}): (1, None, None, None),
     frozenset({0, 1}): (1, False, None, None),
     frozenset({3, 6}): (1, None, None, None),
     frozenset({2, 5}): (1, None, None, None),
     frozenset({0, 4}): (1, None, None, None),
     frozenset({3, 8}): (1, None, None, None),
     frozenset({0, 14}): (0.1, None, None, None),
     frozenset({1, 3}): (1, None, None, None),
     frozenset({8, 13}): (1, None, None, None),
     frozenset({14, 15}): (1, None, None, None),
     frozenset({3, 11}): (1, None, None, None),
     frozenset({2, 7}): (1, None, None, None),
     frozenset({0, 9}): (1, None, None, None)})

C4H5F3O2_REV_TSG = (
    {0: ('C', 0, False, None, None),
     1: ('C', 0, None, None, None),
     2: ('C', 0, False, False, None),
     3: ('C', 0, True, True, None),
     4: ('F', 0, None, None, None),
     5: ('F', 0, None, None, None),
     6: ('F', 0, None, None, None),
     7: ('O', 0, None, None, None),
     8: ('O', 0, None, None, None),
     9: ('H', 0, None, None, None),
     10: ('H', 0, None, None, None),
     11: ('H', 0, None, None, None),
     12: ('H', 0, None, None, None),
     13: ('H', 0, None, None, None),
     14: ('O', 0, None, None, None),
     15: ('H', 0, None, None, None)},
    {frozenset({7, 12}): (1, None, None, None),
     frozenset({2, 10}): (1, None, None, None),
     frozenset({1, 2}): (1, None, None, None),
     frozenset({0, 1}): (1, None, False, None),
     frozenset({3, 6}): (1, None, None, None),
     frozenset({2, 5}): (1, None, None, None),
     frozenset({0, 4}): (1, None, None, None),
     frozenset({3, 8}): (1, None, None, None),
     frozenset({0, 14}): (0.9, None, None, None),
     frozenset({1, 3}): (1, None, None, None),
     frozenset({8, 13}): (1, None, None, None),
     frozenset({14, 15}): (1, None, None, None),
     frozenset({3, 11}): (1, None, None, None),
     frozenset({2, 7}): (1, None, None, None),
     frozenset({0, 9}): (1, None, None, None)})

# Fleeting Atom Stereo + (Reactant Atom Stereo => Product Bond Stereo)
# CC[C@H](O[O])C => C/C=C/C + O[O]
#  *
# [* marks a fleeting stereo site]
# An interesting case -- fleeting TS atom stereochemistry becomes bond
# stereochemistry
C4H9O2_TSG = (
    {0: ('C', 0, None, None, None), 1: ('C', 0, None, None, None),
     2: ('C', 0, None, None, None), 3: ('C', 0, True, None, None),
     4: ('O', 0, None, None, None), 5: ('O', 0, None, None, None),
     6: ('H', 0, None, None, None), 7: ('H', 0, None, None, None),
     8: ('H', 0, None, None, None), 9: ('H', 0, None, None, None),
     10: ('H', 0, None, None, None), 11: ('H', 0, None, None, None),
     12: ('H', 0, None, None, None), 13: ('H', 0, None, None, None),
     14: ('H', 0, None, None, None)},
    {frozenset({2, 13}): (1, None, None, None),
     frozenset({1, 9}): (1, None, None, None),
     frozenset({0, 6}): (1, None, None, None),
     frozenset({2, 3}): (1, None, True, None),
     frozenset({1, 11}): (1, None, None, None),
     frozenset({4, 5}): (1, None, None, None),
     frozenset({0, 2}): (1, None, None, None),
     frozenset({4, 12}): (0.1, None, None, None),
     frozenset({2, 12}): (0.9, None, None, None),
     frozenset({3, 14}): (1, None, None, None),
     frozenset({3, 5}): (0.9, None, None, None),
     frozenset({1, 3}): (1, None, None, None),
     frozenset({0, 7}): (1, None, None, None),
     frozenset({1, 10}): (1, None, None, None),
     frozenset({0, 8}): (1, None, None, None)})

C4H9O2_REV_TSG = (
    {0: ('C', 0, None, None, None), 1: ('C', 0, None, None, None),
     2: ('C', 0, None, None, None), 3: ('C', 0, None, True, None),
     4: ('O', 0, None, None, None), 5: ('O', 0, None, None, None),
     6: ('H', 0, None, None, None), 7: ('H', 0, None, None, None),
     8: ('H', 0, None, None, None), 9: ('H', 0, None, None, None),
     10: ('H', 0, None, None, None), 11: ('H', 0, None, None, None),
     12: ('H', 0, None, None, None), 13: ('H', 0, None, None, None),
     14: ('H', 0, None, None, None)},
    {frozenset({2, 13}): (1, None, None, None),
     frozenset({1, 9}): (1, None, None, None),
     frozenset({0, 6}): (1, None, None, None),
     frozenset({2, 3}): (1, True, None, None),
     frozenset({1, 11}): (1, None, None, None),
     frozenset({4, 5}): (1, None, None, None),
     frozenset({0, 2}): (1, None, None, None),
     frozenset({4, 12}): (0.9, None, None, None),
     frozenset({2, 12}): (0.1, None, None, None),
     frozenset({3, 14}): (1, None, None, None),
     frozenset({3, 5}): (0.1, None, None, None),
     frozenset({1, 3}): (1, None, None, None),
     frozenset({0, 7}): (1, None, None, None),
     frozenset({1, 10}): (1, None, None, None),
     frozenset({0, 8}): (1, None, None, None)})


# def test__ts__expand_reaction_stereo():
#     """ test graph.ts.stereo_expand_reverse_graphs
#     """
#     gra = C4H5F2O_TSG
#     assert len(graph.ts.expand_reaction_stereo(gra, enant=True)) == 16

#     gra = C4H5F3O2_TSG
#     assert len(graph.ts.expand_reaction_stereo(gra, enant=True)) == 8

#     # CC(OO)C(O[O])C(OO)C => CC(OO)C(OO)C(OO)[CH2]
#     gra = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
#             3: ('C', 0, None), 4: ('C', 0, None), 5: ('O', 0, None),
#             6: ('O', 0, None), 7: ('O', 0, None), 8: ('O', 0, None),
#             9: ('O', 0, None), 10: ('O', 0, None), 11: ('H', 0, None),
#             12: ('H', 0, None), 13: ('H', 0, None), 14: ('H', 0, None),
#             15: ('H', 0, None), 16: ('H', 0, None), 17: ('H', 0, None),
#             18: ('H', 0, None), 19: ('H', 0, None), 20: ('H', 0, None),
#             21: ('H', 0, None)},
#            {frozenset({10, 4}): (1, None), frozenset({8, 2}): (1, None),
#             frozenset({9, 3}): (1, None), frozenset({1, 15}): (1, None),
#             frozenset({0, 12}): (1, None), frozenset({18, 3}): (1, None),
#             frozenset({0, 11}): (0.9, None), frozenset({16, 1}): (1, None),
#             frozenset({11, 7}): (0.1, None), frozenset({0, 13}): (1, None),
#             frozenset({1, 14}): (1, None), frozenset({3, 4}): (1, None),
#             frozenset({9, 6}): (1, None), frozenset({21, 6}): (1, None),
#             frozenset({10, 7}): (1, None), frozenset({19, 4}): (1, None),
#             frozenset({0, 2}): (1, None), frozenset({17, 2}): (1, None),
#             frozenset({2, 4}): (1, None), frozenset({8, 5}): (1, None),
#             frozenset({20, 5}): (1, None), frozenset({1, 3}): (1, None)})
#     assert len(graph.ts.expand_reaction_stereo(gra, enant=False)) == 4
