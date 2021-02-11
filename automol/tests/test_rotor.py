"""
    test automol.rotor
"""

import automol


# Species Z-Matrix
C2H5OH_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(
        automol.smiles.inchi('CCO')))
C3H7OH_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(
        automol.smiles.inchi('CCCO')))
C6H13OH_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(
        automol.smiles.inchi('CCC(C)CCO')))

# C2H5OH+CH3: H Abstraction
TS_ZMA = (
    ('C', (None, None, None), (None, None, None),
     (None, None, None)),
    ('C', (0, None, None), ('R1', None, None),
     (2.8527305589911376, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (2.0650927099003145, 1.9970598353387239, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
     (2.4921708142537837, 1.85108224203117, 4.282842507517861)),
    ('H', (0, 1, 2), ('R4', 'A4', 'D4'),
     (2.06395887422504, 1.9832682435894646, 2.283384589786898)),
    ('O', (1, 0, 2), ('R5', 'A5', 'D5'),
     (2.672072741397369, 1.964896907882972, 3.08662176346898)),
    ('H', (1, 0, 5), ('R6', 'A6', 'D6'),
     (2.0828561354796182, 1.9173227231321108, 4.109370742503641)),
    ('H', (1, 0, 5), ('R7', 'A7', 'D7'),
     (2.0724626417896, 1.9299676335628095, 2.054542254253157)),
    ('H', (5, 1, 0), ('R8', 'A8', 'D8'),
     (1.8233967384542584, 1.8720541183231338, 5.246932715722245)),
    ('X', (3, 0, 1), ('R9', 'A9', 'D9'),
     (1.8897261254578281, 1.5707963267948966, 3.141592653589793)),
    ('C', (3, 9, 0), ('R10', 'A10', 'D10'),
     (2.5978065046668766, 1.6376022945734834, 3.13944589860984)),
    ('H', (10, 3, 9), ('R11', 'A11', 'D11'),
     (2.0588566136863036, 1.8250157496526347, 5.830364868737413)),
    ('H', (10, 11, 3), ('R12', 'A12', 'D12'),
     (2.061313257649399, 1.9735851568994, 1.9984421361063034)),
    ('H', (10, 11, 12), ('R13', 'A13', 'D13'),
     (2.0605573671992157, 1.9730615581238018, 2.284154279987027)))


def test__tors():
    """ Various torsion builders
    """

    gra, lin_keys = automol.rotor.graph_with_keys(C2H5OH_ZMA)

    all_axes = automol.rotor.all_torsion_axes(gra, lin_keys)
    print(all_axes)
    assert all_axes == frozenset({frozenset({0, 1}), frozenset({1, 5})})

    all_groups = automol.rotor.all_torsion_groups(gra, lin_keys)
    assert all_groups == (
        ((2, 3, 4), (5, 6, 7, 8)),
        ((0, 2, 3, 4, 6, 7), (8,)))

    all_symms = automol.rotor.all_torsion_symmetries(gra, lin_keys)
    assert all_symms == (3, 1)

    axes1 = list(all_axes)[0]
    assert automol.rotor.torsion_groups(gra, axes1) == all_groups[0]
    assert automol.rotor.torsion_symmetry(gra, axes1, lin_keys) == all_symms[0]


def test__rotor():
    """ rotor
    """

    rotors = automol.rotor.from_zmatrix(C3H7OH_ZMA)

    for rotor in rotors:
        for tors in rotor:
            print(tors.span)
            print(tors.indices)

    rotor_names1 = automol.rotor.names(rotors)
    assert rotor_names1 == (('D5',), ('D8',), ('D11',))

    rotor_names2 = automol.rotor.names(rotors, flat=True)
    assert rotor_names2 == ('D5', 'D8', 'D11')

    rotor_axes = automol.rotor.axes(rotors)
    assert rotor_axes == (
        (frozenset({0, 1}),),
        (frozenset({1, 5}),),
        (frozenset({8, 5}),))

    rotor_groups = automol.rotor.groups(rotors)
    assert rotor_groups == (
        (((2, 3, 4), (5, 6, 7, 8, 9, 10, 11)),),
        (((0, 2, 3, 4, 6, 7), (8, 9, 10, 11)),),
        (((11,), (0, 1, 2, 3, 4, 6, 7, 9, 10)),))

    rotor_symms = automol.rotor.symmetries(rotors)
    assert rotor_symms == ((3,), (1,), (1,))


def test__ts():
    """ test a transition state build of rotor obj
    """

    # Replace with a zrxn object or just steal from test_reac

    rotors = automol.rotor.from_zmatrix(TS_ZMA)

    rotor_names1 = automol.rotor.names(rotors)
    assert rotor_names1 == (('D5',), ('D8',), ('D10',))

    rotor_names2 = automol.rotor.names(rotors, flat=True)
    assert rotor_names2 == ('D5', 'D8', 'D10')

    rotor_axes = automol.rotor.axes(rotors)
    assert rotor_axes == (
        (frozenset({0, 1}),),
        (frozenset({1, 5}),),
        (frozenset({9, 3}),))

    rotor_groups = automol.rotor.groups(rotors)
    assert rotor_groups == (
        (((2, 3, 4, 9, 10, 11, 12), (5, 6, 7, 8)),),
        (((0, 2, 3, 4, 6, 7, 9, 10, 11, 12), (8,)),),
        (((10, 11, 12), (0, 1, 2, 4, 5, 6, 7, 8)),))

    rotor_symms = automol.rotor.symmetries(rotors)
    assert rotor_symms == ((1,), (1,), (3,))


def test__string():
    """ test.automol.rotor.tors.string
        test.automol.rotor.tors.from_string
    """

    rotors = automol.rotor.from_zmatrix(C3H7OH_ZMA)
    tors_str = automol.rotor.string(rotors)
    print(tors_str)
    tors_dct = automol.rotor.from_string(tors_str)
    print(tors_dct)

    # assert rotors == automol.rotor.from_data(C3H7OH_ZMA, tors_dct)


def test__name_input():
    """ test various ways to do the build
    """

    inp_names1 = (('D5',), ('D11',), ('D8', 'D14', 'D17', 'D20'))
    rotors1 = automol.rotor.from_zmatrix(C6H13OH_ZMA, tors_names=inp_names1)
    assert automol.rotor.names(rotors1) == (
        ('D5',), ('D11',), ('D8', 'D14', 'D17', 'D20'))
    assert automol.rotor.axes(rotors1) == (
        (frozenset({0, 1}),),
        (frozenset({8, 5}),),
        (frozenset({1, 5}), frozenset({9, 5}),
         frozenset({9, 14}), frozenset({17, 14})))
    assert automol.rotor.symmetries(rotors1) == (
        (3,), (3,), (1, 1, 1, 1))

    inp_names2 = (('D5',), ('D11',))
    rotors2 = automol.rotor.from_zmatrix(C6H13OH_ZMA, tors_names=inp_names2)
    assert automol.rotor.names(rotors2) == (
        ('D5',), ('D11',))
    assert automol.rotor.axes(rotors2) == (
        (frozenset({0, 1}),), (frozenset({8, 5}),))
    assert automol.rotor.symmetries(rotors2) == (
        (3,), (3,))

    inp_names3 = (('D8', 'D14', 'D17', 'D20'),)
    rotors3 = automol.rotor.from_zmatrix(C6H13OH_ZMA, tors_names=inp_names3)
    assert automol.rotor.names(rotors3) == (
        ('D8', 'D14', 'D17', 'D20'),)
    assert automol.rotor.axes(rotors3) == (
        (frozenset({1, 5}), frozenset({9, 5}),
         frozenset({9, 14}), frozenset({17, 14})),)
    assert automol.rotor.symmetries(rotors3) == (
        (1, 1, 1, 1),)


def test__mdhr():
    """ building mdhr with the sorting
    """

    rotors = automol.rotor.from_zmatrix(C6H13OH_ZMA, multi=True)
    assert automol.rotor.names(rotors) == (
        ('D8', 'D14', 'D17', 'D20'), ('D5',), ('D11',))
    assert automol.rotor.axes(rotors) == (
        (frozenset({1, 5}), frozenset({9, 5}),
         frozenset({9, 14}), frozenset({17, 14})),
        (frozenset({0, 1}),),
        (frozenset({8, 5}),))
    assert automol.rotor.groups(rotors) == (
        (((0, 2, 3, 4, 6, 7),
          (8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
         ((14, 15, 16, 17, 18, 19, 20),
          (0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13)),
         ((0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16),
          (17, 18, 19, 20)),
         ((20,),
          (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 18, 19))),
        (((2, 3, 4),
          (5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),),
        (((11, 12, 13),
          (0, 1, 2, 3, 4, 6, 7, 9, 10, 14, 15, 16, 17, 18, 19, 20)),))
    assert automol.rotor.symmetries(rotors) == (
        (1, 1, 1, 1), (3,), (3,))

    # test break up when multi is not given but names are


def test__relabel():
    """
    """

    rotors = automol.rotor.from_zmatrix(TS_ZMA)
    geo, rotors = automol.rotor.relabel_for_geometry(rotors)

    print(automol.zmat.string(TS_ZMA, one_indexed=False))
    print()
    # print(automol.geom.string(geo))
    # print()

    rotor_axes = automol.rotor.axes(rotors)
    print(rotor_axes)
    # assert rotor_axes == (
    #     (frozenset({0, 1}),),
    #     (frozenset({1, 5}),),
    #     (frozenset({9, 3}),))

    rotor_groups = automol.rotor.groups(rotors)
    print(rotor_groups)
    # assert rotor_groups == (
    #     (((2, 3, 4, 9, 10, 11, 12), (5, 6, 7, 8)),),
    #     (((0, 2, 3, 4, 6, 7, 9, 10, 11, 12), (8,)),),
    #     (((10, 11, 12), (0, 1, 2, 4, 5, 6, 7, 8)),))

    rotor_symms = automol.rotor.symmetries(rotors)
    print(rotor_symms)
    # assert rotor_symms == ((1,), (1,), (3,))


if __name__ == '__main__':
    test__tors()
    test__rotor()
    # test__ts()
    # test__string()
    # test__name_input()
    # test__mdhr()
    # test__relabel()
