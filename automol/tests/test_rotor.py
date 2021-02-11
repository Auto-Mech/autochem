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
C4H5OH_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(
        automol.smiles.inchi('C#CCCO')))
C6H13OH_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(
        automol.smiles.inchi('CCC(C)CCO')))


def test__tors():
    """ Various torsion builders
    """

    gra, lin_keys = automol.rotor.graph_with_keys(C2H5OH_ZMA)

    all_axes = automol.rotor.all_torsion_axes(gra, lin_keys)
    assert all_axes == ((0, 1), (1, 5))

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

    rotor_names1 = automol.rotor.names(rotors)
    rotor_names2 = automol.rotor.names(rotors, flat=True)
    assert rotor_names1 == (('D5',), ('D8',), ('D11',))
    assert rotor_names2 == ('D5', 'D8', 'D11')

    rotor_axes1 = automol.rotor.axes(rotors)
    rotor_axes2 = automol.rotor.axes(rotors, flat=True)
    assert rotor_axes1 == (((0, 1),), ((1, 5),), ((8, 5),))
    assert rotor_axes2 == ((0, 1), (1, 5), (8, 5))

    rotor_groups1 = automol.rotor.groups(rotors)
    rotor_groups2 = automol.rotor.groups(rotors, flat=True)
    assert rotor_groups1 == (
        (((2, 3, 4), (5, 6, 7, 8, 9, 10, 11)),),
        (((0, 2, 3, 4, 6, 7), (8, 9, 10, 11)),),
        (((11,), (0, 1, 2, 3, 4, 6, 7, 9, 10)),))
    assert rotor_groups2 == (
        ((2, 3, 4), (5, 6, 7, 8, 9, 10, 11)),
        ((0, 2, 3, 4, 6, 7), (8, 9, 10, 11)),
        ((11,), (0, 1, 2, 3, 4, 6, 7, 9, 10)))

    rotor_symms1 = automol.rotor.symmetries(rotors)
    rotor_symms2 = automol.rotor.symmetries(rotors, flat=True)
    assert rotor_symms1 == ((3,), (1,), (1,))
    assert rotor_symms2 == (3, 1, 1)

    rotor_grids1 = automol.rotor.grids(rotors)
    rotor_grids2 = automol.rotor.grids(rotors, flat=True)
    assert rotor_grids1 == (
        ((1.0274898979189644, 1.5510886735172633,
          2.074687449115562, 2.5982862247138607),),
        ((1.0913723623432035, 1.6149711379415024,
          2.138569913539801, 2.6621686891381,
          3.185767464736399, 3.709366240334697,
          4.232965015932996, 4.7565637915312955,
          5.280162567129594, 5.803761342727893,
          6.327360118326191, 6.85095889392449),),
        ((1.0207703162606525, 1.5443690918589512,
          2.0679678674572504, 2.591566643055549,
          3.1151654186538478, 3.6387641942521465,
          4.162362969850445, 4.685961745448744,
          5.2095605210470435, 5.733159296645342,
          6.25675807224364, 6.78035684784194),))
    assert rotor_grids2 == (
        (1.0274898979189644, 1.5510886735172633,
         2.074687449115562, 2.5982862247138607),
        (1.0913723623432035, 1.6149711379415024,
         2.138569913539801, 2.6621686891381,
         3.185767464736399, 3.709366240334697,
         4.232965015932996, 4.7565637915312955,
         5.280162567129594, 5.803761342727893,
         6.327360118326191, 6.85095889392449),
        (1.0207703162606525, 1.5443690918589512,
         2.0679678674572504, 2.591566643055549,
         3.1151654186538478, 3.6387641942521465,
         4.162362969850445, 4.685961745448744,
         5.2095605210470435, 5.733159296645342,
         6.25675807224364, 6.78035684784194))

    # Test the  geometry labeling
    _, grotors = automol.rotor.relabel_for_geometry(rotors)

    grotor_axes = automol.rotor.axes(grotors)
    assert grotor_axes == (((0, 1),), ((1, 5),), ((8, 5),))

    grotor_groups = automol.rotor.groups(grotors)
    assert grotor_groups == (
        (((2, 3, 4), (5, 6, 7, 8, 9, 10, 11)),),
        (((0, 2, 3, 4, 6, 7), (8, 9, 10, 11)),),
        (((11,), (0, 1, 2, 3, 4, 6, 7, 9, 10)),))


def test__rotor_wdummy():
    """ rotor
    """

    # Build with ZMA
    print(automol.zmat.string(C4H5OH_ZMA))
    rotors = automol.rotor.from_zmatrix(C4H5OH_ZMA)

    rotor_names = automol.rotor.names(rotors)
    print(rotor_names)
    rotor_axes = automol.rotor.axes(rotors)
    print(rotor_axes)

    rotor_groups = automol.rotor.groups(rotors)
    print(rotor_groups)

    rotor_symms = automol.rotor.symmetries(rotors)
    print(rotor_symms)

    # Test the  geometry labeling
    _, grotors = automol.rotor.relabel_for_geometry(rotors)

    grotor_axes = automol.rotor.axes(grotors)
    grotor_groups = automol.rotor.groups(grotors)
    print(grotor_axes)
    print(grotor_groups)


def test__name_input():
    """ test various ways to do the build
    """

    inp_names1 = (('D5',), ('D11',), ('D8', 'D14', 'D17', 'D20'))
    rotors1 = automol.rotor.from_zmatrix(C6H13OH_ZMA, tors_names=inp_names1)
    assert automol.rotor.names(rotors1) == (
        ('D5',), ('D11',), ('D8', 'D14', 'D17', 'D20'))
    assert automol.rotor.axes(rotors1) == (
        ((0, 1),),
        ((8, 5),),
        ((1, 5), (9, 5), (9, 14), (17, 14)))
    assert automol.rotor.symmetries(rotors1) == (
        (3,), (3,), (1, 1, 1, 1))

    inp_names2 = (('D5',), ('D11',))
    rotors2 = automol.rotor.from_zmatrix(C6H13OH_ZMA, tors_names=inp_names2)
    assert automol.rotor.names(rotors2) == (('D5',), ('D11',))
    assert automol.rotor.axes(rotors2) == (((0, 1),), ((8, 5),))
    assert automol.rotor.symmetries(rotors2) == ((3,), (3,))

    inp_names3 = (('D8', 'D14', 'D17', 'D20'),)
    rotors3 = automol.rotor.from_zmatrix(C6H13OH_ZMA, tors_names=inp_names3)
    assert automol.rotor.names(rotors3) == (('D8', 'D14', 'D17', 'D20'),)
    assert automol.rotor.axes(rotors3) == (
        ((1, 5), (9, 5), (9, 14), (17, 14)),)
    assert automol.rotor.symmetries(rotors3) == ((1, 1, 1, 1),)


def test__mdhr():
    """ building mdhr with the sorting
    """

    rotors = automol.rotor.from_zmatrix(C6H13OH_ZMA, multi=True)
    assert automol.rotor.names(rotors) == (
        ('D8', 'D14', 'D17', 'D20'), ('D5',), ('D11',))
    assert automol.rotor.axes(rotors) == (
        ((1, 5), (9, 5), (9, 14), (17, 14)), ((0, 1),), ((8, 5),))
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
    assert automol.rotor.dimensions(rotors) == (4, 1, 1)


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


if __name__ == '__main__':
    # test__tors()
    # test__rotor()
    test__rotor_wdummy()
    # test__name_input()
    # test__mdhr()
    # test__string()
