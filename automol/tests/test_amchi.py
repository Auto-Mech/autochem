""" test automol.amchi
"""
from automol import amchi

AR_CHI = 'AMChI=1/Ar'
CH4O_CH_CHI = 'AMChI=1/CH4O.CH/c1-2;/h2H,1H3;h1H'
CH2O2_CHI = 'AMChI=1/CH2O2/c2-1-3/h1-2H/q+1/p+1'
C2H6O_CHI = 'AMChI=1/C2H6O/c1-2-3/h3H,2H2,1H3/i2D/t2-/m1/s1'
C2H6O_CHI_NO_STEREO = 'AMChI=1/C2H6O/c1-2-3/h3H,2H2,1H3/i2D'

C4H10ZN_CHI = 'AMChI=1/2C2H5.Zn/c2*1-2;/h2*1H2,2H3;'
C4H10ZN_CHIS = (
    'AMChI=1/C2H5/c1-2/h1H2,2H3', 'AMChI=1/C2H5/c1-2/h1H2,2H3',
    'AMChI=1/Zn')

C4H5F2O_CHI = 'AMChI=1/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H/b3-1+,4-2+;'
C4H5F2O_CHI_NO_STEREO = 'AMChI=1/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H'
C4H5F2O_CHIS = (
    'AMChI=1/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1+,4-2+', 'AMChI=1/HO/h1H')

C2H2F2_CHI = 'AMChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1+'
C2H2F2_CHI_NO_STEREO = 'AMChI=1/C2H2F2/c3-1-2-4/h1-2H'
C2H2F2_CHI_STEREO_UNKNOWN = 'InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1?'

C2H4F2O2_CHI = 'AMChI=1/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2-/m0/s1'

C8H13O_CHI = (
    'AMChI=1/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-/m0/s1')
C8H13O_CHI_NO_STEREO = 'AMChI=1/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3'

C3H8FNO2_CHI = 'InChI=1S/C3H8FNO2/c1-3(4,2-5)7-6/h6H,2,5H2,1H3'
C10H14CLFO_CHI = (
    'AMChI=1/C10H14ClFO/c1-7(9(6-12)10(13)5-11)8-3-2-4-8'
    '/h2-4,7,9-10,13H,5-6H2,1H3')


def test__from_data():
    """ test getters
    """
    assert AR_CHI == amchi.from_data(
        fml_str=amchi.formula_string(AR_CHI),
    )

    assert CH2O2_CHI == amchi.from_data(
        fml_str=amchi.formula_string(CH2O2_CHI),
        main_lyr_dct=amchi.main_layers(CH2O2_CHI),
        char_lyr_dct=amchi.charge_layers(CH2O2_CHI),
    )

    assert C2H6O_CHI == amchi.from_data(
        fml_str=amchi.formula_string(C2H6O_CHI),
        main_lyr_dct=amchi.main_layers(C2H6O_CHI),
        iso_lyr_dct=amchi.isotope_layers(C2H6O_CHI),
    )

    assert C2H2F2_CHI == amchi.from_data(
        fml_str=amchi.formula_string(C2H2F2_CHI),
        main_lyr_dct=amchi.main_layers(C2H2F2_CHI),
        ste_lyr_dct=amchi.stereo_layers(C2H2F2_CHI),
    )

    assert C8H13O_CHI == amchi.from_data(
        fml_str=amchi.formula_string(C8H13O_CHI),
        main_lyr_dct=amchi.main_layers(C8H13O_CHI),
        ste_lyr_dct=amchi.stereo_layers(C8H13O_CHI),
    )


def test__version():
    """ amchi.version
    """
    assert amchi.version(C2H2F2_CHI) == '1'


def test__formula_string():
    """ amchi.formula_string
    """

    assert amchi.formula_string(AR_CHI) == 'Ar'
    assert amchi.formula_string(CH4O_CH_CHI) == 'CH4O.CH'
    assert amchi.formula_string(CH2O2_CHI) == 'CH2O2'
    assert amchi.formula_string(C2H6O_CHI) == 'C2H6O'


def test__formula():
    """ amchi.formula
    """
    assert amchi.formula(AR_CHI) == {'Ar': 1}
    assert amchi.formula(CH4O_CH_CHI) == {'C': 2, 'H': 5, 'O': 1}
    assert amchi.formula(CH2O2_CHI) == {'C': 1, 'H': 2, 'O': 2}
    assert amchi.formula(C2H6O_CHI) == {'C': 2, 'H': 6, 'O': 1}


def test__is_enantiomer():
    """ test amchi.is_enantiomer
    """
    chi1 = 'AMChI=1/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m1/s1'
    chi2 = 'AMChI=1/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3'
    assert amchi.is_enantiomer(chi1)
    assert not amchi.is_enantiomer(chi2)


def test__has_stereo():
    """ test amchi.has_stereo
    """
    assert amchi.has_stereo(C2H6O_CHI)
    assert amchi.has_stereo(C4H5F2O_CHI)
    assert amchi.has_stereo(C2H2F2_CHI)
    assert amchi.has_stereo(C8H13O_CHI)
    assert not amchi.has_stereo(AR_CHI)
    assert not amchi.has_stereo(C2H6O_CHI_NO_STEREO)


def test__symbols():
    """ test amchi.symbols
    """
    symb_dct = amchi.symbols(C3H8FNO2_CHI)
    print(symb_dct)
    assert symb_dct == {0: 'C', 1: 'C', 2: 'C', 3: 'F', 4: 'N', 5: 'O', 6: 'O'}

    symb_dct = amchi.symbols(C10H14CLFO_CHI)
    print(symb_dct)
    assert symb_dct == {0: 'C', 1: 'C', 2: 'C', 3: 'C', 4: 'C', 5: 'C', 6: 'C',
                        7: 'C', 8: 'C', 9: 'C', 10: 'Cl', 11: 'F', 12: 'O'}


def test__bonds():
    """ test amchi.bonds
    """
    bnds = amchi.bonds(C3H8FNO2_CHI, one_indexed=True)
    print(bnds)
    assert bnds == {
        frozenset({1, 3}), frozenset({3, 4}), frozenset({3, 2}),
        frozenset({2, 5}), frozenset({3, 7}), frozenset({7, 6})}

    bnds = amchi.bonds(C10H14CLFO_CHI, one_indexed=True)
    print(bnds)
    assert bnds == {
        frozenset({1, 7}), frozenset({7, 9}), frozenset({7, 8}),
        frozenset({8, 3}), frozenset({3, 2}), frozenset({2, 4}),
        frozenset({4, 8}), frozenset({9, 6}), frozenset({9, 10}),
        frozenset({6, 12}), frozenset({10, 13}), frozenset({10, 5}),
        frozenset({5, 11})}


def test__stereo_atoms():
    """ test amchi.stereo_atoms
    """
    atms = amchi.stereo_atoms(C2H6O_CHI)
    print(atms)
    assert atms == (1,)

    atms = amchi.stereo_atoms(C2H4F2O2_CHI)
    print(atms)
    assert atms == (0, 1)


def test__stereo_bonds():
    """ test amchi.stereo_bonds
    """
    bnds = amchi.stereo_bonds(C8H13O_CHI)
    print(bnds)
    assert bnds == ((4, 2), (5, 3))


def test__join():
    """ test amchi.join
    """
    assert amchi.join(amchi.split(C4H10ZN_CHI)) == C4H10ZN_CHI
    assert amchi.join(amchi.split(C4H5F2O_CHI)) == C4H5F2O_CHI

    chi = ('AMChI=1/C3H7O4.C2H5FO/c1-3(7-5)2-6-4;1-2(3)4/'
           'h3,5H,2H2,1H3;2,4H,1H3/t3-;2-/m01/s1')
    print(chi[:-3])
    print(amchi.join(amchi.split(chi)))
    assert amchi.join(amchi.split(chi)) == chi[:-3]


if __name__ == '__main__':
    test__symbols()
    test__bonds()
