""" test automol.inchi
"""
import numpy
from automol import inchi

AR_ICH = 'InChI=1S/Ar'
# CH_H2_ICH = 'InChI=1S/CH.H2/h1H;h1H'
CH4O_CH_ICH = 'InChI=1S/CH4O.CH/c1-2;/h2H,1H3;h1H'
CH2O2_ICH = 'InChI=1S/CH2O2/c2-1-3/h1-2H/q+1/p+1'
C2H6O_ICH = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i2D/t2-/m1/s1'
C2H6O_ICH_NO_STEREO = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i2D'

C4H10ZN_ICH = 'InChI=1S/2C2H5.Zn/c2*1-2;/h2*1H2,2H3;'
C4H10ZN_ICHS = (
    'InChI=1S/C2H5/c1-2/h1H2,2H3', 'InChI=1S/C2H5/c1-2/h1H2,2H3',
    'InChI=1S/Zn')

C4H5F2O_ICH = 'InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H/b3-1+,4-2+;'
C4H5F2O_ICH_NO_STEREO = 'InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H'
C4H5F2O_ICHS = (
    'InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1+,4-2+', 'InChI=1S/HO/h1H')

C2H2F2_ICH = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+'
C2H2F2_ICH_NO_STEREO = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H'
C2H2F2_ICH_STEREO_UNKNOWN = 'InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1?'

C2H4F2O2_ICH = 'InChI=1S/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2-/m0/s1'

C8H13O_ICH = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-/m0/s1')
C8H13O_ICH_NO_STEREO = 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3'


def test__from_data():
    """ test getters
    """
    assert AR_ICH == inchi.standard_form(inchi.from_data(
        fml_lyr=inchi.formula_layer(AR_ICH),
    ))

    assert CH2O2_ICH == inchi.standard_form(inchi.from_data(
        fml_lyr=inchi.formula_layer(CH2O2_ICH),
        main_lyr_dct=inchi.main_layers(CH2O2_ICH),
        char_lyr_dct=inchi.charge_layers(CH2O2_ICH),
    ))

    assert C2H6O_ICH == inchi.standard_form(inchi.from_data(
        fml_lyr=inchi.formula_layer(C2H6O_ICH),
        main_lyr_dct=inchi.main_layers(C2H6O_ICH),
        iso_lyr_dct=inchi.isotope_layers(C2H6O_ICH),
    ))

    assert C2H2F2_ICH == inchi.standard_form(inchi.from_data(
        fml_lyr=inchi.formula_layer(C2H2F2_ICH),
        main_lyr_dct=inchi.main_layers(C2H2F2_ICH),
        ste_lyr_dct=inchi.stereo_layers(C2H2F2_ICH),
    ))

    assert C8H13O_ICH == inchi.standard_form(inchi.from_data(
        fml_lyr=inchi.formula_layer(C8H13O_ICH),
        main_lyr_dct=inchi.main_layers(C8H13O_ICH),
        ste_lyr_dct=inchi.stereo_layers(C8H13O_ICH),
    ))


def test__formula_string():
    """ inchi.formula_string
    """

    assert inchi.formula_layer(AR_ICH) == 'Ar'
    assert inchi.formula_layer(CH4O_CH_ICH) == 'CH4O.CH'
    assert inchi.formula_layer(CH2O2_ICH) == 'CH2O2'
    assert inchi.formula_layer(C2H6O_ICH) == 'C2H6O'


def test__version():
    """ inchi.version
    """
    assert inchi.version(C2H2F2_ICH) == '1S'
    assert inchi.version(C2H2F2_ICH_NO_STEREO) == '1S'
    assert inchi.version(C2H2F2_ICH_STEREO_UNKNOWN) == '1'


def test__standard_form():
    """ test inchi.standard_form
    """
    assert (inchi.standard_form(C2H6O_ICH, stereo=False) ==
            C2H6O_ICH_NO_STEREO)
    assert (inchi.standard_form(C4H5F2O_ICH, stereo=False) ==
            C4H5F2O_ICH_NO_STEREO)
    assert (inchi.standard_form(C2H2F2_ICH, stereo=False) ==
            C2H2F2_ICH_NO_STEREO)
    assert (inchi.standard_form(C8H13O_ICH, stereo=False) ==
            C8H13O_ICH_NO_STEREO)


def test__has_stereo():
    """ test inchi.has_stereo
    """
    assert inchi.has_stereo(C2H6O_ICH)
    assert inchi.has_stereo(C4H5F2O_ICH)
    assert inchi.has_stereo(C2H2F2_ICH)
    assert inchi.has_stereo(C8H13O_ICH)
    assert not inchi.has_stereo(
        inchi.standard_form(C2H6O_ICH, stereo=False))
    assert not inchi.has_stereo(
        inchi.standard_form(C4H5F2O_ICH, stereo=False))
    assert not inchi.has_stereo(
        inchi.standard_form(C2H2F2_ICH, stereo=False))
    assert not inchi.has_stereo(
        inchi.standard_form(C8H13O_ICH, stereo=False))


def test__split():
    """ test inchi.split
    """
    assert (tuple(map(inchi.standard_form, inchi.split(C4H10ZN_ICH))) ==
            C4H10ZN_ICHS)
    assert (tuple(map(inchi.standard_form, inchi.split(C4H5F2O_ICH))) ==
            C4H5F2O_ICHS)

    ich = ('InChI=1S/C3H7O4.C2H5FO/c1-3(7-5)2-6-4;1-2(3)4/'
           'h3,5H,2H2,1H3;2,4H,1H3/t3-;2-/m01/s1')
    assert tuple(map(inchi.standard_form, inchi.split(ich))) == (
        'InChI=1S/C3H7O4/c1-3(7-5)2-6-4/h3,5H,2H2,1H3/t3-/m0/s1',
        'InChI=1S/C2H5FO/c1-2(3)4/h2,4H,1H3/t2-/m1/s1')


def test__join():
    """ test inchi.join
    """
    assert (inchi.standard_form(inchi.join(inchi.split(C4H10ZN_ICH))) ==
            C4H10ZN_ICH)
    assert (inchi.standard_form(inchi.join(inchi.split(C4H5F2O_ICH))) ==
            C4H5F2O_ICH)

    ich = ('InChI=1S/C3H7O4.C2H5FO/c1-3(7-5)2-6-4;1-2(3)4/'
           'h3,5H,2H2,1H3;2,4H,1H3/t3-;2-/m01/s1')
    assert inchi.standard_form(inchi.join(inchi.split(ich))) == ich


def test__argsort():
    """ test inchi.argsort
    """
    ref_ichs = ['InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3', 'InChI=1S/CH2/h1H2',
                'InChI=1S/BH3/h1H3', 'InChI=1S/H3N/h1H3', 'InChI=1S/H2O/h1H2']

    for _ in range(50):
        ichs = numpy.random.permutation(ref_ichs)
        idxs = inchi.argsort(ichs)
        srt_ichs = [ichs[idx] for idx in idxs]
        assert srt_ichs == ref_ichs


def test__recalculate():
    """ inchi.recalculate
    """
    assert inchi.recalculate(C2H2F2_ICH_NO_STEREO) == C2H2F2_ICH_NO_STEREO
    assert (inchi.recalculate(C2H2F2_ICH_NO_STEREO, stereo=True)
            == C2H2F2_ICH_STEREO_UNKNOWN)


def test__stereo_atoms():
    """ test inchi.stereo_atoms
    """
    atms = inchi.stereo_atoms(C2H6O_ICH)
    print(atms)
    assert atms == (1,)

    atms = inchi.stereo_atoms(C2H4F2O2_ICH)
    print(atms)
    assert atms == (0, 1)


def test__stereo_bonds():
    """ test inchi.stereo_bonds
    """
    bnds = inchi.stereo_bonds(C8H13O_ICH)
    print(bnds)
    assert bnds == ((4, 2), (5, 3))


def test__is_enantiomer():
    """ test inchi.is_enantiomer
    """
    ich1 = 'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m1/s1'
    ich2 = 'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3'
    assert inchi.is_enantiomer(ich1)
    assert not inchi.is_enantiomer(ich2)


def test__reflect():
    """ test inchi.reflect
    """
    ich1 = 'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m1/s1'
    ich2 = 'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3'
    assert (inchi.reflect(ich1) ==
            'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m0/s1')
    assert inchi.reflect(ich2) == ich2


def test__are_diastereomers():
    """ test inchi.diasteromer
    """

    assert not inchi.are_diastereomers(
        'InChI=1S/C4H8O3/c1-3(5)4(2)7-6/h4,6H,1-2H3/t4-/m0/s1',
        'InChI=1S/C4H8O3/c1-4(7-6)2-3-5/h3-4,6H,2H2,1H3/t4-/m0/s1'
    )
    assert not inchi.are_diastereomers(
        'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m1/s1',
        'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m1/s1'
    )
    assert not inchi.are_diastereomers(
        'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m1/s1',
        'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m0/s1'
    )
    assert inchi.are_diastereomers(
        'InChI=1S/C4H8O3/c1-3-4(7-3)2-6-5/h3-5H,2H2,1H3/t3-,4-/m0/s1',
        'InChI=1S/C4H8O3/c1-3-4(7-3)2-6-5/h3-5H,2H2,1H3/t3-,4+/m0/s1'
    )
    assert inchi.are_diastereomers(
        'InChI=1S/C4H8O2/c1-2-3-4-6-5/h3-5H,2H2,1H3/b4-3-',
        'InChI=1S/C4H8O2/c1-2-3-4-6-5/h3-5H,2H2,1H3/b4-3+'
    )
    assert inchi.are_diastereomers(
        'InChI=1S/C4H7O2/c1-4(6)2-3-5/h2-5H,1H3/b3-2-/t4-/m0/s1',
        'InChI=1S/C4H7O2/c1-4(6)2-3-5/h2-5H,1H3/b3-2+/t4-/m0/s1'
    )


def test__stereo():
    """ test inchi.add_stereo
        test inchi.expand_stereo
    """

    # Add and Expand stereo to C8H13O
    c8h13o_ste = (
        ('InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
         'b5-3+,6-4+/t8-/m1/s1'),
        ('InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
         'b5-3+,6-4-/t8-/m1/s1'),
        ('InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
         'b5-3-,6-4+/t8-/m1/s1'),
        ('InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
         'b5-3-,6-4-/t8-/m1/s1'),
        ('InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
         'b5-3+,6-4+/t8-/m0/s1'),
        ('InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
         'b5-3+,6-4-/t8-/m0/s1'),
        ('InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
         'b5-3-,6-4+/t8-/m0/s1'),
        ('InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
         'b5-3-,6-4-/t8-/m0/s1')
    )
    assert inchi.add_stereo(C8H13O_ICH_NO_STEREO) in c8h13o_ste
    assert set(inchi.expand_stereo(C8H13O_ICH_NO_STEREO)) == set(c8h13o_ste)

    # some cases that were breaking
    assert set(inchi.expand_stereo('InChI=1S/H2N2/c1-2/h1-2H')) == {
        'InChI=1S/H2N2/c1-2/h1-2H/b2-1+',
        'InChI=1S/H2N2/c1-2/h1-2H/b2-1-'}
    assert set(inchi.expand_stereo('InChI=1S/CH2N/c1-2/h1-2H')) == {
        'InChI=1S/CH2N/c1-2/h1-2H'}
    assert set(inchi.expand_stereo('InChI=1S/C2/c1-2')) == {
        'InChI=1S/C2/c1-2'}
    assert set(inchi.expand_stereo('InChI=1S/C3H3/c1-3-2/h1-3H')) == {
        'InChI=1S/C3H3/c1-3-2/h1-3H'}


def test__racemic():
    """ test amchi.racemic
    """
    chi = 'InChI=1S/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2-/m0/s1'
    assert (inchi.racemic(chi) ==
            'InChI=1/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2-/s3')

    assert inchi.racemic(AR_ICH) == AR_ICH
    assert inchi.racemic(C8H13O_ICH_NO_STEREO) == C8H13O_ICH_NO_STEREO


if __name__ == '__main__':
    # test__stereo_atoms()
    # test__stereo_bonds()
    # test__stereo()
    # test__are_diastereomers()
    # test__is_enantiomer()
    # test__reflect()
    # test__filter_enantiomer_reactions()
    # test__stereo()
    test__racemic()
