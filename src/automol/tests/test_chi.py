""" test automol.chi
"""
import numpy
from automol import chi

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
    assert AR_ICH == chi.standard_form(chi.from_data(
        fml_slyr=chi.formula_sublayer(AR_ICH),
    ))

    assert CH2O2_ICH == chi.standard_form(chi.from_data(
        fml_slyr=chi.formula_sublayer(CH2O2_ICH),
        main_lyr_dct=chi.main_sublayers(CH2O2_ICH),
        char_lyr_dct=chi.charge_sublayers(CH2O2_ICH),
    ))

    assert C2H6O_ICH == chi.standard_form(chi.from_data(
        fml_slyr=chi.formula_sublayer(C2H6O_ICH),
        main_lyr_dct=chi.main_sublayers(C2H6O_ICH),
        iso_lyr_dct=chi.isotope_sublayers(C2H6O_ICH),
    ))

    assert C2H2F2_ICH == chi.standard_form(chi.from_data(
        fml_slyr=chi.formula_sublayer(C2H2F2_ICH),
        main_lyr_dct=chi.main_sublayers(C2H2F2_ICH),
        ste_lyr_dct=chi.stereo_sublayers(C2H2F2_ICH),
    ))

    assert C8H13O_ICH == chi.standard_form(chi.from_data(
        fml_slyr=chi.formula_sublayer(C8H13O_ICH),
        main_lyr_dct=chi.main_sublayers(C8H13O_ICH),
        ste_lyr_dct=chi.stereo_sublayers(C8H13O_ICH),
    ))


def test__formula_string():
    """ chi.formula_string
    """

    assert chi.formula_layer(AR_ICH) == 'Ar'
    assert chi.formula_layer(CH4O_CH_ICH) == 'CH4O.CH'
    assert chi.formula_layer(CH2O2_ICH) == 'CH2O2'
    assert chi.formula_layer(C2H6O_ICH) == 'C2H6O'


def test__version():
    """ chi.version
    """
    assert chi.version(C2H2F2_ICH) == '1S'
    assert chi.version(C2H2F2_ICH_NO_STEREO) == '1S'
    assert chi.version(C2H2F2_ICH_STEREO_UNKNOWN) == '1'


def test__standard_form():
    """ test chi.standard_form
    """
    print(chi.standard_form(C2H6O_ICH, stereo=False))
    print(C2H6O_ICH_NO_STEREO)
    assert (chi.standard_form(C2H6O_ICH, stereo=False) ==
            C2H6O_ICH_NO_STEREO)
    assert (chi.standard_form(C4H5F2O_ICH, stereo=False) ==
            C4H5F2O_ICH_NO_STEREO)
    assert (chi.standard_form(C2H2F2_ICH, stereo=False) ==
            C2H2F2_ICH_NO_STEREO)
    assert (chi.standard_form(C8H13O_ICH, stereo=False) ==
            C8H13O_ICH_NO_STEREO)


def test__is_enantiomer():
    """ test chi.is_enantiomer
    """
    ich1 = 'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m1/s1'
    ich2 = 'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3'
    assert chi.is_enantiomer(ich1)
    assert not chi.is_enantiomer(ich2)


def test__reflect():
    """ test chi.reflect
    """
    ich1 = 'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m1/s1'
    ich2 = 'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3'
    print(chi.reflect(ich1))
    print('InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m0/s1')
    assert (chi.reflect(ich1) ==
            'InChI=1S/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m0/s1')
    assert chi.reflect(ich2) == ich2


def test__inchi_key():
    """ test chi.inchi_key
    """
    ich = 'InChI=1S/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    ach = 'AMChI=1/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    ick = chi.inchi_key(ich)
    ack = chi.inchi_key(ach)
    print(ick)
    print(ack)

    # AMChI leads to a non-standard InChI key
    assert ick[-4:] == 'SA-N'
    assert ack[-4:] == 'NA-N'
    # Otherwise, they are identical
    assert ick[:-4] == ack[:-4]


def test__formula():
    """ test chi.formula
    """
    ich = 'InChI=1S/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    ach = 'AMChI=1/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    fml1 = chi.formula(ich)
    fml2 = chi.formula(ach)
    print(fml1)
    print(fml2)

    assert fml1 == fml2 == {'C': 4, 'H': 7, 'F': 1, 'O': 1}


def test__is_standard_form():
    """ test chi.is_standard_form
    """
    ich = 'InChI=1S/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    ach = 'AMChI=1/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    print(chi.is_standard_form(ich))
    print(chi.is_standard_form(ach))
    assert chi.is_standard_form(ich)
    assert chi.is_standard_form(ach)


def test__has_multiple_components():
    """ test chi.has_multiple_components
    """
    ich = 'InChI=1S/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    ach = 'AMChI=1/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    print(chi.has_multiple_components(ich))
    print(chi.has_multiple_components(ach))
    print(chi.has_multiple_components(CH4O_CH_ICH))
    print(chi.has_multiple_components(C4H10ZN_ICH))

    assert not chi.has_multiple_components(ich)
    assert not chi.has_multiple_components(ach)
    assert chi.has_multiple_components(CH4O_CH_ICH)
    assert chi.has_multiple_components(C4H10ZN_ICH)


def test__has_stereo():
    """ test chi.has_stereo
    """
    assert chi.has_stereo(C2H6O_ICH)
    assert chi.has_stereo(C4H5F2O_ICH)
    assert chi.has_stereo(C2H2F2_ICH)
    assert chi.has_stereo(C8H13O_ICH)
    assert not chi.has_stereo(
        chi.standard_form(C2H6O_ICH, stereo=False))
    assert not chi.has_stereo(
        chi.standard_form(C4H5F2O_ICH, stereo=False))
    assert not chi.has_stereo(
        chi.standard_form(C2H2F2_ICH, stereo=False))
    assert not chi.has_stereo(
        chi.standard_form(C8H13O_ICH, stereo=False))

    ach = 'AMChI=1/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    print(chi.has_stereo(ach))


def test__same_connectivity():
    """ test chi.same_connectivity
    """
    ich1 = 'InChI=1S/CH2/h1H2'
    ach1 = 'AMChI=1S/CH2/h1H2'
    ich2 = 'InChI=1S/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    ach2 = 'AMChI=1/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    ich3 = 'InChI=1S/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2-/t4-/m0/s1'
    ach3 = 'AMChI=1/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    print(chi.same_connectivity(ich1, ach1))
    print(chi.same_connectivity(ich2, ach2))
    print(chi.same_connectivity(ich3, ach3))

    assert chi.same_connectivity(ich1, ach1)
    assert chi.same_connectivity(ich2, ach2)
    assert chi.same_connectivity(ich3, ach3)


def test__equivalent():
    """ test chi.equivalent
    """
    ich1 = 'InChI=1S/CH2/h1H2'
    ach1 = 'AMChI=1S/CH2/h1H2'
    ich2 = 'InChI=1S/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    ach2 = 'AMChI=1/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    ich3 = 'InChI=1S/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2-/t4-/m0/s1'
    ach3 = 'AMChI=1/C4H7FO/c1-4(6)2-3-5/h2-4,6H,1H3/b3-2+/t4-/m1/s1'
    print(chi.equivalent(ich1, ach1))
    print(chi.equivalent(ich2, ach2))
    print(chi.equivalent(ich3, ach3))

    assert chi.equivalent(ich1, ach1)
    assert chi.equivalent(ich2, ach2)
    assert not chi.equivalent(ich3, ach3)


def test__split():
    """ test chi.split
    """
    assert (tuple(map(chi.standard_form, chi.split(C4H10ZN_ICH))) ==
            C4H10ZN_ICHS)
    assert (tuple(map(chi.standard_form, chi.split(C4H5F2O_ICH))) ==
            C4H5F2O_ICHS)

    ich = ('InChI=1S/C3H7O4.C2H5FO/c1-3(7-5)2-6-4;1-2(3)4/'
           'h3,5H,2H2,1H3;2,4H,1H3/t3-;2-/m01/s1')
    assert tuple(map(chi.standard_form, chi.split(ich))) == (
        'InChI=1S/C3H7O4/c1-3(7-5)2-6-4/h3,5H,2H2,1H3/t3-/m0/s1',
        'InChI=1S/C2H5FO/c1-2(3)4/h2,4H,1H3/t2-/m1/s1')


def test__join():
    """ test chi.join
    """
    assert (chi.standard_form(chi.join(chi.split(C4H10ZN_ICH))) ==
            C4H10ZN_ICH)
    assert (chi.standard_form(chi.join(chi.split(C4H5F2O_ICH))) ==
            C4H5F2O_ICH)

    ich = ('InChI=1S/C3H7O4.C2H5FO/c1-3(7-5)2-6-4;1-2(3)4/'
           'h3,5H,2H2,1H3;2,4H,1H3/t3-;2-/m01/s1')
    assert chi.standard_form(chi.join(chi.split(ich))) == ich


def test__sorted_():
    """ test chi.sorted_
    """
    ref_ichs = ('InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3', 'InChI=1S/CH2/h1H2',
                'InChI=1S/BH3/h1H3', 'InChI=1S/H3N/h1H3', 'InChI=1S/H2O/h1H2')

    for _ in range(10):
        ichs = numpy.random.permutation(ref_ichs)
        srt_ichs = chi.sorted_(ichs)
        assert srt_ichs == ref_ichs

    ref_achs = (
        'AMChI=1/C6H14/c1-4-5-6(2)3/h6H,4-5H2,1-3H3',
        'AMChI=1/C6H14/c1-4-6(3)5-2/h6H,4-5H2,1-3H3',
        'AMChI=1/C6H14/c1-3-5-6-4-2/h3-6H2,1-2H3',
        'AMChI=1/C5H11/c1-4-5(2)3/h4H2,1-3H3',
        'AMChI=1/C5H11/c1-4-5(2)3/h4-5H,1-3H3',
        'AMChI=1/C5H11/c1-4-5(2)3/h5H,2,4H2,1,3H3',
        'AMChI=1/C5H11/c1-4-5(2)3/h5H,1,4H2,2-3H3',
        'AMChI=1/C4H10/c1-4(2)3/h4H,1-3H3',
        'AMChI=1/C4H10/c1-3-4-2/h3-4H2,1-2H3',
        'AMChI=1/C3H8O2/c1-2-3-5-4/h4H,2-3H2,1H3',
        'AMChI=1/C3H8/c1-3-2/h3H2,1-2H3',
        'AMChI=1/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2+',
        'AMChI=1/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2-/m1',
        'AMChI=1/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2-/m1',
        'AMChI=1/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2-/m0',
        'AMChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1+',
        'AMChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1-',
        'AMChI=1/C2H6/c1-2/h1-2H3',
        'AMChI=1/CH4/h1H4',
        'AMChI=1/CH3/h1H3',
        'AMChI=1/BH3/h1H3',
        'AMChI=1/H3N/h1H3'
    )

    for _ in range(10):
        achs = numpy.random.permutation(ref_achs)
        print(achs)
        print()
        srt_achs = chi.sorted_(achs)
        assert srt_achs == ref_achs


def test__recalculate():
    """ chi.recalculate
    """
    assert chi.recalculate(C2H2F2_ICH_NO_STEREO) == C2H2F2_ICH_NO_STEREO
    assert (chi.recalculate(C2H2F2_ICH_NO_STEREO, stereo=True)
            == C2H2F2_ICH_STEREO_UNKNOWN)


def test__stereo():
    """ test chi.add_stereo
        test chi.expand_stereo
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
    assert chi.add_stereo(C8H13O_ICH_NO_STEREO) in c8h13o_ste
    print(set(chi.expand_stereo(C8H13O_ICH_NO_STEREO)))
    print(len(set(chi.expand_stereo(C8H13O_ICH_NO_STEREO))))
    print(len(c8h13o_ste))
    assert set(chi.expand_stereo(C8H13O_ICH_NO_STEREO)) == set(c8h13o_ste)

    # some cases that were breaking
    print(chi.expand_stereo('InChI=1S/H2N2/c1-2/h1-2H'))
    assert set(chi.expand_stereo('InChI=1S/H2N2/c1-2/h1-2H')) == {
        'InChI=1S/H2N2/c1-2/h1-2H/b2-1+',
        'InChI=1S/H2N2/c1-2/h1-2H/b2-1-'}
    print(chi.expand_stereo('InChI=1S/CH2N/c1-2/h1-2H'))
    assert set(chi.expand_stereo('InChI=1S/CH2N/c1-2/h1-2H')) == {
        'AMChI=1/CH2N/c1-2/h1-2H/b2-1-',
        'AMChI=1/CH2N/c1-2/h1-2H/b2-1+'}
    assert set(chi.expand_stereo('InChI=1S/C2/c1-2')) == {
        'InChI=1S/C2/c1-2'}
    print(chi.expand_stereo('InChI=1S/C3H3/c1-3-2/h1-3H'))
    assert set(chi.expand_stereo('InChI=1S/C3H3/c1-3-2/h1-3H')) == {
        'AMChI=1/C3H3/c1-3-2/h1-3H/b3-1-,3-2-',
        'AMChI=1/C3H3/c1-3-2/h1-3H/b3-1-,3-2+',
        'AMChI=1/C3H3/c1-3-2/h1-3H/b3-1+,3-2+'}


def test__canonical_enantiomer():
    """ test chi.canonical_enantiomer
             chi.is_canonical_enantiomer
    """
    # 1. SPECIES
    ref_ich = 'InChI=1S/C5H11O6/c1-3(9-6)5(11-8)4(2)10-7/h3-7H,1-2H3'

    # 1A. Full expansion -- includes non-canonical enantiomers
    print("Full species expansion:")
    for ich in chi.expand_stereo(ref_ich, enant=True):
        print(f"InChI:\n{ich}")
        is_can = chi.is_canonical_enantiomer(ich)
        print(f"Canonical? {is_can}")

        # Convert it to a canonical enantiomer like this
        ich = chi.canonical_enantiomer(ich)
        print(f"Canonical enantiomer:\n{ich}\n")
        assert chi.is_canonical_enantiomer(ich)

    # 1B. Restricted expansion -- includes only canonical enantiomers
    print("Restricted species expansion:")
    for ich in chi.expand_stereo(ref_ich, enant=False):
        print(ich)
        assert chi.is_canonical_enantiomer(ich)


def test__racemic():
    """ test amchi.racemic
    """
    ich = 'InChI=1S/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2-/m0/s1'
    assert (chi.racemic(ich) ==
            'InChI=1/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2-/s3')

    ach = 'AMChI=1/C6H11O/c1-3-4-6-5(2)7-6/h5-6H,1,3-4H2,2H3/t5-,6+/m0/s1'
    assert (chi.racemic(ach) ==
            'AMChI=1/C6H11O/c1-3-4-6-5(2)7-6/h5-6H,1,3-4H2,2H3/t5-,6+/m2/s3')

    assert chi.racemic(AR_ICH) == AR_ICH
    assert chi.racemic(C8H13O_ICH_NO_STEREO) == C8H13O_ICH_NO_STEREO


def test__is_complete():
    """ test inchi.is_complete
    """
    ich = 'InChI=1S/C4H6O2/c1-3-4(2,5-3)6-3/h1-2H3/t3-,4+'
    assert chi.is_complete(ich)


if __name__ == '__main__':
    # test__from_data()
    # test__formula_string()
    # test__version()
    # test__standard_form()
    # test__is_enantiomer()
    # test__reflect()
    # test__inchi_key()
    # test__formula()
    # test__is_standard_form()
    # test__has_multiple_components()
    # test__has_stereo()
    # test__same_connectivity()
    # test__equivalent()
    # test__sorted_()
    # test__recalculate()
    # test__filter_enantiomer_reactions()
    test__stereo()
    # test__join()
    # test__canonical_enantiomer()
    # test__racemic()
    test__is_complete()
