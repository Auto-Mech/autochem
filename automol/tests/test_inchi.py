""" test automol.inchi
"""
import numpy
from automol import inchi

AR_ICH = 'InChI=1S/Ar'
CH2O2_ICH = 'InChI=1S/CH2O2/c2-1-3/h1-2H/q+1/p+1'
# we can't handle isotopes because OpenBabel...
# C2H6O_ICH = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i2D/t2-/m1/s1'
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

C8H13O_ICH = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-/m0/s1')
C8H13O_ICH_NO_STEREO = 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3'


def test__from_data():
    """ test getters
    """
    assert AR_ICH == inchi.standard_form(inchi.from_data(
        fml_slyr=inchi.formula_sublayer(AR_ICH),
    ))

    assert CH2O2_ICH == inchi.standard_form(inchi.from_data(
        fml_slyr=inchi.formula_sublayer(CH2O2_ICH),
        main_dct=inchi.main_sublayers(CH2O2_ICH),
        char_dct=inchi.charge_sublayers(CH2O2_ICH),
    ))

    assert C2H2F2_ICH == inchi.standard_form(inchi.from_data(
        fml_slyr=inchi.formula_sublayer(C2H2F2_ICH),
        main_dct=inchi.main_sublayers(C2H2F2_ICH),
        ste_dct=inchi.stereo_sublayers(C2H2F2_ICH),
    ))

    assert C8H13O_ICH == inchi.standard_form(inchi.from_data(
        fml_slyr=inchi.formula_sublayer(C8H13O_ICH),
        main_dct=inchi.main_sublayers(C8H13O_ICH),
        ste_dct=inchi.stereo_sublayers(C8H13O_ICH),
    ))


def test__version():
    """ inchi.version
    """
    assert inchi.version(C2H2F2_ICH) == '1S'
    assert inchi.version(C2H2F2_ICH_NO_STEREO) == '1S'
    assert inchi.version(C2H2F2_ICH_STEREO_UNKNOWN) == '1'


def test__standard_form():
    """ test inchi.standard_form
    """
    assert (inchi.standard_form(C4H5F2O_ICH, remove_stereo=True) ==
            C4H5F2O_ICH_NO_STEREO)
    assert (inchi.standard_form(C2H2F2_ICH, remove_stereo=True) ==
            C2H2F2_ICH_NO_STEREO)
    assert (inchi.standard_form(C8H13O_ICH, remove_stereo=True) ==
            C8H13O_ICH_NO_STEREO)


def test__has_stereo():
    """ test inchi.has_stereo
    """
    assert inchi.has_stereo(C4H5F2O_ICH)
    assert inchi.has_stereo(C2H2F2_ICH)
    assert inchi.has_stereo(C8H13O_ICH)
    assert not inchi.has_stereo(
        inchi.standard_form(C4H5F2O_ICH, remove_stereo=True))
    assert not inchi.has_stereo(
        inchi.standard_form(C2H2F2_ICH, remove_stereo=True))
    assert not inchi.has_stereo(
        inchi.standard_form(C8H13O_ICH, remove_stereo=True))


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
    assert (inchi.recalculate(C2H2F2_ICH_NO_STEREO, force_stereo=True)
            == C2H2F2_ICH_STEREO_UNKNOWN)


if __name__ == '__main__':
    # test__from_data()
    # test__version()
    # test__join()
    # test__split()
    # test__standard_form()
    # test__has_stereo()
    test__argsort()
