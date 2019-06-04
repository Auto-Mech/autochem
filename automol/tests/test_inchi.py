""" test automol.inchi
"""
from automol import inchi

AR_ICH = 'InChI=1S/Ar'

C2H2F2_ICH = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+'
C2H2F2_ICH_NO_STEREO = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H'
C2H2F2_ICH_STEREO_UNKNOWN = 'InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1?'

C8H13O_ICH = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-/m0/s1')
C8H13O_ICH_NO_ENANTIOMER = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-')
C8H13O_ICH_PARTIAL_STEREO = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/b5-3-/t8-/m0/s1')
C8H13O_ICH_NO_STEREO = 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3'
C8H13O_ICHS = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-/m0/s1',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4+/t8-/m0/s1',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3+,6-4-/t8-/m0/s1',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3+,6-4+/t8-/m0/s1'
)


def test__recalculate():
    """ inchi.recalculate
    """
    assert inchi.recalculate(C2H2F2_ICH_NO_STEREO) == C2H2F2_ICH_NO_STEREO
    assert (inchi.recalculate(C2H2F2_ICH_NO_STEREO, force_stereo=True)
            == C2H2F2_ICH_STEREO_UNKNOWN)


def test__is_closed():
    """ inchi.is_closed
    """
    assert inchi.is_closed(C8H13O_ICH) is True
    assert inchi.is_closed(C8H13O_ICH_PARTIAL_STEREO) is False
    assert inchi.is_closed(C8H13O_ICH_NO_STEREO) is True
    assert inchi.is_closed(C8H13O_ICH_NO_ENANTIOMER) is False


def test__prefix():
    """ inchi.prefix
    """
    assert inchi.prefix(C2H2F2_ICH) == 'InChI=1S'
    assert inchi.prefix(C2H2F2_ICH_NO_STEREO) == 'InChI=1S'
    assert inchi.prefix(C2H2F2_ICH_STEREO_UNKNOWN) == 'InChI=1'


def test__version():
    """ inchi.version
    """
    assert inchi.version(C2H2F2_ICH) == '1S'
    assert inchi.version(C2H2F2_ICH_NO_STEREO) == '1S'
    assert inchi.version(C2H2F2_ICH_STEREO_UNKNOWN) == '1'


def test__formula_layer():
    """ inchi.formula_layer
    """
    assert inchi.formula_layer(C2H2F2_ICH) == 'C2H2F2'
    assert (inchi.formula_layer('InChI=1S/2C2H5.Zn/c2*1-2;/h2*1H2,2H3;')
            == '2C2H5.Zn')


def test__key_layer():
    """ inchi.key_layer
    """
    assert inchi.key_layer(C2H2F2_ICH, 'c') == 'c3-1-2-4'
    assert inchi.key_layer(C2H2F2_ICH, 'h') == 'h1-2H'
    assert inchi.key_layer(C2H2F2_ICH, 'b') == 'b2-1+'
    assert inchi.key_layer(C2H2F2_ICH_STEREO_UNKNOWN, 'c') == 'c3-1-2-4'
    assert inchi.key_layer(C2H2F2_ICH_STEREO_UNKNOWN, 'h') == 'h1-2H'
    assert inchi.key_layer(C2H2F2_ICH_STEREO_UNKNOWN, 'b') == 'b2-1?'
    assert inchi.key_layer(C2H2F2_ICH_NO_STEREO, 'c') == 'c3-1-2-4'
    assert inchi.key_layer(C2H2F2_ICH_NO_STEREO, 'h') == 'h1-2H'
    assert inchi.key_layer(C2H2F2_ICH_NO_STEREO, 'b') is None


def test__key_layer_content():
    """ inchi.key_layer_content
    """
    assert inchi.key_layer_content(C2H2F2_ICH, 'c') == '3-1-2-4'
    assert inchi.key_layer_content(C2H2F2_ICH, 'h') == '1-2H'
    assert inchi.key_layer_content(C2H2F2_ICH, 'b') == '2-1+'
    assert (inchi.key_layer_content(C2H2F2_ICH_STEREO_UNKNOWN, 'c')
            == '3-1-2-4')
    assert (inchi.key_layer_content(C2H2F2_ICH_STEREO_UNKNOWN, 'h')
            == '1-2H')
    assert (inchi.key_layer_content(C2H2F2_ICH_STEREO_UNKNOWN, 'b')
            == '2-1?')
    assert inchi.key_layer_content(C2H2F2_ICH_NO_STEREO, 'c') == '3-1-2-4'
    assert inchi.key_layer_content(C2H2F2_ICH_NO_STEREO, 'h') == '1-2H'
    assert inchi.key_layer_content(C2H2F2_ICH_NO_STEREO, 'b') is None


def test__core_parent():
    """ inchi.core_parent
    """
    assert inchi.core_parent(AR_ICH) == AR_ICH
    assert (inchi.core_parent(C2H2F2_ICH)
            == 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H')
    assert (inchi.core_parent(C2H2F2_ICH_NO_STEREO)
            == 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H')
    assert (inchi.core_parent(C2H2F2_ICH_STEREO_UNKNOWN)
            == 'InChI=1/C2H2F2/c3-1-2-4/h1-2H')
    assert (inchi.core_parent(C8H13O_ICH_NO_STEREO)
            == 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3')
    assert (inchi.core_parent(C8H13O_ICH)
            == 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3')


def test__atom_stereo_elements():
    """ inchi.atom_stereo_elements
    """
    assert inchi.atom_stereo_elements(C8H13O_ICH_NO_STEREO) == ()
    assert inchi.atom_stereo_elements(C8H13O_ICH) == (('8', '-'),)


def test__bond_stereo_elements():
    """ inchi.bond_stereo_elements
    """
    assert inchi.bond_stereo_elements(C8H13O_ICH_NO_STEREO) == ()
    assert (inchi.bond_stereo_elements(C8H13O_ICH)
            == (('5-3', '-'), ('6-4', '-')))


def test__has_unknown_stereo_elements():
    """ inchi.has_unknown_stereo_elements
    """
    assert (inchi.has_unknown_stereo_elements(C8H13O_ICH)
            is False)
    assert (inchi.has_unknown_stereo_elements(C8H13O_ICH_PARTIAL_STEREO)
            is True)
    assert (inchi.has_unknown_stereo_elements(C8H13O_ICH_NO_STEREO)
            is True)
    assert (inchi.has_unknown_stereo_elements(C8H13O_ICH_NO_ENANTIOMER)
            is False)


def test__substereomers():
    """ inchi.substereomers
    """
    assert (inchi.substereomers(C8H13O_ICH_NO_STEREO)
            == C8H13O_ICHS)
    assert inchi.substereomers(C8H13O_ICH) == (C8H13O_ICH,)
