""" test inchi_key
"""

from automol import inchi
from automol import inchi_key


C2H2F2_ICH = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+'
C2H2F2_ICH_NO_STEREO = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H'
C2H2F2_ICH_STEREO_UNKNOWN = 'InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1?'


def test__first_hash():
    """ inchi.key.first_hash()
    """
    print(inchi.inchi_key(C2H2F2_ICH))
    assert (inchi_key.first_hash(
        inchi.inchi_key(C2H2F2_ICH)) == 'WFLOTYSKFUPZQB')
    assert (inchi_key.first_hash(
        inchi.inchi_key(C2H2F2_ICH_NO_STEREO)) == 'WFLOTYSKFUPZQB')
    assert (inchi_key.first_hash(
        inchi.inchi_key(C2H2F2_ICH_STEREO_UNKNOWN)) == 'WFLOTYSKFUPZQB')


def test__second_hash():
    """ inchi.key.second_hash()
    """
    assert (inchi_key.second_hash(
        inchi.inchi_key(C2H2F2_ICH)) == 'OWOJBTED')
    assert (inchi_key.second_hash(
        inchi.inchi_key(C2H2F2_ICH_NO_STEREO)) == 'UHFFFAOY')
    assert (inchi_key.second_hash(
        inchi.inchi_key(C2H2F2_ICH_STEREO_UNKNOWN)) == 'HXYFBOIP')


def test__version_indicator():
    """ inchi.key.version_indicator()
    """
    assert (inchi_key.version_indicator(
        inchi.inchi_key(C2H2F2_ICH)) == 'SA')
    assert (inchi_key.version_indicator(
        inchi.inchi_key(C2H2F2_ICH_NO_STEREO)) == 'SA')
    assert (inchi_key.version_indicator(
        inchi.inchi_key(C2H2F2_ICH_STEREO_UNKNOWN)) == 'NA')


def test__protonation_indicator():
    """ inchi.key.protonation_indicator()
    """
    ich1 = 'InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)'
    ich2 = 'InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)/p-1'
    ich3 = 'InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)/p+1'
    assert inchi_key.protonation_indicator(inchi.inchi_key(ich1)) == 'N'
    assert inchi_key.protonation_indicator(inchi.inchi_key(ich2)) == 'M'
    assert inchi_key.protonation_indicator(inchi.inchi_key(ich3)) == 'O'


if __name__ == "__main__":
    print("Hello, world!")
    test__first_hash()