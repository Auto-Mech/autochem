"""
    test automol.rotor
"""

import automol


ZMA1 = ()  # CH2OOH
ZMA2 = ()  # CH2ONO2


def test__instab():
    """ test.automol.instab
    """

    instab_zmas1 = automol.instab.product_zmas(ZMA1)
    instab_zmas2 = automol.instab.product_zmas(ZMA2)

    assert instab_zmas1 == ref_instab_zmas1
    assert instab_zmas2 == ref_instab_zmas2


if __name__ == '__main__':
    test__instab()
