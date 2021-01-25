"""
    test automol.rotor
"""

import automol


# CH2OOH
ZMA1 = (
    ('C', (None, None, None), (None, None, None),
     (None, None, None)),
    ('O', (0, None, None), ('R1', None, None),
     (2.638561158117497, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (2.0966833415410435, 1.9181232230723193, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
     (2.0974908160702483, 1.9169003262790676, 2.132511227259572)),
    ('O', (1, 0, 2), ('R4', 'A4', 'D4'),
     (2.4914580601748972, 1.8299443073232646, 1.797247663068249)),
    ('H', (4, 1, 0), ('R5', 'A5', 'D5'),
     (1.8737560578471502, 1.8329273093157175, 1.8571985995692384)))

# CH2ONO2
ZMA2 = (
    ('C', (None, None, None), (None, None, None),
     (None, None, None)),
    ('O', (0, None, None), ('R1', None, None),
     (2.602130817471561, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (2.0967335093840496, 1.9102333367742663, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
     (2.0963463040538937, 1.9113187477402638, 2.0944804672275423)),
    ('N', (1, 0, 2), ('R4', 'A4', 'D4'),
     (2.4982729475739274, 2.1411283443347364, 2.1041435884947064)),
    ('O', (4, 1, 0), ('R5', 'A5', 'D5'),
     (2.2538855452609683, 2.0850904561190813, 3.1497888417306883)),
    ('O', (4, 1, 5), ('R6', 'A6', 'D6'),
     (2.2586314858549517, 2.1117278699480875, 3.1466204524158154)))


def test__instab():
    """ test.automol.instab
    """

    instab_zmas1 = automol.instab.product_zmas(ZMA1)
    instab_zmas2 = automol.instab.product_zmas(ZMA2)

    print(instab_zmas1)
    print(instab_zmas2)

    assert instab_zmas1 == ref_instab_zmas1
    assert instab_zmas2 == ref_instab_zmas2


if __name__ == '__main__':
    test__instab()
