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
PROD_ZMAS1 = (
    # OH
    (('O', (None, None, None), (None, None, None),
      (None, None, None)),
     ('H', (0, None, None), ('R1', None, None),
      (1.847788877632985, None, None))),
    # CH2O
    (('C', (None, None, None), (None, None, None),
      (None, None, None)),
     ('O', (0, None, None), ('R1', None, None),
      (2.3142184670910955, None, None)),
     ('H', (0, 1, None), ('R2', 'A2', None),
      (2.081910294028233, 2.1334159256991865, None)),
     ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
      (2.0819105104503435, 2.1334159293636286, 3.141592653589793)))
)
BAD_PROD_ZMAS = (
    # O
    (('O', (None, None, None), (None, None, None),
     (None, None, None)),),
    # CH2OH
    (('C', (None, None, None), (None, None, None),
      (None, None, None)),
     ('O', (0, None, None), ('R1', None, None),
      (2.5598206437904127, None, None)),
     ('H', (0, 1, None), ('R2', 'A2', None),
      (2.0350667091575714, 2.008340705389423, None)),
     ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
      (2.0327423747809856, 2.0493848092212823, 3.141588744001008)),
     ('H', (1, 0, 2), ('R4', 'A4', 'D4'),
      (1.840826770401138, 1.9676117503333632, 3.141594966632109)))
)

INSTAB_ZRXN_STR = """
reaction class: beta scission
forward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  2-5: {order: 0.9, stereo_parity: null}
  5-6: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6]
backward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
backward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  2-5: {order: 0.1, stereo_parity: null}
  5-6: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3, 4]
- [5, 6]
"""

INSTAB_ZRXN_ZMA = (
    ('C', (None, None, None), (None, None, None),
     (None, None, None)),
    ('O', (0, None, None), ('R1', None, None),
     (2.6086764535623344, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (2.0966833415410435, 1.9070806550803825, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
     (2.0974908160702483, 1.8832236285512667, 2.104638003820193)),
    ('O', (1, 0, 2), ('R4', 'A4', 'D4'),
     (2.7476437083330785, 1.7824634600747093, 1.8136801951838313)),
    ('H', (4, 1, 0), ('R5', 'A5', 'D5'),
     (1.8447154354308446, 1.7658112135783564, 1.8571985995692388)))

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
PROD_ZMAS2 = (
    # NO2
    (('O', (None, None, None), (None, None, None),
      (None, None, None)),
     ('N', (0, None, None), ('R1', None, None),
      (2.6908226786788956, None, None)),
     ('O', (1, 0, None), ('R2', 'A2', None),
      (2.690824057320017, 1.8901908487803016, None))),
    # CH2O
    (('C', (None, None, None), (None, None, None),
      (None, None, None)),
     ('O', (0, None, None), ('R1', None, None),
      (2.314218121713856, None, None)),
     ('H', (0, 1, None), ('R2', 'A2', None),
      (2.081910657025832, 2.133415907619935, None)),
     ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
      (2.081910076584056, 2.133416208019412, 3.141594388471524)))
)

# CH3OO
ZMA3 = (
    ('C', (None, None, None), (None, None, None),
     (None, None, None)),
    ('O', (0, None, None), ('R1', None, None),
     (2.599568669917405, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (2.0720614763737584, 2.0034137221826973, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
     (2.1033532434977693, 1.96846586552775, 2.1568237825269154)),
    ('H', (0, 1, 2), ('R4', 'A4', 'D4'),
     (2.102894293957403, 1.8826666760922017, 4.18922576686547)),
    ('O', (1, 0, 2), ('R5', 'A5', 'D5'),
     (2.4633037728539113, 1.9322884438952062, 5.722447660200159)))


def test__prod_zmas():
    """ test.automol.reac.instability_product_zmas
    """

    instab_zmas1 = automol.reac.instability_product_zmas(ZMA1)
    instab_zmas2 = automol.reac.instability_product_zmas(ZMA2)
    instab_zmas3 = automol.reac.instability_product_zmas(ZMA3)

    for zma, ref_zma in zip(instab_zmas1, PROD_ZMAS1):
        assert automol.zmat.almost_equal(zma, ref_zma)
    for zma, ref_zma in zip(instab_zmas2, PROD_ZMAS2):
        assert automol.zmat.almost_equal(zma, ref_zma)
    assert not instab_zmas3


def test__transformation():
    """ test automol.zmat
    """

    ref_zrxn = automol.reac.from_string(INSTAB_ZRXN_STR)

    zrxn, zma = automol.reac.instability_transformation(ZMA1, PROD_ZMAS1)
    assert zrxn == ref_zrxn
    assert automol.zmat.almost_equal(zma, INSTAB_ZRXN_ZMA)


if __name__ == '__main__':
    # test__prod_zmas()
    test__transformation()
