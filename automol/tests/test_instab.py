"""
    test automol instability functions
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
  1: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  2: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  5: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  2-5: {order: 0.9, stereo_parity: null}
  5-6: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6]
backward TS atoms:
  1: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  2: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  5: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
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


# STE
STE_ZMA = (
    ('C', (None, None, None), (None, None, None),
        (None, None, None)),
    ('C', (0, None, None), ('R1', None, None),
        (2.879329461537935, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
        (2.0685853010041964, 1.936203621625536, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
        (2.0681998385698437, 1.9437991078452497, 2.1080659988667785)),
    ('H', (0, 1, 2), ('R4', 'A4', 'D4'),
        (2.068569230988053, 1.9221145780688327, 4.196383631606691)),
    ('C', (1, 0, 2), ('R5', 'A5', 'D5'),
        (2.8951853965121304, 1.9564542588074854, 5.252745506231529)),
    ('H', (1, 0, 5), ('R6', 'A6', 'D6'),
        (2.0727280210037247, 1.8998644928532522, 2.130145952014642)),
    ('H', (1, 0, 5), ('R7', 'A7', 'D7'),
        (2.0740453471105984, 1.8851505787673934, 4.161901847872003)),
    ('C', (5, 1, 0), ('R8', 'A8', 'D8'),
        (2.815645112515456, 1.92909277095781, 3.1002400398684897)),
    ('O', (5, 1, 8), ('R9', 'A9', 'D9'),
        (2.710554414591679, 1.8963614225749115, 2.0628190070391694)),
    ('H', (5, 1, 8), ('R10', 'A10', 'D10'),
        (2.0738054178739627, 1.9061743903687773, 4.152460126257847)),
    ('O', (8, 5, 1), ('R11', 'A11', 'D11'),
        (2.57957319683146, 2.0630955394242707, 3.9898373226106236)),
    ('H', (8, 5, 11), ('R12', 'A12', 'D12'),
        (2.037337829840982, 2.161046773185505, 3.342493863112154)),
    ('O', (11, 8, 5), ('R13', 'A13', 'D13'),
        (2.7396264771381125, 1.8404255371622773, 3.200314265734847)),
    ('H', (13, 11, 8), ('R14', 'A14', 'D14'),
        (1.8396656471039372, 1.6762185407776191, 4.304713884757464)),
    ('O', (9, 5, 1), ('R15', 'A15', 'D15'),
        (2.759974384617291, 1.8588964031436905, 2.682937554897634)),
    ('H', (15, 9, 5), ('R16', 'A16', 'D16'),
        (1.8427096812519266, 1.6838436964284405, 1.7402981108559878))
)

STE_PROD_ZMAS = (
    # OH
    (('O', (None, None, None), (None, None, None),
        (None, None, None)),
     ('H', (0, None, None), ('R1', None, None),
        (1.8477888298795644, None, None))),
    # stereo CHO
    (('C', (None, None, None), (None, None, None),
        (None, None, None)),
     ('C', (0, None, None), ('R1', None, None),
        (2.874942617600433, None, None)),
     ('H', (0, 1, None), ('R2', 'A2', None),
        (2.068486178501242, 1.9368605648050443, None)),
     ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
        (2.068867675123655, 1.920464980301257, 2.084268267429751)),
     ('H', (0, 1, 2), ('R4', 'A4', 'D4'),
        (2.06670996549009, 1.9490219370816204, 4.169403811281169)),
     ('C', (1, 0, 2), ('R5', 'A5', 'D5'),
        (2.892054198885529, 1.978109466699197, 0.9931015869119763)),
     ('H', (1, 0, 5), ('R6', 'A6', 'D6'),
        (2.0731098941442894, 1.895113111124007, 4.156030631730728)),
     ('H', (1, 0, 5), ('R7', 'A7', 'D7'),
        (2.0735983863560614, 1.9083973887518564, 2.148304203029353)),
     ('C', (5, 1, 0), ('R8', 'A8', 'D8'),
        (2.8781935635948996, 1.9361198955628287, 1.2783798744140946)),
     ('O', (5, 1, 8), ('R9', 'A9', 'D9'),
        (2.6986014630716118, 1.8964799675479287, 2.0837938021682145)),
     ('H', (5, 1, 8), ('R10', 'A10', 'D10'),
        (2.07270883018819, 1.9463407836628417, 4.2023952224973335)),
     ('O', (8, 5, 1), ('R11', 'A11', 'D11'),
        (2.3208847722366746, 2.189986883591254, 2.1757813432302395)),
     ('H', (8, 5, 11), ('R12', 'A12', 'D12'),
        (2.0841000255210202, 2.0077189801794337, 3.0900070546866)),
     ('O', (9, 5, 1), ('R13', 'A13', 'D13'),
        (2.7550592578760584, 1.8564133553322275, 1.5438089821509182)),
     ('H', (13, 9, 5), ('R14', 'A14', 'D14'),
        (1.8433348334954938, 1.6767577456458034, 2.2296868763866917)))
)

STE_INSTAB_ZRXN_STR = """
reaction class: beta scission
forward TS atoms:
  1: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  6: {symbol: C, implicit_hydrogens: 0, stereo_parity: true}
  7: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  9: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  10: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  12: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  14: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  15: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  16: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  17: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 1, stereo_parity: null}
  2-6: {order: 1, stereo_parity: null}
  2-7: {order: 1, stereo_parity: null}
  2-8: {order: 1, stereo_parity: null}
  6-9: {order: 1, stereo_parity: null}
  6-10: {order: 1, stereo_parity: null}
  6-11: {order: 1, stereo_parity: null}
  9-12: {order: 1, stereo_parity: null}
  9-13: {order: 1, stereo_parity: null}
  10-16: {order: 1, stereo_parity: null}
  12-14: {order: 0.9, stereo_parity: null}
  14-15: {order: 1, stereo_parity: null}
  16-17: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
backward TS atoms:
  1: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  6: {symbol: C, implicit_hydrogens: 0, stereo_parity: true}
  7: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  9: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  10: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  12: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  14: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  15: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  16: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  17: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
backward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 1, stereo_parity: null}
  2-6: {order: 1, stereo_parity: null}
  2-7: {order: 1, stereo_parity: null}
  2-8: {order: 1, stereo_parity: null}
  6-9: {order: 1, stereo_parity: null}
  6-10: {order: 1, stereo_parity: null}
  6-11: {order: 1, stereo_parity: null}
  9-12: {order: 1, stereo_parity: null}
  9-13: {order: 1, stereo_parity: null}
  10-14: {order: 1, stereo_parity: null}
  12-16: {order: 0.1, stereo_parity: null}
  14-15: {order: 1, stereo_parity: null}
  16-17: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
- [16, 17]

"""


def test__prod_zmas():
    """ test.automol.reac.instability_product_zmas
    """

    instab_zmas1 = automol.reac.instability_product_zmas(ZMA1)
    instab_zmas2 = automol.reac.instability_product_zmas(ZMA2)
    instab_zmas3 = automol.reac.instability_product_zmas(ZMA3)

    for zma, ref_zma in zip(instab_zmas1, PROD_ZMAS1):
        assert (automol.zmat.vmatrix(zma) ==
                automol.zmat.vmatrix(ref_zma))
        assert (automol.zmat.graph(zma) ==
                automol.zmat.graph(ref_zma))
    for zma, ref_zma in zip(instab_zmas2, PROD_ZMAS2):
        assert (automol.zmat.vmatrix(zma) ==
                automol.zmat.vmatrix(ref_zma))
        assert (automol.zmat.graph(zma) ==
                automol.zmat.graph(ref_zma))
    assert not instab_zmas3


def test__transformation():
    """ test automol.zmat
    """

    ref_zrxn = automol.reac.from_old_string(INSTAB_ZRXN_STR)

    zrxn, zma = automol.reac.instability_transformation(ZMA1, PROD_ZMAS1)
    assert ref_zrxn is not None
    assert zrxn is not None
    assert zma is not None
    # assert zrxn == ref_zrxn
    # assert automol.zmat.almost_equal(zma, INSTAB_ZRXN_ZMA)


def test__stereo():
    """ test the stereo functions
    """

    ref_zrxn = automol.reac.from_old_string(STE_INSTAB_ZRXN_STR)
    zrxn, _ = automol.reac.instability_transformation(STE_ZMA, STE_PROD_ZMAS)
    print(automol.reac.string(zrxn))
    assert ref_zrxn is not None
    assert zrxn is not None
    # assert zrxn == ref_zrxn


if __name__ == '__main__':
    # test__prod_zmas()
    test__transformation()
    test__stereo()
