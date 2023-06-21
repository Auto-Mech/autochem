"""
    test automol.rotor
"""

import automol


# Transition state ZRXN object
C2H5OH_CH3_ZMA = (
    ('O', (None, None, None), (None, None, None),
     (None, None, None)),
    ('C', (0, None, None),
     ('R1', None, None), (2.684898037064155, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (1.8395472100737587, 1.8705001697540788, None)),
    ('C', (1, 0, 2), ('R3', 'A3', 'D3'),
     (2.852996269125268, 1.9268699770127682, 5.212206498615698)),
    ('H', (1, 0, 3), ('R4', 'A4', 'D4'),
     (2.064438713799504, 1.895091359334669, 2.135134545272428)),
    ('H', (1, 0, 3), ('R5', 'A5', 'D5'),
     (2.069338237845188, 1.888820446792396, 4.191772450395824)),
    ('H', (3, 1, 0), ('R6', 'A6', 'D6'),
     (2.0729166580744267, 1.9615649653998994, 0.9466768427806197)),
    ('H', (3, 1, 6), ('R7', 'A7', 'D7'),
     (2.083231604937746, 1.9334212394620052, 4.303178526038503)),
    ('H', (3, 1, 6), ('R8', 'A8', 'D8'),
     (2.0598428624094454, 1.9460020384003762, 2.149779853439449)),
    ('X', (6, 3, 1), ('R9', 'A9', 'D9'),
     (1.8897261254578284, 1.5707963267948966, 0.0)),
    ('C', (6, 9, 3), ('R10', 'A10', 'D10'),
     (3.1052664634613643, 1.6260531775487805, 3.4000422999151843)),
    ('H', (10, 6, 9), ('R11', 'A11', 'D11'),
     (2.0344920615786286, 1.3824840782104857, 4.5160641435691)),
    ('H', (10, 6, 11), ('R12', 'A12', 'D12'),
     (2.04136042741288, 1.1552191390225406, 2.5679911687678985)),
    ('H', (10, 6, 11), ('R13', 'A13', 'D13'),
     (2.1129899240028305, 2.7774451991784437, 4.005114953706984)))

C2H5OH_CH3_ZRXN_STR = """
reaction class: hydrogen abstraction
forward TS atoms:
  1: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  4: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  10: {symbol: X, implicit_hydrogens: 0, stereo_parity: null}
  11: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  14: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  2-4: {order: 1, stereo_parity: null}
  2-5: {order: 1, stereo_parity: null}
  2-6: {order: 1, stereo_parity: null}
  4-7: {order: 0.9, stereo_parity: null}
  4-8: {order: 1, stereo_parity: null}
  4-9: {order: 1, stereo_parity: null}
  7-10: {order: 0, stereo_parity: null}
  7-11: {order: 0.1, stereo_parity: null}
  11-12: {order: 1, stereo_parity: null}
  11-13: {order: 1, stereo_parity: null}
  11-14: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
- [11, 12, 13, 14]
backward TS atoms:
  1: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  2: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  6: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  7: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  8: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
backward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 0.9, stereo_parity: null}
  5-6: {order: 0.1, stereo_parity: null}
  6-7: {order: 1, stereo_parity: null}
  6-9: {order: 1, stereo_parity: null}
  6-10: {order: 1, stereo_parity: null}
  7-8: {order: 1, stereo_parity: null}
  7-11: {order: 1, stereo_parity: null}
  7-12: {order: 1, stereo_parity: null}
  8-13: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3, 4, 5]
- [6, 7, 8, 9, 10, 11, 12, 13]
"""
C2H5OH_CH3_ZRXN = automol.reac.from_string(C2H5OH_CH3_ZRXN_STR)


def test__ts():
    """ build rotors for a transition state
    """

    rotors = automol.rotor.from_zmatrix(
        C2H5OH_CH3_ZMA, zrxn=C2H5OH_CH3_ZRXN)

    assert automol.rotor.names(rotors) == (('D3',), ('D6',), ('D11',))
    assert automol.rotor.axes(rotors) == (((0, 1),), ((1, 3),), ((6, 10),))
    assert automol.rotor.groups(rotors) == (
        (((2,), (3, 4, 5, 6, 7, 8, 10, 11, 12, 13)),),
        (((0, 2, 4, 5), (6, 7, 8, 10, 11, 12, 13)),),
        (((0, 1, 2, 3, 4, 5, 7, 8), (11, 12, 13)),))
    assert automol.rotor.symmetries(rotors) == ((1,), (1,), (3,))
    assert automol.rotor.dimensions(rotors) == (1, 1, 1)


if __name__ == '__main__':
    test__ts()
