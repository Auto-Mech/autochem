"""
    test automol.rotor
"""

import numpy
import automol


SPC_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(
        automol.smiles.inchi('CCCO')))

TS_ZMA = (
    (('C', (None, None, None), (None, None, None)),
     ('C', (0, None, None), ('R1', None, None)),
     ('H', (0, 1, None), ('R2', 'A2', None)),
     ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
     ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
     ('O', (1, 0, 2), ('R5', 'A5', 'D5')),
     ('H', (1, 0, 5), ('R6', 'A6', 'D6')),
     ('H', (1, 0, 5), ('R7', 'A7', 'D7')),
     ('H', (5, 1, 0), ('R8', 'A8', 'D8')),
     ('X', (7, 1, 0), ('R9', 'A9', 'D9')),
     ('C', (7, 9, 1), ('R10', 'A10', 'D10')),
     ('H', (10, 7, 9), ('R11', 'A11', 'D11')),
     ('H', (10, 11, 7), ('R12', 'A12', 'D12')),
     ('H', (10, 11, 12), ('R13', 'A13', 'D13'))),
    {'R1': 2.86047843610551,
     'R2': 2.06868318953868, 'A2': 1.92742643417190,
     'R3': 2.07189572395196, 'A3': 1.94002247538354, 'D3': 4.17755377506205,
     'R4': 2.07454134052760, 'A4': 1.93385797246550, 'D4': 2.07811118047208,
     'R5': 2.63465616411330, 'A5': 2.00629786306952, 'D5': 3.13352050579931,
     'R6': 2.06509270990031, 'A6': 1.97849825874376, 'D6': 4.11924756074067,
     'R7': 2.41053464563400, 'A7': 1.85996596792382, 'D7': 2.11338428465489,
     'R8': 1.82415262890444, 'A8': 1.88251213120108, 'D8': 0.89541150213840,
     'R11': 2.0599904493615, 'R12': 2.0581007232361, 'A12': 1.9886403670271,
     'R13': 2.0584786684612, 'A13': 1.9847552641121, 'D13': 2.3199789082134,
     'D12': 1.9638130584174, 'D10': 3.2229703752932, 'A11': 1.8333968207207,
     'A10': 1.5739204661559, 'D11': 2.3565871892740, 'R10': 2.7297093882238,
     'D9': 3.14159265358979, 'A9': 1.57079632679489, 'R9': 1.88972612545782})
TS_TORS_NAMES = ('D5', 'D8', 'D11')
TS_TRANS = ('', frozenset({7, 10}), frozenset({1, 7}))


TORS_DCT = {
    'D12': {
        'group1': (5, 6, 7),
        'group2': (1, 2),
        'axis': (3, 4),
        'span': 6.28,
        'sym_num': 1,
    },
    'D10': {
        'group1': (1, 2, 3, 4, 5, 6, 7, 8),
        'group2': (11, 12),
        'axis': (9, 10),
        'span': 3.14,
        'sym_num': 2
    }
}


def test__tors_string():
    """ test.automol.rotor.tors.string
        test.automol.rotor.tors.from_string
    """

    tors_str = automol.rotor.tors.string(TORS_DCT)
    print(tors_str)
    tors_dct = automol.rotor.tors.from_string(tors_str)
    print(tors_dct)


def test__():
    """ rotor
    """

    rotors = automol.rotor.from_zma(SPC_ZMA)
    print(rotors)


if __name__ == '__main__':
    # test__tors_inf()
    # test__tors_string()
    test__()