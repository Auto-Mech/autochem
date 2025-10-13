"""
  Tests the functional group
"""

import automol

SPCS_CHECKS_SMI = {
    'C6H5': 'C=1C=C[C]=CC=1',
    'AMN': 'Cc1ccc2ccccc2c1',
    'INDENYL': 'C=1C=CC=2C=C[CH]C=2C=1',
    'C12H8': 'C=1C=C2C=CC=C3C=CC(C=1)=C32',
    'C9H7O': '[O]C=1CC=C2C=CC=CC2=1',
    'C5H4O': 'O=C1C=CC=C1',
    'OC6H4CH3': 'CC1=CC=CC=C1[O]',
    'HOC6H4CH3': 'CC1=CC=CC(=C1)O',
    'C6H5CH2OOH': 'OOCC=1C=CC=CC=1',
    'BZFUR': 'C=1C=CC2=C(C=1)C=CO2',
    'C10H7CH2': '[CH2]C1=CC=CC=2C=CC=CC1=2',
    'C6H5C2H3': 'C=CC1=CC=CC=C1',
    'C10H9': 'C=1C=CC=2[CH]CC=CC=2C=1',
    'CH3C6H4': 'CC1=CC=CC=[C]1',
    'C6H5CCC6H5': 'C=1C=CC(=CC=1)C#CC2=CC=CC=C2',
    'C6H5C3H3-A': 'C=C=CC1=CC=CC=C1',
    'C6H5C3H3-P': 'C#CCC1=CC=CC=C1',
    'CYC5H7': 'C1=CCC[CH]1',
    'MEINDENYL': 'C[C]1C=CC=2C=CC=CC=21',
    'BENZOFLUORENE':'C=1C=CC2=C(C=1)C=CC=3CC=4C=CC=CC=4C=32',
    'C6H5OCH3': 'COC1=CC=CC=C1',
    'CATECHOL': 'OC=1C=CC=CC=1O',
    'SALICALD': 'O=CC1=CC=CC=C1O',
}


for spc, smi in SPCS_CHECKS_SMI.items():
    SPCS_CHECKS_SMI[spc] = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi(smi)))

SPCS_GRPS = {
    'BZFUR': {'FUR-M': ((4, 5, 8, 7, 6),) ,
        'A1-M': ((0, 1, 3, 7, 6, 2),) ,},
    'AMN': {'A1-M': ((5, 6, 9, 10, 7, 8), (1, 2, 4, 10, 9, 3)) ,
            'A1,CH3-M': ((5, 6, 9, 10, 7, 8, 0),)},
    'C6H5': {'A1-R': ((0, 1, 3, 5, 4, 2),)},
    'INDENYL': {
        'C5-RSR': ((2, 5, 7, 8, 6),) ,
     },
    'C12H8': {
        'A1-M': ((0, 2, 8, 11, 9, 4), (1, 3, 8, 11, 10, 5)) ,
        'A1,C2H3-M': ((0, 2, 8, 11, 9, 4, 6, 7),),
     },
    'C9H7O': {'C5O-RSR': ((4, 5, 8, 7, 6, 9),) ,},
    'C5H4O': {'C5O-M': ((0, 1, 3, 4, 2, 5),) ,},
    'OC6H4CH3': {'A1O-RSR': ((1, 2, 4, 6, 5, 3, 7),) ,
                 'A1O,CH3-RSR': ((1, 2, 4, 6, 5, 3, 7, 0),),
            },
    'HOC6H4CH3': {
        'A1-M': ((1, 2, 5, 4, 6, 3),) ,
        'A1,CH3-M':  ((1, 2, 5, 4, 6, 3, 0),),
        'A1,OH-M': ((1, 2, 5, 4, 6, 3, 7),),
        'A1,OH,CH3-M': ((1, 2, 5, 4, 6, 3, 7, 0),),
        },     
    'C6H5CH2OOH': {'A1-M': ((0, 1, 3, 6, 4, 2),) ,
        },
    'C10H7CH2': {'A1CH2-RSR': ((3, 4, 8, 10, 9, 6, 0),) ,},
    'C6H5C2H3': {'A1-M': ((2, 3, 5, 7, 6, 4),) ,
                 'A1,C2H3-M': ((2, 3, 5, 7, 6, 4, 0, 1),),
        },
    'C10H9': {'A1CH2-RSR': ((0, 1, 5, 9, 8, 4, 7), (0, 1, 5, 9, 8, 4, 6)) ,
              },
    'CH3C6H4': {'A1-R': ((1, 2, 4, 6, 5, 3),) ,
                'A1,CH3-R': ((1, 2, 4, 6, 5, 3, 0),),
        },
    'C6H5CCC6H5': {
        'A1-M': ((1, 4, 8, 13, 9, 5), (0, 2, 6, 12, 7, 3)) ,
        'A1,C2H-M': ((1, 4, 8, 13, 9, 5, 10, 11),),
        },
    'C6H5C3H3-A': {'A1-M': ((2, 3, 6, 8, 7, 4),) ,
        'A1,C3.DD-M': ((2, 3, 6, 8, 7, 4, 0, 1, 5),)},
    'C6H5C3H3-P': {'A1-M': ((2, 3, 6, 8, 7, 4),) ,
        'A1,C3.ST-M': ((2, 3, 6, 8, 7, 4, 0, 1, 5),)
        },
    'CYC5H7': {'C5H2-RSR': ((0, 1, 3, 4, 2),) ,},
    'MEINDENYL': {'C5-RSR': ((5, 6, 8, 9, 7),) ,
                  'C5,CH3-RSR': ((5, 6, 8, 9, 7, 0),)
        },
    'BENZOFLUORENE': {'C5-M': ((10, 12, 15, 16, 13),) ,
        'A1-M': ((0, 2, 6, 14, 11, 4), (8, 9, 13, 16, 14, 11), (1, 3, 7, 15, 12, 5)) ,
        },
    'C6H5OCH3': {'A1-M': ((1, 2, 4, 6, 5, 3),) ,
        'A1,OCH3-M': ((1, 2, 4, 6, 5, 3, 7, 0),)},
    'CATECHOL': {'A1-M': ((0, 1, 3, 5, 4, 2),) ,
        'A1,OH-M': ((0, 1, 3, 5, 4, 2, 6), (0, 1, 3, 5, 4, 2, 7)),
        'A1,OH,OH-M': ((0, 1, 3, 5, 4, 2, 6, 7),)},
    'SALICALD': {
        'A1-M': ((0, 1, 3, 6, 5, 2),) ,
        'A1,OH-M': ((0, 1, 3, 6, 5, 2, 8),) ,
        'A1,CHO-M': ((0, 1, 3, 6, 5, 2, 4, 7),), 
        'A1,OH,CHO-M': ((0, 1, 3, 6, 5, 2, 8, 4, 7),), 
        },
}

def test_super_functional_group_dct():
    
    for SPC, DCT in SPCS_GRPS.items():
        print(SPC)
        gra = SPCS_CHECKS_SMI[SPC]
        fgrps = automol.graph.SuperFunctionalGroup()
        fgrps.assign_grps(gra)
        for key, val in fgrps.sup_grps.items():
            if len(val) > 0:
                if key not in DCT.keys():
                    print(key, val)
                assert val == DCT[key]

            
if __name__ == '__main__':
    test_super_functional_group_dct()
