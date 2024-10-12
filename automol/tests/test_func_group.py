"""
  Tests the functional group
"""

import automol
from automol.graph import kekule
from automol.graph import explicit

SPCS_CHECKS_SMI = {
    'C2H6': 'CC',
    'C4H10': 'CCCC',
    'C4H6': 'C=CC=C',
    'C4H8': 'CC=CC',
    'C4H8O': 'CC1OC1C',
    'C4H2': 'C#CC#C',
    'C2H5OH': 'CCO',
    'C2H5CL': 'CCCl',
    'C2H5SH': 'CCS',
    'CH3OCH3': 'COC',
    'CYC_ETHER': 'CC1CO1',
    'C2H5OOH': 'CCOO',
    'C2H5OO': 'CCO[O]',
    'CCOCCO': 'CC(=O)CC=O',
    'C2H5CO_OH': 'CC(=O)O',
    'CCOOC': 'CC(=O)OC',
    'C3H4-A': 'C=C=C',
    'C3H4-P': 'C#CC',
    'C3H5-A': 'C=C[CH2]',
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
    'C2H6': {'alkane': ((0, 1),) ,
        'methyl': ((1, 5, 7, 6), (0, 3, 2, 4)) ,},
    'C4H10': {'alkane': ((2, 3), (0, 2), (1, 3)) ,
        'methyl': ((1, 9, 7, 8), (0, 6, 4, 5)) ,},
    'C4H6': {'alkene': ((0, 2), (1, 3)) ,},
    'C4H8': {'alkane': ((0, 2), (1, 3)),
             'alkene': ((2, 3),),
             'methyl': ((1, 9, 7, 8), (0, 6, 4, 5))},
    'C4H8O': {'alkane': ((2, 3), (0, 2), (1, 3)),
              'cyclic_ether': ((2, 4, 3),),
              'alkoxy_oc': ((4, 3), (4, 2)),
              'methyl': ((1, 9, 10, 8), (0, 6, 5, 7)),},
    'C4H2': {'alkyne': ((0, 2), (1, 3)) ,
        },
    'C2H5OH': {'alkane': ((0, 1),) ,
        'alcohol': ((2, 8),) ,
        'methyl': ((0, 3, 4, 5),) ,},
    'C2H5CL': {'alkane': ((0, 1),) ,
        'halide': ((1, 2),) ,
        'methyl': ((0, 3, 4, 5),) ,},
    'C2H5SH': {'alkane': ((0, 1),) ,
        'thiol': ((1, 2, 8),) ,
        'methyl': ((0, 3, 4, 5),) ,},
    'CH3OCH3': {'alkoxy_oc': ((2, 1), (2, 0)) ,
        'ether': ((0, 2, 1),) ,
        'methyl': ((1, 7, 8, 6), (0, 3, 4, 5)) ,},
    'CYC_ETHER': {'alkane': ((1, 2), (0, 2)) ,
        'alkoxy_oc': ((3, 2), (3, 1)) ,
        'cyclic_ether': ((1, 3, 2),) ,
        'methyl': ((0, 6, 4, 5),) ,},
    'C2H5OOH': {'alkane': ((0, 1),) ,
        'alkoxy_oc': ((3, 1),) ,
        'hydroperoxy': ((1, 3, 2, 9),) ,
        'methyl': ((0, 6, 4, 5),) ,},
    'C2H5OO': {'alkane': ((0, 1),) ,
        'peroxy': ((1, 3, 2),) ,
        'methyl': ((0, 6, 4, 5),) ,},
    'CCOCCO': {'alkane': ((1, 2), (0, 3), (1, 3)) ,
        'aldehyde': ((2, 4),) ,
        'ketone': ((3, 5),) ,
        'methyl': ((0, 6, 7, 8),) ,},
    'C2H5CO_OH': {'alkane': ((0, 1),) ,
        'alcohol': ((3, 7),) ,
        'ketone': ((1, 2),) ,
        'methyl': ((0, 6, 4, 5),) ,},
    'CCOOC': {'alkane': ((0, 2),) ,
        'alkoxy_oc': ((4, 1),) ,
        'ester': ((3, 2, 4, 1),) ,
        'methyl': ((1, 9, 8, 10), (0, 6, 5, 7)) ,},
    'C3H4-A': {'allene': ((0, 1, 2),) ,},
    'C3H4-P': {'alkyne': ((0, 2),) , 'propyne': ((0, 1, 2),) ,},
    'C3H5-A': {'allyl': ((2, 0, 1),) ,},
    'C6H5': {
         'aromatic': ((0, 1, 3, 5, 4, 2),),
         'phenyl':  ((0, 1, 3, 5, 4, 2),),
     },
    'AMN':  {
        'methyl': ((0, 12, 11, 13),),
        'benzene': ((5, 6, 9, 10, 7, 8), (1, 2, 4, 10, 9, 3)) ,
        'aromatic': ((5, 6, 9, 10, 7, 8), (1, 2, 4, 10, 9, 3)),
     },
    'INDENYL': {'aromatic': ((0, 1, 4, 8, 7, 3),) ,
        'cyclopentadienyl': ((2, 5, 7, 8, 6),) ,
     },
    'C12H8': {
        'alkene': ((6, 7),), 
        'aromatic': ((0, 2, 8, 11, 9, 4), (1, 3, 8, 11, 10, 5)),
        'benzene': ((0, 2, 8, 11, 9, 4), (1, 3, 8, 11, 10, 5)),
        
     },
    'C9H7O': {'alkoxy': ((8, 9),) ,
        'aromatic': ((0, 1, 3, 7, 6, 2),) ,
        'cyclopentadienonyl': ((4, 5, 8, 7, 6, 9),) ,},
    'C5H4O': {'alkene': ((0, 2), (1, 3)) ,
        'ketone': ((4, 5),) , #exclude?
        'cyclopentadienone': ((0, 1, 3, 4, 2, 5),) ,},
    'OC6H4CH3': {
        'alkoxy': ((6, 7),) ,
        'methyl': ((0, 10, 8, 9),) ,
        'aromatic': ((1, 2, 4, 6, 5, 3),) ,
        'phenoxy': ((1, 2, 4, 6, 5, 3, 7),) ,
        },
    'HOC6H4CH3': {
        'alcohol': ((7, 15),) ,
        'methyl': ((0, 10, 8, 9),) ,
        'benzene': ((1, 2, 5, 4, 6, 3),) ,
        'aromatic': ((1, 2, 5, 4, 6, 3),) ,},
    'C6H5CH2OOH': {'alkoxy_oc': ((8, 5),) ,
        'hydroperoxy': ((5, 8, 7, 16),) ,
        'benzene': ((0, 1, 3, 6, 4, 2),) ,
        'aromatic': ((0, 1, 3, 6, 4, 2),) ,},
    'BZFUR': {'alkene': ((4, 5),) ,
        'ether': ((5, 8, 7),) ,        #exclude?
        'cyclic_ether': ((5, 8, 7),) , #exclude?
        'benzene': ((0, 1, 3, 7, 6, 2),) ,
        'aromatic': ((0, 1, 3, 7, 6, 2),) ,
        'furan': ((4, 5, 8, 7, 6),) ,
        },
    'C10H7CH2': {'allyl': ((8, 0, 4),) ,
        'aromatic': ((3, 4, 8, 10, 9, 6), (1, 2, 7, 10, 9, 5)) ,
        'benzyl': ((3, 4, 8, 10, 9, 6, 0),) ,
        },
    'C6H5C2H3': {'alkene': ((0, 1),) ,
        'benzene': ((2, 3, 5, 7, 6, 4),) ,
        'aromatic': ((2, 3, 5, 7, 6, 4),) ,
        },
    'C10H9': {'allyl': ((9, 8, 7),) ,
        'aromatic': ((0, 1, 5, 9, 8, 4),) ,
        'benzyl': ((0, 1, 5, 9, 8, 4, 7), (0, 1, 5, 9, 8, 4, 6)) ,
        },
    'CH3C6H4':{'methyl': ((0, 7, 8, 9),) ,
        'phenyl': ((1, 2, 4, 6, 5, 3),) ,
        'aromatic': ((1, 2, 4, 6, 5, 3),) ,},
    'C6H5CCC6H5': {'alkyne': ((10, 11),) ,
        'benzene': ((1, 4, 8, 13, 9, 5), (0, 2, 6, 12, 7, 3)) ,
        'aromatic': ((1, 4, 8, 13, 9, 5), (0, 2, 6, 12, 7, 3)) ,
        },
    'C6H5C3H3-A': {'allene': ((0, 1, 5),) ,
        'benzene': ((2, 3, 6, 8, 7, 4),) ,
        'aromatic': ((2, 3, 6, 8, 7, 4),) ,},
    'C6H5C3H3-P': {'alkyne': ((0, 1),) ,
        'propyne': ((0, 1, 5),) ,
        'benzene': ((2, 3, 6, 8, 7, 4),) ,
        'aromatic': ((2, 3, 6, 8, 7, 4),) ,},
    'CYC5H7': {'allyl': ((0, 1, 2),) ,
        'cyclopentenyl': ((0, 1, 3, 4, 2),) ,
        },
    'MEINDENYL': {
        'methyl': ((0, 10, 12, 11),) ,
        'aromatic': ((1, 2, 4, 9, 8, 3),) ,
        'cyclopentadienyl': ((5, 6, 8, 9, 7),) ,
        },
    'BENZOFLUORENE': {'benzene': ((0, 2, 6, 14, 11, 4), (8, 9, 13, 16, 14, 11), (1, 3, 7, 15, 12, 5)) ,
        'aromatic': ((0, 2, 6, 14, 11, 4), (8, 9, 13, 16, 14, 11), (1, 3, 7, 15, 12, 5)) ,
        'cyclopentadiene': ((10, 12, 15, 16, 13),) ,
        },
    'C6H5OCH3': {'alkoxy_oc': ((7, 0),) ,
        'ether': ((0, 7, 6),) ,
        'methyl': ((0, 10, 8, 9),) ,
        'benzene': ((1, 2, 4, 6, 5, 3),) ,
        'aromatic': ((1, 2, 4, 6, 5, 3),) ,
        },
    'CATECHOL': {'alcohol': ((6, 12), (7, 13)) ,
        'benzene': ((0, 1, 3, 5, 4, 2),) ,
        'aromatic': ((0, 1, 3, 5, 4, 2),) ,},
    'SALICALD': {'alcohol': ((8, 14),) ,
        'aldehyde': ((4, 7),) ,
        'benzene': ((0, 1, 3, 6, 5, 2),) ,
        'aromatic': ((0, 1, 3, 6, 5, 2),) ,
        },
}


def test_functional_group_dct():
    """ test automol.graph.functional_group_dct
    """
    for SPC, DCT in SPCS_GRPS.items():
        gra = SPCS_CHECKS_SMI[SPC]
        fgrps = automol.graph.functional_group_dct(gra)
        for key, val in DCT.items():
            assert val == fgrps[key]


def __sites():
    """ test automol.graph._func_group.'alkene'_sites
        test automol.graph._func_group.'alkene'_sites
    """
    assert automol.graph.alkene_sites(kekule(explicit(SPCS_CHECKS_SMI['C4H6']))) == ((0, 2), (1, 3))
    assert automol.graph.alkyne_sites(kekule(explicit(SPCS_CHECKS_SMI['C4H2']))) == ((0, 2), (1, 3))
    assert not automol.graph.alkene_sites(kekule(explicit(SPCS_CHECKS_SMI['C4H2'])))
    assert not automol.graph.alkyne_sites(kekule(explicit(SPCS_CHECKS_SMI['C4H6'])))


def test_species_types():
    """ test automol.graph._func_group.hydrocarbon_species
             automol.graph._func_group.radical_species
    """

    assert automol.graph.is_hydrocarbon_species(SPCS_CHECKS_SMI['C2H6'])
    assert not automol.graph.is_hydrocarbon_species(SPCS_CHECKS_SMI['C2H5OH'])

    assert automol.graph.is_radical_species(SPCS_CHECKS_SMI['C2H5OO'])
    assert not automol.graph.is_radical_species(SPCS_CHECKS_SMI['C2H5OH'])


def test_unique_atoms():
    """ test automol.graph._func_group.chem_unique_atoms_of_type
    """

    uni_idxs1 = automol.graph.atom_equivalence_class_reps(SPCS_CHECKS_SMI['C4H10'], symb='H')
    assert uni_idxs1 == {4, 10}

    uni_idxs2 = automol.graph.atom_equivalence_class_reps(SPCS_CHECKS_SMI['C4H10'], symb='C')
    assert uni_idxs2 == {0, 2}


if __name__ == '__main__':
    __sites()
    test_functional_group_dct()
    test_unique_atoms()
    test_species_types()