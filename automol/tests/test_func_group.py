"""
  Tests the functional group
"""

import automol


# Graphs used to inspect for functional groups
C2H6_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CC')))

C4H10_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CCCC')))

C4H6_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('C=CC=C')))

C4H2_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('C#CC#C')))

C2H5OH_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CCO')))

C2H5CL_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CCCl')))

C2H5SH_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CCS')))

CH3OCH3_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('COC')))

CYC_ETHER_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CC1CO1')))

C2H5OOH_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CCOO')))

C2H5OO_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CCO[O]')))

CCOCCO_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CC(=O)CC=O')))

C2H5CO_OH_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CC(=O)O')))

CCOOC_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('CC(=O)OC')))

C6H5_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('c1ccccc1')))

AMN_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('Cc1ccc2ccccc2c1')))

INDENYL_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('[CH]1C=Cc2ccccc12')))

C12H8_GRA = automol.geom.graph(
    automol.chi.geometry(
        automol.smiles.chi('c1cc2cccc3C=Cc(c1)c23')))

# GRAP = automol.geom.graph(
#     automol.chi.geometry(
#         automol.smiles.chi('C[N+](=O)[O-]')))
# GRAP = automol.geom.graph(
#     automol.chi.geometry(
#         automol.smiles.chi('CC(=O)N')))

INI_FGRP_DCT = {
    automol.graph.FunctionalGroup.ALKENE: (),
    automol.graph.FunctionalGroup.ALKOXY: (),
    automol.graph.FunctionalGroup.PEROXY: (),
    automol.graph.FunctionalGroup.HYDROPEROXY: (),
    automol.graph.FunctionalGroup.ETHER: (),
    automol.graph.FunctionalGroup.CYCLIC_ETHER: (),
    automol.graph.FunctionalGroup.CARBOX_ACID: (),
    automol.graph.FunctionalGroup.ESTER: (),
    automol.graph.FunctionalGroup.ALCOHOL: (),
    automol.graph.FunctionalGroup.ALDEHYDE: (),
    automol.graph.FunctionalGroup.KETONE: (),
    automol.graph.FunctionalGroup.AMIDE: (),
    automol.graph.FunctionalGroup.NITRO: (),
    automol.graph.FunctionalGroup.HALIDE: (),
    automol.graph.FunctionalGroup.THIOL: (),
    automol.graph.FunctionalGroup.METHYL: (),
    automol.graph.FunctionalGroup.PHENYL: (),
    automol.graph.FunctionalGroup.AMINE: (),
    automol.graph.FunctionalGroup.AROMATIC: (),
}


def test_functional_group_dct():
    """ test automol.graph.functional_group_dct
    """
    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.ALKENE: ((6, 7), (10, 5), (1, 3)), 
        automol.graph.FunctionalGroup.AROMATIC: ((0, 2, 8, 11, 9, 4), (1, 3, 8, 11, 10, 5)),
        automol.graph.FunctionalGroup.PHENYL: ((8, 2, 0, 4, 9, 11),),
     })
    fgrps = automol.graph.functional_group_dct(C12H8_GRA)
    assert fgrps == ref_fgrps
    
    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.AROMATIC: ((0, 1, 3, 5, 4, 2),),
        automol.graph.FunctionalGroup.PHENYL: ((4, 5, 0, 2, 1, 3),),
     })
    fgrps = automol.graph.functional_group_dct(C6H5_GRA)
    
    assert fgrps == ref_fgrps
    
    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.METHYL: ((0, 12, 11, 13),),
        automol.graph.FunctionalGroup.PHENYL: ((8, 7, 9, 10, 5, 6), (9, 10, 2, 4, 1, 3)),
        automol.graph.FunctionalGroup.AROMATIC: ((5, 6, 9, 10, 7, 8), (1, 2, 4, 10, 9, 3)),
     })
    fgrps = automol.graph.functional_group_dct(AMN_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.ALKENE: ((2, 6),), 
        automol.graph.FunctionalGroup.PHENYL: ((1, 4, 8, 7, 0, 3),),
        automol.graph.FunctionalGroup.AROMATIC: ((0, 1, 4, 8, 7, 3),),
     })
    fgrps = automol.graph.functional_group_dct(INDENYL_GRA)

    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.ALCOHOL: ((1, 2, 8),),
        automol.graph.FunctionalGroup.METHYL: ((0, 3, 4, 5),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5OH_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.HALIDE: ((1, 2),),
        automol.graph.FunctionalGroup.METHYL: ((0, 3, 4, 5),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5CL_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.THIOL: ((1, 2, 8),),
        automol.graph.FunctionalGroup.METHYL: ((0, 3, 4, 5),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5SH_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.ETHER: ((0, 2, 1),),
        automol.graph.FunctionalGroup.METHYL: ((1, 7, 8, 6), (0, 3, 4, 5))
    })
    fgrps = automol.graph.functional_group_dct(CH3OCH3_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.CYCLIC_ETHER: ((1, 3, 2),),
        automol.graph.FunctionalGroup.METHYL: ((0, 6, 4, 5),),
    })
    fgrps = automol.graph.functional_group_dct(CYC_ETHER_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.HYDROPEROXY: ((1, 3, 2, 9),),
        automol.graph.FunctionalGroup.METHYL: ((0, 6, 4, 5),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5OOH_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.PEROXY: ((1, 3, 2),),
        automol.graph.FunctionalGroup.METHYL: ((0, 6, 4, 5),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5OO_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.ALDEHYDE: ((2, 4),),
        automol.graph.FunctionalGroup.KETONE: ((3, 5),),
        automol.graph.FunctionalGroup.METHYL: ((0, 6, 7, 8),)
    })
    fgrps = automol.graph.functional_group_dct(CCOCCO_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.CARBOX_ACID: ((2, 1, 3, 7),),
        automol.graph.FunctionalGroup.METHYL: ((0, 6, 4, 5),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5CO_OH_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.ESTER: ((3, 2, 4, 1),),
        automol.graph.FunctionalGroup.METHYL: ((1, 9, 8, 10), (0, 6, 5, 7))
    })
    fgrps = automol.graph.functional_group_dct(CCOOC_GRA)
    assert fgrps == ref_fgrps


def __sites():
    """ test automol.graph._func_group.alkene_sites
        test automol.graph._func_group.alkene_sites
    """

    assert automol.graph.alkene_sites(C4H6_GRA) == ((0, 2), (1, 3))
    assert automol.graph.alkyne_sites(C4H2_GRA) == ((0, 2), (1, 3))
    assert not automol.graph.alkene_sites(C4H2_GRA)
    assert not automol.graph.alkyne_sites(C4H6_GRA)


def test_species_types():
    """ test automol.graph._func_group.hydrocarbon_species
             automol.graph._func_group.radical_species
    """

    assert automol.graph.is_hydrocarbon_species(C2H6_GRA)
    assert not automol.graph.is_hydrocarbon_species(C2H5OH_GRA)

    assert automol.graph.is_radical_species(C2H5OO_GRA)
    assert not automol.graph.is_radical_species(C2H5OH_GRA)


def test_unique_atoms():
    """ test automol.graph._func_group.chem_unique_atoms_of_type
    """

    uni_idxs1 = automol.graph.atom_equivalence_class_reps(C4H10_GRA, symb='H')
    print(uni_idxs1)
    assert uni_idxs1 == {4, 10}

    uni_idxs2 = automol.graph.atom_equivalence_class_reps(C4H10_GRA, symb='C')
    print(uni_idxs2)
    assert uni_idxs2 == {0, 2}


if __name__ == '__main__':
    test_functional_group_dct()
    test_unique_atoms()
