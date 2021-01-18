"""
  Tests the functional group
"""

import automol


# Graphs used to inspect for functional groups
C2H6_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CC')))

C4H10_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CCCC')))

C2H5OH_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CCO')))

C2H5CL_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CCCl')))

C2H5SH_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CCS')))

CH3OCH3_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('COC')))

CYC_ETHER_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CC1CO1')))

C2H5OOH_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CCOO')))

C2H5OO_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CCO[O]')))

CCOCCO_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CC(=O)CC=O')))

C2H5CO_OH_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CC(=O)O')))

CCOOC_GRA = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CC(=O)OC')))

GRAP = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('C[N+](=O)[O-]')))

GRAP = automol.geom.graph(
    automol.inchi.geometry(
        automol.smiles.inchi('CC(=O)N')))

INI_FGRP_DCT = {
    automol.graph.FunctionalGroup.PEROXY: (),
    automol.graph.FunctionalGroup.HYDROPEROXY: (),
    automol.graph.FunctionalGroup.ETHER: (),
    automol.graph.FunctionalGroup.EPOXIDE: (),
    automol.graph.FunctionalGroup.CARBOX_ACID: (),
    automol.graph.FunctionalGroup.ESTER: (),
    automol.graph.FunctionalGroup.ALCOHOL: (),
    automol.graph.FunctionalGroup.ALDEHYDE: (),
    automol.graph.FunctionalGroup.KETONE: (),
    automol.graph.FunctionalGroup.AMIDE: (),
    automol.graph.FunctionalGroup.NITRO: (),
    automol.graph.FunctionalGroup.HALIDE: (),
    automol.graph.FunctionalGroup.THIOL: ()
}


def test_functional_group_dct():
    """ test automol.graph.functional_group_dct
    """

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.ALCOHOL: ((1, 2, 8),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5OH_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.HALIDE: ((1, 2),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5CL_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.THIOL: ((1, 2, 8),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5SH_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.ETHER: ((0, 2, 1),)
    })
    fgrps = automol.graph.functional_group_dct(CH3OCH3_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.EPOXIDE: ((1, 3, 2),)
    })
    fgrps = automol.graph.functional_group_dct(CYC_ETHER_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.HYDROPEROXY: ((1, 3, 2, 9),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5OOH_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.PEROXY: ((1, 3, 2),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5OO_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.ALDEHYDE: ((2, 4),),
        automol.graph.FunctionalGroup.KETONE: ((3, 5),)
    })
    fgrps = automol.graph.functional_group_dct(CCOCCO_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.CARBOX_ACID: ((2, 1, 3, 7),)
    })
    fgrps = automol.graph.functional_group_dct(C2H5CO_OH_GRA)
    assert fgrps == ref_fgrps

    ref_fgrps = INI_FGRP_DCT.copy()
    ref_fgrps.update({
        automol.graph.FunctionalGroup.ESTER: ((3, 2, 4, 1),)
    })
    fgrps = automol.graph.functional_group_dct(CCOOC_GRA)
    assert fgrps == ref_fgrps


def test_species_types():
    """ test automol.graph._func_group.hydrocarbon_species
             automol.graph._func_group.radical_species
    """

    assert automol.graph.hydrocarbon_species(C2H6_GRA)
    assert not automol.graph.hydrocarbon_species(C2H5OH_GRA)

    assert automol.graph.radical_species(C2H5OO_GRA)
    assert not automol.graph.radical_species(C2H5OH_GRA)


def test_unique_atoms():
    """ test automol.graph._func_group.chem_unique_atoms_of_type
    """

    uni_idxs1 = automol.graph.chem_unique_atoms_of_type(C4H10_GRA, 'H')
    assert uni_idxs1 == (4, 10)

    uni_idxs2 = automol.graph.chem_unique_atoms_of_type(C4H10_GRA, 'C')
    assert uni_idxs2 == (0, 2)

test_functional_group_dct()
test_species_types()
test_unique_atoms()

