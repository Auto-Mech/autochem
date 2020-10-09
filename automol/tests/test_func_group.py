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


def test_functional_group_dct():
    """ test automol.graph.functional_group_dct
    """

    print('\nethanol')
    ref_fgrps = {
        automol.graph.FUNC_GROUP.ALCOHOL: ((1, 2, 8),)
    }
    fgrps = automol.graph.functional_group_dct(C2H5OH_GRA)
    assert fgrps == ref_fgrps

    print('\nethyl halide')
    ref_fgrps = {
        automol.graph.FUNC_GROUP.HALIDE: ((1, 2),)
    }
    fgrps = automol.graph.functional_group_dct(C2H5CL_GRA)
    assert fgrps == ref_fgrps

    print('\nethyl thiol')
    ref_fgrps = {
        automol.graph.FUNC_GROUP.THIOL: ((1, 2, 8),)
    }
    fgrps = automol.graph.functional_group_dct(C2H5SH_GRA)
    assert fgrps == ref_fgrps

    print('\ndimethyl ether')
    ref_fgrps = {
        automol.graph.FUNC_GROUP.ETHER: ((0, 2, 1),)
    }
    fgrps = automol.graph.functional_group_dct(CH3OCH3_GRA)
    assert fgrps == ref_fgrps

    print('\n1,2-epoxypropane')
    ref_fgrps = {
        automol.graph.FUNC_GROUP.EPOXIDE: ((1, 3, 2),)
    }
    fgrps = automol.graph.functional_group_dct(CYC_ETHER_GRA)
    assert fgrps == ref_fgrps

    print('\nperoxide')
    ref_fgrps = {
        automol.graph.FUNC_GROUP.HYDROPEROXY: ((1, 3, 2, 9),)
    }
    fgrps = automol.graph.functional_group_dct(C2H5OOH_GRA)
    assert fgrps == ref_fgrps

    print('\nperoxy')
    ref_fgrps = {
        automol.graph.FUNC_GROUP.PEROXY: ((1, 3, 2),)
    }
    fgrps = automol.graph.functional_group_dct(C2H5OO_GRA)
    assert fgrps == ref_fgrps

    print('\nald and ketone')
    ref_fgrps = {
        automol.graph.FUNC_GROUP.ALDEHYDE: ((2, 4),),
        automol.graph.FUNC_GROUP.KETONE: ((3, 5),)
    }
    fgrps = automol.graph.functional_group_dct(CCOCCO_GRA)
    assert fgrps == ref_fgrps

    print('\ncarbox')
    ref_fgrps = {
        automol.graph.FUNC_GROUP.CARBOX_ACID: ((2, 1, 3, 7),)
    }
    fgrps = automol.graph.functional_group_dct(C2H5CO_OH_GRA)
    assert fgrps == ref_fgrps

    print('\nester')
    ref_fgrps = {
        automol.graph.FUNC_GROUP.ESTER: ((3, 2, 4, 1),)
    }
    fgrps = automol.graph.functional_group_dct(CCOOC_GRA)
    assert fgrps == ref_fgrps

    # print('\nnitromethane')
    # ref_fgrps = {
    # }
    # FGRPS = functional_group_dct(GRAP)
    # assert fgrps == ref_fgrps

    # print('\namide')
    # ref_fgrps = {
    # }
    # FGRPS = functional_group_dct(GRAP)
    # assert fgrps == ref_fgrps


def test_species_types():
    """ test automol.graph._func_group.hydrocarbon_species
             automol.graph._func_group.radical_species
    """

    print('\n is hydrocarbon')
    assert automol.graph.hydrocarbon_species(C2H6_GRA)
    assert not automol.graph.hydrocarbon_species(C2H5OH_GRA)

    print('\n is radical')
    assert automol.graph.radical_species(C2H5OH_GRA)
    assert not automol.graph.radical_species(C2H5OO_GRA)


def test_unique_atoms():
    """ test automol.graph._func_group.chem_unique_atoms_of_type
    """
    
    print(automol.geom.string(
        automol.inchi.geometry(
            automol.smiles.inchi('CCCC'))))

    uni_idxs = automol.graph.chem_unique_atoms_of_type(C4H10_GRA, 'H')
    print('uni ids', uni_idxs)
    uni_idxs = automol.graph.chem_unique_atoms_of_type(C4H10_GRA, 'C')
    print('uni ids', uni_idxs)


if __name__ == '__main__':
    # test_functional_group_dct()
    # test_species_types()
    test_unique_atoms()
