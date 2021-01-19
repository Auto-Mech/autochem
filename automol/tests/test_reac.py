""" test automol.reac
"""
import automol


def test__species__demo():
    """ doesn't really belong here, but demonstrates equivalent functionality
    for species
    """
    ich = automol.smiles.inchi('CC#CC#CCCCC#CC')
    geo = automol.inchi.geometry(ich)
    gra = automol.geom.graph(geo)

    # graph aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.convert.geom.zmatrix(geo)
    gra = automol.graph.relabel_for_zmatrix(gra, zma_keys, dummy_key_dct)

    lin_keys = sorted(automol.graph.dummy_atom_neighbor_keys(gra).values())
    bnd_keys = automol.graph.rotational_bond_keys(gra, lin_keys=lin_keys)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D9', 'D12', 'D15', 'D26'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    ggra = automol.graph.relabel_for_geometry(gra)

    # Check that the geometry graph can be converted back, if needed
    zgra = automol.graph.insert_dummy_atoms(ggra, gdummy_key_dct)
    assert zgra == gra

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(ggra, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.graph.rotational_groups(ggra, *a) for a in axes]
    sym_nums = [
        automol.graph.rotational_symmetry_number(ggra, *a, lin_keys=lin_keys)
        for a in axes]
    assert sym_nums == [3, 1, 1, 3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__hydrogen_migration():
    """ test hydrogen migration functionality
    """
    rct_smis = ['CCCO[O]']
    prd_smis = ['C[CH]COO']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=True)
    print(automol.geom.string(geo))

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    rxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)
    tsg = rxn.forward_ts_graph

    lin_keys = sorted(automol.graph.dummy_atom_neighbor_keys(tsg).values())
    bnd_keys = automol.graph.rotational_bond_keys(tsg, lin_keys=lin_keys)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D9'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_bnd_key = automol.reac.scan_coordinate(rxn)
    scan_name = automol.zmat.distance_coordinate_name(zma, *scan_bnd_key)
    assert scan_name == 'R5'
    print(scan_name)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(rxn)
    gtsg = grxn.forward_ts_graph

    # Check that the reaction object can be converted back, if needed
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == rxn

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(gtsg, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.graph.rotational_groups(gtsg, *a) for a in axes]
    sym_nums = [
        automol.graph.rotational_symmetry_number(gtsg, *a, lin_keys=lin_keys)
        for a in axes]
    assert sym_nums == [3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__beta_scission():
    """ test beta scission functionality
    """
    rct_smis = ['CCCO[O]']
    prd_smis = ['[O][O]', 'CC[CH2]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=True)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    rxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)
    tsg = rxn.forward_ts_graph

    lin_keys = sorted(automol.graph.dummy_atom_neighbor_keys(tsg).values())
    bnd_keys = automol.graph.rotational_bond_keys(tsg, lin_keys=lin_keys)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D8', 'D11', 'D5'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_bnd_key = automol.reac.scan_coordinate(rxn)
    scan_name = automol.zmat.distance_coordinate_name(zma, *scan_bnd_key)
    assert scan_name == 'R8'
    print(scan_name)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(rxn)
    gtsg = grxn.forward_ts_graph

    # Check that the reaction object can be converted back, if needed
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == rxn

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(gtsg, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.graph.rotational_groups(gtsg, *a) for a in axes]
    sym_nums = [
        automol.graph.rotational_symmetry_number(gtsg, *a, lin_keys=lin_keys)
        for a in axes]
    assert sym_nums == [3, 1, 1]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__ring_forming_scission():
    """ test ring-forming scission functionality
    """
    rct_smis = ['[CH2]COO']
    prd_smis = ['C1CO1', '[OH]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=True)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    rxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)
    tsg = rxn.forward_ts_graph

    lin_keys = sorted(automol.graph.dummy_atom_neighbor_keys(tsg).values())
    bnd_keys = automol.graph.rotational_bond_keys(tsg, lin_keys=lin_keys)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D8'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_bnd_key = automol.reac.scan_coordinate(rxn)
    scan_name = automol.zmat.distance_coordinate_name(zma, *scan_bnd_key)
    assert scan_name == 'R7'
    print(scan_name)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(rxn)
    gtsg = grxn.forward_ts_graph

    # Check that the reaction object can be converted back, if needed
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == rxn

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(gtsg, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.graph.rotational_groups(gtsg, *a) for a in axes]
    sym_nums = [
        automol.graph.rotational_symmetry_number(gtsg, *a, lin_keys=lin_keys)
        for a in axes]
    assert sym_nums == [1]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__elimination():
    """ test elimination functionality
    """
    rct_smis = ['CCCO[O]']
    prd_smis = ['CC=C', 'O[O]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=True)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    rxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)
    tsg = rxn.forward_ts_graph

    lin_keys = sorted(automol.graph.dummy_atom_neighbor_keys(tsg).values())
    bnd_keys = automol.graph.rotational_bond_keys(tsg, lin_keys=lin_keys)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D9'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_bnd_key = automol.reac.scan_coordinate(rxn)
    scan_name = automol.zmat.distance_coordinate_name(zma, *scan_bnd_key)
    assert scan_name == 'R2'
    print(scan_name)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(rxn)
    gtsg = grxn.forward_ts_graph

    # Check that the reaction object can be converted back, if needed
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == rxn

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(gtsg, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.graph.rotational_groups(gtsg, *a) for a in axes]
    sym_nums = [
        automol.graph.rotational_symmetry_number(gtsg, *a, lin_keys=lin_keys)
        for a in axes]
    assert sym_nums == [3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__hydrogen_abstraction():
    """ test hydrogen abstraction functionality
    """
    rct_smis = ['C(C)(C)C', '[OH]']
    prd_smis = ['[C](C)(C)C', 'O']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=True)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    rxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)
    tsg = rxn.forward_ts_graph

    lin_keys = sorted(automol.graph.dummy_atom_neighbor_keys(tsg).values())
    bnd_keys = automol.graph.rotational_bond_keys(tsg, lin_keys=lin_keys)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D11', 'D5', 'D8', 'D16'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_bnd_key = automol.reac.scan_coordinate(rxn)
    scan_name = automol.zmat.distance_coordinate_name(zma, *scan_bnd_key)
    assert scan_name == 'R15'
    print(scan_name)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(rxn)
    gtsg = grxn.forward_ts_graph

    # Check that the reaction object can be converted back, if needed
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == rxn

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(gtsg, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.graph.rotational_groups(gtsg, *a) for a in axes]
    sym_nums = [
        automol.graph.rotational_symmetry_number(gtsg, *a, lin_keys=lin_keys)
        for a in axes]
    assert sym_nums == [3, 3, 3, 1]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__addition():
    """ test addition functionality
    """
    rct_smis = ['CC[CH2]', '[O][O]']
    prd_smis = ['CCCO[O]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=True)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    rxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)
    tsg = rxn.forward_ts_graph

    lin_keys = sorted(automol.graph.dummy_atom_neighbor_keys(tsg).values())
    bnd_keys = automol.graph.rotational_bond_keys(tsg, lin_keys=lin_keys)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D11', 'D8', 'D5'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_bnd_key = automol.reac.scan_coordinate(rxn)
    scan_name = automol.zmat.distance_coordinate_name(zma, *scan_bnd_key)
    assert scan_name == 'R10'
    print(scan_name)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(rxn)
    gtsg = grxn.forward_ts_graph

    # Check that the reaction object can be converted back, if needed
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == rxn

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(gtsg, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.graph.rotational_groups(gtsg, *a) for a in axes]
    sym_nums = [
        automol.graph.rotational_symmetry_number(gtsg, *a, lin_keys=lin_keys)
        for a in axes]
    assert sym_nums == [3, 1, 1]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__insertion():
    """ test insertion functionality
    """
    rct_smis = ['CC=C', 'O[O]']
    prd_smis = ['CCCO[O]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=True)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    rxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)
    tsg = rxn.forward_ts_graph

    lin_keys = sorted(automol.graph.dummy_atom_neighbor_keys(tsg).values())
    bnd_keys = automol.graph.rotational_bond_keys(tsg, lin_keys=lin_keys)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D9'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_bnd_key = automol.reac.scan_coordinate(rxn)
    scan_name = automol.zmat.distance_coordinate_name(zma, *scan_bnd_key)
    assert scan_name == 'R3'
    print(scan_name)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(rxn)
    gtsg = grxn.forward_ts_graph

    # Check that the reaction object can be converted back, if needed
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == rxn

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(gtsg, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.graph.rotational_groups(gtsg, *a) for a in axes]
    sym_nums = [
        automol.graph.rotational_symmetry_number(gtsg, *a, lin_keys=lin_keys)
        for a in axes]
    assert sym_nums == [3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__substitution():
    """ test substitution functionality
    """
    rct_smis = ['CO', '[CH2]C']
    prd_smis = ['CCC', '[OH]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=True)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    rxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)
    tsg = rxn.forward_ts_graph

    lin_keys = sorted(automol.graph.dummy_atom_neighbor_keys(tsg).values())
    bnd_keys = automol.graph.rotational_bond_keys(tsg, lin_keys=lin_keys)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D3', 'D8', 'D11'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_bnd_key = automol.reac.scan_coordinate(rxn)
    scan_name = automol.zmat.distance_coordinate_name(zma, *scan_bnd_key)
    assert scan_name == 'R7'
    print(scan_name)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(rxn)
    gtsg = grxn.forward_ts_graph

    # Check that the reaction object can be converted back, if needed
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == rxn

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(gtsg, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.graph.rotational_groups(gtsg, *a) for a in axes]
    sym_nums = [
        automol.graph.rotational_symmetry_number(gtsg, *a, lin_keys=lin_keys)
        for a in axes]
    assert sym_nums == [1, 1, 3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def _from_smiles(rct_smis, prd_smis):
    rct_ichs = list(map(automol.smiles.inchi, rct_smis))
    prd_ichs = list(map(automol.smiles.inchi, prd_smis))

    rct_geos = list(map(automol.inchi.geometry, rct_ichs))
    prd_geos = list(map(automol.inchi.geometry, prd_ichs))

    rct_gras = list(map(automol.graph.without_stereo_parities,
                        map(automol.geom.graph, rct_geos)))
    prd_gras = list(map(automol.graph.without_stereo_parities,
                        map(automol.geom.graph, prd_geos)))

    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = automol.reac.find(rct_gras, prd_gras)
    rxn = rxns[0]

    rxn, rct_geos, prd_geos = (
        automol.reac.standard_keys_with_sorted_geometries(
            rxn, rct_geos, prd_geos))
    return rxn, rct_geos, prd_geos


if __name__ == '__main__':
    # test__species__demo()
    # test__reac__hydrogen_migration()
    # test__reac__beta_scission()
    # test__reac__ring_forming_scission()
    # test__reac__elimination()
    # test__reac__hydrogen_abstraction()
    # test__reac__addition()
    # test__reac__insertion()
    test__reac__substitution()
