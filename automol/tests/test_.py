""" test automol.reac
"""

import sys
import automol


def test__reac__hydrogen_abstraction():
    """ test hydrogen abstraction functionality
    """
    rct_smis = ['C', '[H]']
    prd_smis = ['[CH3]', '[HH]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=True)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    rxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)
    tsg = rxn.forward_ts_graph

    # lin_keys = sorted(
    #     automol.graph.dummy_atoms_neighbor_atom_key(tsg).values())
    # bnd_keys = automol.graph.rotational_bond_keys(tsg, lin_keys=lin_keys)
    # names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    # assert names == {'D11', 'D5', 'D8', 'D16'}
    # print(automol.zmat.string(zma, one_indexed=False))
    # print(names)

    # scan_name = automol.reac.scan_coordinate(rxn, zma)
    # const_names = automol.reac.constraint_coordinates(rxn, zma)
    # assert scan_name == 'R15'
    # assert const_names == ()
    # print(scan_name)
    # print(const_names)

    # # graph aligned to geometry keys
    # # (for getting rotational groups and symmetry numbers)
    # geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    # grxn = automol.reac.relabel_for_geometry(rxn)
    # gtsg = grxn.forward_ts_graph

    # # Check that the reaction object can be converted back, if needed
    # zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    # assert zrxn == rxn

    # lin_keys = sorted(gdummy_key_dct.keys())
    # gbnd_keys = automol.graph.rotational_bond_keys(gtsg, lin_keys=lin_keys)
    # assert len(gbnd_keys) == len(bnd_keys)

    # axes = sorted(map(sorted, gbnd_keys))
    # groups_lst = [automol.graph.rotational_groups(gtsg, *a) for a in axes]
    # sym_nums = [
    #     automol.graph.rotational_symmetry_number(gtsg, *a, lin_keys=lin_keys)
    #     for a in axes]
    # assert sym_nums == [3, 3, 3, 1]
    # for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
    #     print('axis:', axis)
    #     print('\tgroup 1:', groups[0])
    #     print('\tgroup 2:', groups[1])
    #     print('\tsymmetry number:', sym_num)


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
    test__reac__hydrogen_abstraction()
