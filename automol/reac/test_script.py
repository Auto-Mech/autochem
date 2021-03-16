""" test script
"""
import automol

rct_smis = ['CCC[CH2]']
prd_smis = ['CC[CH]C']

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

geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

scan_name = automol.reac.scan_coordinate(zrxn, zma)
print(scan_name)
