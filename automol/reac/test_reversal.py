""" test script
"""
import automol

rct_smis = ['CC', '[O][O]']
prd_smis = ['C[CH2]', 'O[O]']

rct_ichs = list(map(automol.smiles.chi, rct_smis))
prd_ichs = list(map(automol.smiles.chi, prd_smis))

rct_geos = list(map(automol.chi.geometry, rct_ichs))
prd_geos = list(map(automol.chi.geometry, prd_ichs))

rct_gras = list(map(automol.graph.without_stereo,
                    map(automol.geom.graph, rct_geos)))
prd_gras = list(map(automol.graph.without_stereo,
                    map(automol.geom.graph, prd_geos)))

rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

rxns = automol.reac.find(rct_gras, prd_gras)
rxn = rxns[-1]
print(automol.reac.string(rxn))

rxn, rct_geos, prd_geos = (
    automol.reac.standard_keys_with_sorted_geometries(
        rxn, rct_geos, prd_geos))

# First, generate a TS z-matrix and Reaction object for the forward direction
#   1. Generate the TS geometry
geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)
print("1. Geometry:")
print(automol.geom.string(geo))

#   2. Generate z-matrix
zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
print("2. Z-matrix:")
print(automol.zmat.string(zma, one_indexed=False))
zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

#   3. Identify the scan coordinate
scan_name = automol.reac.scan_coordinate(zrxn, zma)
print("3. Scan coordinate:")
print(scan_name)

#   4. Identify the torsional coordinates
bnd_keys = automol.reac.rotational_bond_keys(zrxn)
print("4. Rotational bonds and torsion coordinates:")
print(bnd_keys)
names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
print(names)


# Now, do exactly the same thing for the reverse direction
#   0. Reverse the direction of the reaction
print("\nNow in reverse...")
rxn_rev = automol.reac.reverse(rxn)
print(automol.reac.string(rxn_rev))

#   1. Generate the TS geometry (use the product geoms since this is reversed)
geo = automol.reac.ts_geometry(rxn_rev, prd_geos, log=False)
print("1. Geometry:")
print(automol.geom.string(geo))

#   2. Generate z-matrix
zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn_rev, geo)
print("2. Z-matrix:")
print(automol.zmat.string(zma, one_indexed=False))
zrxn_rev = automol.reac.relabel_for_zmatrix(rxn_rev, zma_keys, dummy_key_dct)

#   3. Identify the scan coordinate
scan_name = automol.reac.scan_coordinate(zrxn_rev, zma)
print("3. Scan coordinate:")
print(scan_name)

#   4. Identify the torsional coordinates
bnd_keys = automol.reac.rotational_bond_keys(zrxn_rev)
print("4. Rotational bonds and torsion coordinates:")
print(bnd_keys)
names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
print(names)
