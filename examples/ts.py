""" Basic demo of distance geometry for transition states
"""
import numpy
import automol

# 1. Choose reaction

#    a. bimolecular: HC(CH3)3 + OH => C(CH3)3 + H2O
# RCT_ICHS = list(map(automol.smiles.inchi, ['C(C)(C)C', '[OH]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['[C](C)(C)C', 'O']))
# UNIMOL = False

#    b. unimolecular: (CH3)2[CH]CH2CH2O[O] => (CH3)2[C]CH2CH2O[OH]
RCT_ICHS = list(map(automol.smiles.inchi, ['C(C)(C)CCO[O]']))
PRD_ICHS = list(map(automol.smiles.inchi, ['[C](C)(C)CCOO']))
UNIMOL = True

# 2. Generate reactant/product graphs
RCT_GEOS = list(map(automol.inchi.geometry, RCT_ICHS))
PRD_GEOS = list(map(automol.inchi.geometry, PRD_ICHS))

RCT_GRAS = list(map(automol.graph.without_stereo_parities,
                    map(automol.geom.graph, RCT_GEOS)))
PRD_GRAS = list(map(automol.graph.without_stereo_parities,
                    map(automol.geom.graph, PRD_GEOS)))

RCT_GRAS, _ = automol.graph.standard_keys_for_sequence(RCT_GRAS)
PRD_GRAS, _ = automol.graph.standard_keys_for_sequence(PRD_GRAS)

# 3. Classify reaction and get bonds broken/formed
TRAS, _, _ = automol.graph.reac.classify(RCT_GRAS, PRD_GRAS)
TRA = TRAS[0]
FRM_BND_KEYS = list(automol.graph.trans.formed_bond_keys(TRA))
BRK_BND_KEYS = list(automol.graph.trans.broken_bond_keys(TRA))
print("Bonds broken:", BRK_BND_KEYS)
print("Bonds formed:", FRM_BND_KEYS)

# 4. Generate reactant graph and sorted list of atom keys
GRA = automol.graph.union_from_sequence(RCT_GRAS)
KEYS = sorted(automol.graph.atom_keys(GRA))
SYMS = list(map(automol.graph.atom_symbols(GRA).__getitem__, KEYS))

# 5. Generate bounds matrices
FRM_BNDS_DCT = {bnd: (1.6, 1.6) for bnd in FRM_BND_KEYS}
LMAT, UMAT = automol.graph.embed.ts_distance_bounds_matrices(
        GRA, KEYS, FRM_BNDS_DCT, rct_geos=RCT_GEOS, relax_torsions=UNIMOL)
print("Lower bounds matrix:")
print(numpy.round(LMAT, 1))
print("Upper bounds matrix:")
print(numpy.round(UMAT, 1))
print()

# 6. Sample a geometry from the bounds matrices
XMAT = automol.embed.sample_raw_distance_coordinates(LMAT, UMAT, dim4=True)
GEO_INIT = automol.embed.geometry_from_coordinates(XMAT, SYMS)

# 7. Clean up the sample's coordinates
XMAT, CONV = automol.embed.cleaned_up_coordinates(XMAT, LMAT, UMAT)
GEO = automol.embed.geometry_from_coordinates(XMAT, SYMS)

# 8. Print geometries
print("Reactant 1 geometry:")
print(automol.geom.string(RCT_GEOS[0]))
print()

print("Sample geometry:")
print(automol.geom.string(GEO_INIT))
print()

print("Cleaned up geometry:")
print(automol.geom.string(GEO))
print()

# 9. Check the connectivity
GRA = automol.graph.without_stereo_parities(GRA)
GRA2 = automol.geom.connectivity_graph(GEO)
print("Is the connectivity consistent?", 'Yes' if GRA == GRA2 else 'No')

# 10. Generate a z-matrix
ZMA = automol.geom.zmatrix(GEO, FRM_BND_KEYS+BRK_BND_KEYS)
print(automol.zmatrix.string(ZMA))
