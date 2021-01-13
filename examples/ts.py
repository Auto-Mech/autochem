""" Basic demo of distance geometry for transition states
"""
import sys
# import numpy
import itertools
import automol

# 1. Choose reaction

#    a. hydrogen migration: (CH3)2[CH]CH2CH2O[O] => (CH3)2[C]CH2CH2O[OH]
# RCT_ICHS = list(map(automol.smiles.inchi, ['C(C)(C)CCO[O]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['[C](C)(C)CCOO']))

#    b. beta-scission: CH3CH2CH2OO => CH3CH2[CH2] + OO
# RCT_ICHS = list(map(automol.smiles.inchi, ['CCCO[O]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['[O][O]', 'CC[CH2]']))

#    c. ring-forming scission: [CH2]CH2OOH => cCH2CH2O + OH
# RCT_ICHS = list(map(automol.smiles.inchi, ['[CH2]COO']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['C1CO1', '[OH]']))

#    d. elimination: CH3CH2OO => C2H4 + HOO
# RCT_ICHS = list(map(automol.smiles.inchi, ['CCO[O]']))
# PRD_ICHS = reversed(list(map(automol.smiles.inchi, ['C=C', 'O[O]'])))

#    e. hydrogen abstraction: HC(CH3)3 + OH => C(CH3)3 + H2O
# RCT_ICHS = list(map(automol.smiles.inchi, ['C(C)(C)C', '[OH]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['[C](C)(C)C', 'O']))

#    f. addition: CH3CH2[CH2] + OO => CH3CH2CH2OO
# RCT_ICHS = list(map(automol.smiles.inchi, ['CC[CH2]', '[O][O]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['CCCO[O]']))

#    g. insertion 1: CH3CH3 + CH2 => C2H4 + HOO
# RCT_ICHS = reversed(list(map(automol.smiles.inchi, ['C=C', 'O[O]'])))
# PRD_ICHS = list(map(automol.smiles.inchi, ['CCO[O]']))

#    g. insertion 2: CH3CH3 + CH2 => C2H4 + HOO
# RCT_ICHS = list(map(automol.smiles.inchi, ['CC', '[CH2]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['CCC']))

#    h. substitution: H2O2 + H => H2O + OH
RCT_ICHS = list(map(automol.smiles.inchi, ['CO', '[CH2]C']))
PRD_ICHS = list(map(automol.smiles.inchi, ['CCC', '[OH]']))

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
RXNS = automol.reac.classify(RCT_GRAS, PRD_GRAS)
RXN = RXNS[0]
# Sort and standardize keys (must go together)
RCT_GEOS, PRD_GEOS = RXN.standardize_keys_and_sort_geometries_(
    RCT_GEOS, PRD_GEOS)

# 4. Generate embedding info
XMAT, LMAT, UMAT, CHI_DCT, PLA_DCT = automol.reac.ts_embedding_info(RXN,
                                                                    RCT_GEOS)

# 7. Optimize the coordinates to satisfy the appropriate constraints
XMAT, CONV = automol.embed.cleaned_up_coordinates(
    XMAT, LMAT, UMAT, chi_dct=CHI_DCT, pla_dct=PLA_DCT, max_dist_err=2e-1,
    dim4=False, log=True)

SYMS = list(itertools.chain(*map(automol.geom.symbols, RCT_GEOS)))
GEO = automol.embed.geometry_from_coordinates(XMAT, SYMS)

# 8. Print some stuff
print("Forward TS graph:")
TSG = RXN.forward_ts_graph
print(automol.graph.string(TSG, one_indexed=False))

print("Reactant 1 geometry:")
print(automol.geom.string(RCT_GEOS[0]))
print()

print("Cleaned up geometry:")
print(automol.geom.string(GEO))
print()

# 9. Check the connectivity
GRA = automol.graph.ts.reactants_graph(TSG)
GRA2 = automol.geom.graph(GEO)
print("Is the graph consistent?", 'Yes' if GRA == GRA2 else 'No')
print("Is the geometry converged?", 'Yes.' if CONV else 'No.')
print()

# 10. Generate a z-matrix
FRM_BND_KEYS = automol.graph.ts.forming_bond_keys(RXN.forward_ts_graph)
BRK_BND_KEYS = automol.graph.ts.breaking_bond_keys(RXN.forward_ts_graph)
TS_BNDS = list(FRM_BND_KEYS | BRK_BND_KEYS)
ZMA = automol.convert.geom.zmatrix_x2z(GEO, TS_BNDS)
print(automol.zmat.string(ZMA))

sys.exit()
