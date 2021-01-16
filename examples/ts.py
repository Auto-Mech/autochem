""" Basic demo of distance geometry for transition states
"""
import sys
import automol

# 1. Choose reaction

#    a. hydrogen migration: (CH3)2[CH]CH2CH2O[O] => (CH3)2[C]CH2CH2O[OH]
# RCT_ICHS = list(map(automol.smiles.inchi, ['C1CCC1C(CC2)C2CO[O]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['C1CCC1[C](CC2)C2COO']))

#    b. beta-scission: CH3CH2CH2OO => CH3CH2[CH2] + OO
# RCT_ICHS = list(map(automol.smiles.inchi, ['CCCO[O]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['[O][O]', 'CC[CH2]']))

#    c. ring-forming scission: [CH2]CH2OOH => cCH2CH2O + OH
# RCT_ICHS = list(map(automol.smiles.inchi, ['[CH2]COO']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['C1CO1', '[OH]']))

#    d. elimination 1: CH3CH2OO => C2H4 + HOO
# RCT_ICHS = list(map(automol.smiles.inchi, ['CCO[O]']))
# PRD_ICHS = reversed(list(map(automol.smiles.inchi, ['C=C', 'O[O]'])))

#    d. elimination 2: CH3CH2CH3 => CH3CH3 + CH2
# RCT_ICHS = list(map(automol.smiles.inchi, ['CCC']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['CC', '[CH2]']))

#    e. hydrogen abstraction 1: HC(CH3)3 + OH => C(CH3)3 + H2O
# RCT_ICHS = list(map(automol.smiles.inchi, ['C(C)(C)C', '[OH]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['[C](C)(C)C', 'O']))

#    e. hydrogen abstraction 2: HC(CH3)3 + OH => C(CH3)3 + H2O
# RCT_ICHS = list(map(automol.smiles.inchi, ['C(C)(C)C#C', '[OH]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['[C](C)(C)C#C', 'O']))

#    f. addition: CH3CH2[CH2] + OO => CH3CH2CH2OO
# RCT_ICHS = list(map(automol.smiles.inchi, ['CC[CH2]', '[O][O]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['CCCO[O]']))

#    g. insertion 1: C2H4 + HOO => CH3CH2OO
# RCT_ICHS = reversed(list(map(automol.smiles.inchi, ['C=C', 'O[O]'])))
# PRD_ICHS = list(map(automol.smiles.inchi, ['CCO[O]']))

#    g. insertion 2: CH3CH3 + CH2 => CH3CH2CH3
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
RXNS = automol.reac.find(RCT_GRAS, PRD_GRAS)
RXN = RXNS[0]
# Sort and standardize keys (must go together)
RCT_GEOS, PRD_GEOS = RXN.standardize_keys_and_sort_geometries_(
    RCT_GEOS, PRD_GEOS)

# 4. Generate the geometry
GEO = automol.reac.ts_geometry(RXN, RCT_GEOS, log=True)

# 5. Generate the z-matrix
ZMA, ROW_KEYS, DUMMY_IDX_DCT = automol.reac.ts_zmatrix(RXN, GEO)

# 8. Print some stuff
print("Forward TS graph (lined up with geometry):")
TSG = RXN.forward_ts_graph
print(automol.graph.string(TSG, one_indexed=False))

print("Forward TS graph (lined up with zmatrix):")
RXN.insert_dummy_atoms_(DUMMY_IDX_DCT)
RXN.relabel_(dict(map(reversed, enumerate(ROW_KEYS))))
TSG = RXN.forward_ts_graph
print(automol.graph.string(TSG, one_indexed=False))

print("Reactant 1 geometry:")
print(automol.geom.string(RCT_GEOS[0]))
print()

print("Cleaned up geometry:")
print(automol.geom.string(GEO))
print()

print("Final z-matrix (one-indexed for comparison):")
print(automol.zmat.string(ZMA, one_indexed=False))

print("Row keys for z-matrix:")
print(ROW_KEYS)

print("Dummy indices for z-matrix:")
print(DUMMY_IDX_DCT)

# 9. Check the connectivity
GRA = automol.geom.graph(automol.geom.join(*RCT_GEOS) if len(RCT_GEOS) > 1
                         else RCT_GEOS[0])
GRA2 = automol.geom.graph(GEO)
print("Is the graph consistent?", 'Yes' if GRA == GRA2 else 'No')
print()

sys.exit()
