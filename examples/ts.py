""" Basic demo of distance geometry for transition states
"""
import sys
# import numpy
import automol

# 1. Choose reaction

#    a. bimolecular: HC(CH3)3 + OH => C(CH3)3 + H2O
RCT_ICHS = list(map(automol.smiles.inchi, ['C(C)(C)C', '[OH]']))
PRD_ICHS = list(map(automol.smiles.inchi, ['[C](C)(C)C', 'O']))
UNIMOL = False

#    b. unimolecular: (CH3)2[CH]CH2CH2O[O] => (CH3)2[C]CH2CH2O[OH]
# RCT_ICHS = list(map(automol.smiles.inchi, ['C(C)(C)CCO[O]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['[C](C)(C)CCOO']))
# UNIMOL = True

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
RXNS = automol.graph.reac.hydrogen_abstractions(RCT_GRAS, PRD_GRAS)
RXN = RXNS[0]
RCT_IDXS, PRD_IDXS = RXN.sort_order()
RCT_GEOS = list(map(RCT_GEOS.__getitem__, RCT_IDXS))
PRD_GEOS = list(map(PRD_GEOS.__getitem__, PRD_IDXS))

RXN.standardize_keys()
FRM_BND_KEYS = automol.graph.ts.forming_bond_keys(RXN.forward_ts_graph)
BRK_BND_KEYS = automol.graph.ts.breaking_bond_keys(RXN.forward_ts_graph)

# 4. Generate reactants graph and a sorted list of atom keys
RCT_GRAS = list(map(automol.convert.geom.graph, RCT_GEOS))
RCT_GRAS, _ = automol.graph.standard_keys_for_sequence(RCT_GRAS)
GRA = automol.graph.union_from_sequence(RCT_GRAS)
KEYS = sorted(automol.graph.atom_keys(GRA))
SYMS = list(map(automol.graph.atom_symbols(GRA).__getitem__, KEYS))

# 5. Generate a geometry for this reaction
if RXN.class_ == automol.par.ReactionClass.HYDROGEN_ABSTRACTION:
    # a. Join the reactant geometries together as a starting guess
    FRM_BND_KEY, = FRM_BND_KEYS
    KEY2, KEY3 = sorted(FRM_BND_KEY)
    # Since the keys are standardized, sorting puts reactant 1's atom first
    R23 = 1.6
    A123 = 180.
    A234 = 90.
    GEO_INIT = automol.geom.ts.join(*RCT_GEOS, key2=KEY2, key3=KEY3, r23=R23,
                                    a123=A123, a234=A234, angstrom=True)

    # b. Generate distance ranges for coordinates at the reaction site
    DIST_DCT = {(KEY2, KEY3): R23}
    ANG_DCT = {(None, KEY2, KEY3): A123, (KEY2, KEY3, None): A234}
    DIST_RANGE_DCT = automol.graph.embed.distance_ranges_from_coordinates(
        GRA, DIST_DCT, ang_dct=ANG_DCT, angstrom=True, degree=True, keys=KEYS)

    # c. Generate bounds matrices
    LMAT, UMAT = automol.graph.embed.join_distance_bounds_matrices(
        GRA, KEYS, DIST_RANGE_DCT, geos=RCT_GEOS, relax_torsions=UNIMOL)
    CHI_DCT = automol.graph.embed.chirality_constraint_bounds(GRA, KEYS)
    PLA_DCT = automol.graph.embed.planarity_constraint_bounds(GRA, KEYS)

if RXN.class_ == automol.par.ReactionClass.HYDROGEN_MIGRATION:
    print('HI')

sys.exit()

XMAT = automol.geom.coordinates(GEO_INIT, angstrom=True)

# 6. Optimize the coordinates to satisfy the appropriate constraints
XMAT, CONV = automol.embed.cleaned_up_coordinates(
    XMAT, LMAT, UMAT, chi_dct=CHI_DCT, pla_dct=PLA_DCT, max_dist_err=1e-1,
    dim4=False)
GEO = automol.embed.geometry_from_coordinates(XMAT, SYMS)

# 7. Print geometries
print("Reactant 1 geometry:")
print(automol.geom.string(RCT_GEOS[0]))
print()

print("Sample geometry:")
print(automol.geom.string(GEO_INIT))
print()

print("Cleaned up geometry:")
print(automol.geom.string(GEO))
print()

# 8. Check the connectivity
GRA2 = automol.geom.graph(GEO)
print("Is the graph consistent?", 'Yes' if GRA == GRA2 else 'No')
print("Is the geometry converged?", 'Yes.' if CONV else 'No.')

# 9. Generate a z-matrix
TS_BNDS = list(FRM_BND_KEYS | BRK_BND_KEYS)
ZMA = automol.convert.geom.zmatrix_x2z(GEO, TS_BNDS)
print(automol.zmat.string(ZMA))

sys.exit()
