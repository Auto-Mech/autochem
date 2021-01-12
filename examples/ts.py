""" Basic demo of distance geometry for transition states
"""
import sys
# import numpy
import automol

# 1. Choose reaction

#    a. hydrogen abstraction: HC(CH3)3 + OH => C(CH3)3 + H2O
# RCT_ICHS = list(map(automol.smiles.inchi, ['C(C)(C)C', '[OH]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['[C](C)(C)C', 'O']))

#    b. hydrogen migration: (CH3)2[CH]CH2CH2O[O] => (CH3)2[C]CH2CH2O[OH]
# RCT_ICHS = list(map(automol.smiles.inchi, ['C(C)(C)CCO[O]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['[C](C)(C)CCOO']))

#    c. addition: CH3CH2[CH2] + OO => CH3CH2CH2OO
# RCT_ICHS = list(map(automol.smiles.inchi, ['CC[CH2]', '[O][O]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['CCCO[O]']))

#    d. beta-scission: CH3CH2CH2OO => CH3CH2[CH2] + OO
# RCT_ICHS = list(map(automol.smiles.inchi, ['CCCO[O]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['CC[CH2]', '[O][O]']))

#    e. ring-forming scission: [CH2]CH2OOH =>
# RCT_ICHS = list(map(automol.smiles.inchi, ['[CH2]COO']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['C1CO1', '[OH]']))

#    f. elimination: CH3CH2OO => C2H4 + HOO
# RCT_ICHS = list(map(automol.smiles.inchi, ['CCO[O]']))
# PRD_ICHS = reversed(list(map(automol.smiles.inchi, ['C=C', 'O[O]'])))

#    g. insertion: CH3CH3 + CH2 => C2H4 + HOO
RCT_ICHS = list(map(automol.smiles.inchi, ['CC', '[CH2]']))
PRD_ICHS = list(map(automol.smiles.inchi, ['CCC']))

#    h. substitution: H2O2 + H => H2O + OH
# RCT_ICHS = list(map(automol.smiles.inchi, ['OO', '[H]']))
# PRD_ICHS = list(map(automol.smiles.inchi, ['O', '[OH]']))

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
RXNS = automol.graph.reac.classify(RCT_GRAS, PRD_GRAS)
RXN = RXNS[0]
# Sort and standardize keys (must go together)
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
    # a. Read in the formed bonds and set the reaction site coordinates
    FRM_BND_KEY, = FRM_BND_KEYS
    KEY2, KEY3 = sorted(FRM_BND_KEY)
    R23 = 1.6
    A123 = 170.
    A234 = 85.
    RELAX_TORS = False

    # b. Use the reactant(s) for the initial geometry
    GEO_INIT = automol.geom.ts.join(*RCT_GEOS, key2=KEY2, key3=KEY3, r23=R23,
                                    a123=A123, a234=A234)

    # c. Generate distance ranges for coordinates at the reaction site
    DIST_DCT = {(KEY2, KEY3): R23}
    ANG_DCT = {(None, KEY2, KEY3): A123, (KEY2, KEY3, None): A234}
    DIST_RANGE_DCT = automol.graph.embed.distance_ranges_from_coordinates(
        GRA, DIST_DCT, ang_dct=ANG_DCT, degree=True, keys=KEYS)

    # d. Generate bounds matrices
    LMAT, UMAT = automol.graph.embed.join_distance_bounds_matrices(
        GRA, KEYS, DIST_RANGE_DCT, geos=RCT_GEOS, relax_torsions=RELAX_TORS)
    CHI_DCT = automol.graph.embed.chirality_constraint_bounds(GRA, KEYS)
    PLA_DCT = automol.graph.embed.planarity_constraint_bounds(GRA, KEYS)

if RXN.class_ == automol.par.ReactionClass.HYDROGEN_MIGRATION:
    # a. Read in the formed bonds and set the reaction site coordinates
    FRM_BND_KEY, = FRM_BND_KEYS
    KEY2, KEY3 = sorted(FRM_BND_KEY)
    R23 = 1.6
    RELAX_TORS = True

    # b. Use the reactant(s) for the initial geometry
    GEO_INIT, = RCT_GEOS

    # b. Generate distance ranges for coordinates at the reaction site
    DIST_DCT = {(KEY2, KEY3): R23}
    DIST_RANGE_DCT = automol.graph.embed.distance_ranges_from_coordinates(
        GRA, DIST_DCT, keys=KEYS)

    # d. Generate bounds matrices
    LMAT, UMAT = automol.graph.embed.join_distance_bounds_matrices(
        GRA, KEYS, DIST_RANGE_DCT, geos=RCT_GEOS, relax_torsions=RELAX_TORS)
    CHI_DCT = automol.graph.embed.chirality_constraint_bounds(GRA, KEYS)
    PLA_DCT = automol.graph.embed.planarity_constraint_bounds(GRA, KEYS)

if RXN.class_ == automol.par.ReactionClass.ADDITION:
    # a. Read in the formed bonds and set the reaction site coordinates
    FRM_BND_KEY, = FRM_BND_KEYS
    KEY2, KEY3 = sorted(FRM_BND_KEY)
    R23 = 1.9
    A123 = 85.
    A234 = 85.
    D1234 = 85.
    RELAX_TORS = False

    # c. Generate distance ranges for coordinates at the reaction site
    # NOTE: NEED TO PUT IN GEOMETRY PARAMETERS AS WELL, TO AVOID OVERWRITING
    # THEM WITH HEURISTICS
    DIST_DCT = {(KEY2, KEY3): R23}
    ANG_DCT = {(None, KEY2, KEY3): A123, (KEY2, KEY3, None): A234}
    DIH_DCT = {(None, KEY2, KEY3, None): D1234}
    DIST_RANGE_DCT = automol.graph.embed.distance_ranges_from_coordinates(
        GRA, DIST_DCT, ang_dct=ANG_DCT, dih_dct=DIH_DCT, degree=True,
        keys=KEYS)

    # b. Use the reactant(s) for the initial geometry
    GEO_INIT = automol.geom.ts.join(*RCT_GEOS, key2=KEY2, key3=KEY3, r23=R23)

    # d. Generate bounds matrices
    LMAT, UMAT = automol.graph.embed.join_distance_bounds_matrices(
        GRA, KEYS, DIST_RANGE_DCT, geos=RCT_GEOS, relax_torsions=RELAX_TORS)
    CHI_DCT = automol.graph.embed.chirality_constraint_bounds(GRA, KEYS)
    PLA_DCT = automol.graph.embed.planarity_constraint_bounds(GRA, KEYS)

if RXN.class_ == automol.par.ReactionClass.BETA_SCISSION:
    # do absolutely nothing
    GEO_INIT, = RCT_GEOS
    RELAX_TORS = False

    DIST_RANGE_DCT = {}

    # d. Generate bounds matrices
    LMAT, UMAT = automol.graph.embed.join_distance_bounds_matrices(
        GRA, KEYS, DIST_RANGE_DCT, geos=RCT_GEOS, relax_torsions=RELAX_TORS)
    CHI_DCT = automol.graph.embed.chirality_constraint_bounds(GRA, KEYS)
    PLA_DCT = automol.graph.embed.planarity_constraint_bounds(GRA, KEYS)

if RXN.class_ == automol.par.ReactionClass.RING_FORM_SCISSION:
    # a. Read in the formed bonds and set the reaction site coordinates
    FRM_BND_KEY, = FRM_BND_KEYS
    BRK_BND_KEY, = BRK_BND_KEYS
    KEY2, = FRM_BND_KEY & BRK_BND_KEY
    print('KEY2', KEY2)
    print('FRM_BND_KEY', FRM_BND_KEY)
    print('BRK_BND_KEY', BRK_BND_KEY)
    KEY1, = FRM_BND_KEY - {KEY2}
    KEY3, = BRK_BND_KEY - {KEY2}

    R12 = 1.9  # set the forming bond to 2.2 angstroms
    R23 = 1.5  # set the breaking bond to 1.9 angstroms
    A234 = 85.
    D1234 = 170.

    RELAX_TORS = True

    # c. Generate distance ranges for coordinates at the reaction site
    # NOTE: NEED TO PUT IN GEOMETRY PARAMETERS AS WELL, TO AVOID OVERWRITING
    # THEM WITH HEURISTICS
    DIST_DCT = {(KEY1, KEY2): R12, (KEY2, KEY3): R23}
    ANG_DCT = {(KEY2, KEY3, None): A234}
    DIH_DCT = {(KEY1, KEY2, KEY3, None): D1234}
    DIST_RANGE_DCT = automol.graph.embed.distance_ranges_from_coordinates(
        GRA, DIST_DCT, ang_dct=ANG_DCT, dih_dct=DIH_DCT, degree=True,
        keys=KEYS)

    # b. Use the reactant(s) for the initial geometry
    GEO_INIT, = RCT_GEOS

    # d. Generate bounds matrices
    LMAT, UMAT = automol.graph.embed.join_distance_bounds_matrices(
        GRA, KEYS, DIST_RANGE_DCT, geos=RCT_GEOS, relax_torsions=RELAX_TORS)
    CHI_DCT = automol.graph.embed.chirality_constraint_bounds(GRA, KEYS)
    PLA_DCT = automol.graph.embed.planarity_constraint_bounds(GRA, KEYS)

if RXN.class_ == automol.par.ReactionClass.ELIMINATION:
    # a. Read in the formed bonds and set the reaction site coordinates
    FRM_BND_KEY, = FRM_BND_KEYS
    KEY2, KEY3 = sorted(FRM_BND_KEY)
    R23 = 1.6
    RELAX_ANG = True
    RELAX_TORS = True

    # Also handle the angles for the forming ring
    FRM_RNG_KEYS, = automol.graph.ts.forming_rings_atom_keys(
        RXN.forward_ts_graph)

    # b. Use the reactant(s) for the initial geometry
    GEO_INIT, = RCT_GEOS

    # b. Generate distance ranges for coordinates at the reaction site
    DIST_DCT = {(KEY2, KEY3): R23}
    DIST_RANGE_DCT = automol.graph.embed.distance_ranges_from_coordinates(
        GRA, DIST_DCT, keys=KEYS, rings_keys=[FRM_RNG_KEYS])

    # d. Generate bounds matrices
    LMAT, UMAT = automol.graph.embed.join_distance_bounds_matrices(
        GRA, KEYS, DIST_RANGE_DCT, geos=RCT_GEOS, relax_torsions=RELAX_TORS,
        relax_angles=RELAX_ANG)
    CHI_DCT = automol.graph.embed.chirality_constraint_bounds(GRA, KEYS)
    PLA_DCT = automol.graph.embed.planarity_constraint_bounds(GRA, KEYS)

if RXN.class_ == automol.par.ReactionClass.INSERTION:
    print("HI INSERTION")

if RXN.class_ == automol.par.ReactionClass.SUBSTITUTION:
    print("HI SUBSTITUTION")

sys.exit()

XMAT = automol.geom.coordinates(GEO_INIT, angstrom=True)

# 6. Optimize the coordinates to satisfy the appropriate constraints
XMAT, CONV = automol.embed.cleaned_up_coordinates(
    XMAT, LMAT, UMAT, chi_dct=CHI_DCT, pla_dct=PLA_DCT, max_dist_err=1e-1,
    dim4=False, log=True)

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
