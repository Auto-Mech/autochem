""" Basic demo of distance geometry functionality
"""
import numpy
import automol

# 1. Choose molecule
# ICH = automol.smiles.inchi('C=C')
# ICH = automol.smiles.inchi('C=N')
ICH = automol.smiles.inchi(r'[H]\N=C\\C(=N)\C=N\\[H]')

# 2. Generate graph and sorted list of atom keys
GEO = automol.inchi.geometry(ICH)
GRA = automol.geom.graph(GEO)
KEYS = sorted(automol.graph.atom_keys(GRA))
SYMS = list(map(automol.graph.atom_symbols(GRA).__getitem__, KEYS))

# 3. Generate bounds matrices
LMAT, UMAT = automol.graph.embed.distance_bounds_matrices(GRA, KEYS)
CHIP_DCT = automol.graph.embed.planarity_constraint_bounds(GRA, KEYS)
print("Lower bounds matrix:")
print(numpy.round(LMAT, 1))
print("Upper bounds matrix:")
print(numpy.round(UMAT, 1))
print()
print(CHIP_DCT)

# 4. Sample a geometry from the bounds matrices
XMAT = automol.embed.sample_raw_distance_coordinates(LMAT, UMAT, dim4=True)
GEO_INIT = automol.embed.geometry_from_coordinates(XMAT, SYMS)

# 5. Clean up the sample's coordinates
XMAT, CONV = automol.embed.cleaned_up_coordinates(
    XMAT, LMAT, UMAT, chip_dct=CHIP_DCT)
GEO = automol.embed.geometry_from_coordinates(XMAT, SYMS)

# 6. Print the largest errors
DMAT = automol.embed.distance_matrix_from_coordinates(XMAT)
ERR_DCT = automol.embed.greatest_distance_errors(DMAT, LMAT, UMAT)
SP_DCT = automol.graph.atom_shortest_paths(GRA)
for (KEY1, KEY2), ERR in ERR_DCT.items():
    print('\tError:', ERR)
    if KEY2 not in SP_DCT[KEY1]:
        print('\nNot connected:', KEY1, KEY2)
    else:
        print('\tPath:', SP_DCT[KEY1][KEY2])

# 7. Print geometries
print("Sample geometry:")
print(automol.geom.string(GEO_INIT))
print()

print("Cleaned up geometry:")
print(automol.geom.string(GEO))
print()

# 8. Check the connectivity
GRA2 = automol.geom.graph(GEO)
print("Is the graph the same?", 'Yes' if GRA == GRA2 else 'No')
