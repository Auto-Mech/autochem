""" script to give a basic demo of distance geometry functionality
"""
import numpy
import automol

# 1. Choose molecule
ICH = automol.smiles.inchi('O')  # water
# ICH = automol.smiles.inchi('CO')  # methanol
# ICH = automol.smiles.inchi('C1CCCCC1')  # hexane
# ICH = automol.smiles.inchi('C1C2CC3CC1CC(C2)C3')  # adamantane

# 2. Generate graph and sorted list of atom keys
GEO = automol.inchi.geometry(ICH)
GRA = automol.geom.graph(GEO)
KEYS = sorted(automol.graph.atom_keys(GRA))
SYMS = list(map(automol.graph.atom_symbols(GRA).__getitem__, KEYS))

# 3. Generate bounds matrices
LMAT, UMAT = automol.graph.embed.distance_bounds_matrices(GRA, KEYS)
print("Lower bounds matrix:")
print(numpy.round(LMAT, 1))
print("Upper bounds matrix:")
print(numpy.round(UMAT, 1))
print()

# 4. Sample a geometry from the bounds matrices
XMAT = automol.embed.sample_raw_distance_coordinates(LMAT, UMAT, dim4=True)
GEO_INIT = automol.embed.geometry_from_coordinates(XMAT, SYMS)

# 5. Clean up the sample's coordinates
XMAT, CONV = automol.embed.cleaned_up_coordinates(XMAT, LMAT, UMAT)
GEO = automol.embed.geometry_from_coordinates(XMAT, SYMS)

# 6. Print geometries
print("Sample geometry:")
print(automol.geom.string(GEO_INIT))
print()

print("Cleaned up geometry:")
print(automol.geom.string(GEO))
print()

# 7. Check the connectivity
GRA = automol.graph.without_stereo_parities(GRA)
GRA2 = automol.geom.connectivity_graph(GEO)
print("Is the connectivity consistent?", 'Yes' if GRA == GRA2 else 'No')
