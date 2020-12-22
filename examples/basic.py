""" script to give a basic demo of distance geometry functionality
"""
import automol

# ICH = automol.smiles.inchi('O')  # water
# ICH = automol.smiles.inchi('CO')  # methanol
# ICH = automol.smiles.inchi('C1CCCCC1')  # hexane
ICH = automol.smiles.inchi('C1C2CC3CC1CC(C2)C3')  # adamantane
GEO = automol.inchi.geometry(ICH)
GRA = automol.geom.graph(GEO)

KEYS = sorted(automol.graph.atom_keys(GRA))
LMAT, UMAT = automol.graph.embed.distance_bounds_matrices(GRA, KEYS)
XMAT = automol.embed.sample_raw_distance_coordinates(LMAT, UMAT, dim4=True)

XMAT, CONV = automol.embed.cleaned_up_coordinates(XMAT, LMAT, UMAT)
print(CONV)

SYMS = list(map(automol.graph.atom_symbols(GRA).__getitem__, KEYS))
GEO = automol.embed.geometry_from_coordinates(XMAT, SYMS)
print(automol.geom.string(GEO))

GRA = automol.graph.without_stereo_parities(GRA)
GRA2 = automol.geom.connectivity_graph(GEO)
print(GRA == GRA2)
