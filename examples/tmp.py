""" testing fake_stereo_geometry function
"""
import automol

ICH = automol.smiles.inchi('FC=CC=CCCC(O)(C)')
GEO_IN = automol.inchi.geometry(ICH)
GRA = automol.geom.graph(GEO_IN)
automol.graph.embed.fake_stereo_geometry(GRA)
