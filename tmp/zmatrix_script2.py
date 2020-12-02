""" script to develop new z-matrix code
"""
import automol
from automol.graph._geom import longest_chain

ICH = automol.smiles.inchi('CCCCCCC')
gra = automol.inchi.graph(ICH)

gra = automol.graph.implicit(gra)

chain = longest_chain(gra)
print(chain)
