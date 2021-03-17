""" test ring zmatrix
"""
import automol

SMI = 'C1CCC1C'
ICH = automol.smiles.inchi(SMI)
GEO = automol.inchi.geometry(ICH)
ZMA = automol.geom.zmatrix(GEO)
# print(automol.geom.string(GEO))
print(automol.zmat.string(ZMA))
GEO = automol.zmat.geometry(ZMA)
print(automol.geom.string(GEO))
