""" script to develop new z-matrix code
"""
import automol
from automol.graph._geom import _start_zmatrix_from_atom
from automol.graph._geom import longest_chain
from automol.graph._geom import atom_symbols
# from qcelemental import constants as qcc

# CC_DIST = 1.5 * qcc.conversion_factor('angstrom', 'bohr')

# 0. Copy-pasted code to get things started
ICH = automol.smiles.inchi('CCCCCCC')
gra = automol.graph.explicit(automol.inchi.graph(ICH))

atm_sym_dct = atom_symbols(gra)

chain = longest_chain(gra)

if atm_sym_dct[chain[0]] != 'H':
    chain = list(reversed(chain))

if len(chain) > 1:
    atm_key = chain[1]
else:
    atm_key = chain[0]

zma, zma_key_dct, dummy_atm_key, gra = _start_zmatrix_from_atom(gra, atm_key)

print(automol.zmatrix.string(zma, one_indexed=False))
print(zma_key_dct)

print(dummy_atm_key)
atm1_key, atm2_key, atm3_key = chain[:3]

# start building z-matrix

# end building z-matrix

geo = automol.zmatrix.geometry(zma)
print(automol.geom.xyz_string(geo))
