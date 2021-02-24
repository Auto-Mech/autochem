""" script testing the torsion code for problem case
"""
import autofile
import automol

ZMA_STR = autofile.io_.read_file('1/zmat.zmat')
ZRXN_STR = autofile.io_.read_file('1/zmat.r.yaml')

ZMA = autofile.data_types.sread.zmatrix(ZMA_STR)
ZRXN = autofile.data_types.sread.reaction(ZRXN_STR)
print('zrxn:\n', ZRXN)
print('\nzma:\n', automol.zmat.string(ZMA))

RCT_GRAS = automol.reac.reactant_graphs(ZRXN)
PRD_GRAS = automol.reac.product_graphs(ZRXN)

RCT_ICHS = list(map(automol.graph.inchi, RCT_GRAS))
PRD_ICHS = list(map(automol.graph.inchi, PRD_GRAS))

RCT_SMIS = list(map(automol.inchi.smiles, RCT_ICHS))
PRD_SMIS = list(map(automol.inchi.smiles, PRD_ICHS))
print(RCT_SMIS)
print(PRD_SMIS)

# # You can also do this to determine linear atoms from zmatrix:
# bnd_keys = automol.reac.rotational_bond_keys(ZRXN, zma=ZMA)
ZBND_KEYS = automol.reac.rotational_bond_keys(ZRXN)
# print('rot bond_keys:', ZBND_KEYS)
print('rot bond_keys:')
for key in ZBND_KEYS:
    print(key)
    automol.zmat.torsion_coordinate_name(ZMA, *key)

NAMES = {automol.zmat.torsion_coordinate_name(ZMA, *k) for k in ZBND_KEYS}
print('zma:\n', automol.zmat.string(ZMA, one_indexed=False))
print('torsion coordinate names:', NAMES)

# graph aligned to geometry keys
# (for getting rotational groups and symmetry numbers)
GEO, _ = automol.convert.zmat.geometry(ZMA)
GRXN = automol.reac.relabel_for_geometry(ZRXN)
print('rxn obj')
print(automol.graph.string(ZRXN.forward_ts_graph))
print()
print(automol.graph.string(GRXN.forward_ts_graph))

print('geo:\n', automol.geom.string(GEO))

print('ZRXN')
print(ZRXN)
print('GRXN')
print(GRXN)

GEO, GDUMMY_KEY_DCT = automol.convert.zmat.geometry(ZMA)
ZRXN_NEW = automol.reac.insert_dummy_atoms(GRXN, GDUMMY_KEY_DCT)

GBND_KEYS = automol.reac.rotational_bond_keys(GRXN)

AXES = sorted(map(sorted, GBND_KEYS))
for axis in AXES:
    print('axis:', axis)
    GROUPS = automol.reac.rotational_groups(GRXN, *axis)
    print('\tgroup 1:', GROUPS[0])
    print('\tgroup 2:', GROUPS[1])
    SYM_NUM = automol.reac.rotational_symmetry_number(GRXN, *axis)
    print('\tsymmetry number:', SYM_NUM)


ROTORS = automol.rotor.from_zmatrix(ZMA, zrxn=ZRXN)
print(ROTORS)
print('zma:\n', automol.zmat.string(ZMA))
rnames = automol.rotor.names(ROTORS, flat=True)
raxes = automol.rotor.axes(ROTORS, flat=True)
rgrps = automol.rotor.groups(ROTORS, flat=True)
rsymms = automol.rotor.symmetries(ROTORS, flat=True)
for name, axis, grps, symm in zip(rnames, raxes, rgrps, rsymms):
    print(name)
    print(axis)
    print(grps)
    print(symm)
