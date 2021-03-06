""" script testing the torsion code for problem case
"""
# import autofile
import automol

ICH = automol.smiles.inchi("CCCC")
GEO = automol.inchi.geometry(ICH)
GRA = automol.geom.graph(GEO)

ZMA, ZMA_KEYS, DUMMY_KEY_DCT = automol.convert.geom.zmatrix(GEO)
ZGRA = automol.graph.relabel_for_zmatrix(GRA, ZMA_KEYS, DUMMY_KEY_DCT)

LIN_KEYS = sorted(
    automol.graph.dummy_atoms_neighbor_atom_key(ZGRA).values())
BND_KEYS = automol.graph.rotational_bond_keys(ZGRA, lin_keys=LIN_KEYS)

NAMES = {automol.zmat.torsion_coordinate_name(ZMA, *k) for k in BND_KEYS}
print('zma:\n', automol.zmat.string(ZMA, one_indexed=False))
print('torsion coordinate names:', NAMES)

# graph aligned to geometry keys
# (for getting rotational groups and symmetry numbers)
GEO, _ = automol.convert.zmat.geometry(ZMA)
GGRA = automol.graph.relabel_for_geometry(ZGRA)
print('geo:\n', automol.geom.string(GEO))

# scan_name = automol.reac.scan_coordinate(ZRXN, ZMA)
# const_names = automol.reac.constraint_coordinates(ZRXN, ZMA)
# print('scan name:', scan_name)
# print('constraint names:', const_names)
#
# #
# # GEO, GDUMMY_KEY_DCT = automol.convert.zmat.geometry(ZMA)
# # ZRXN_NEW = automol.reac.insert_dummy_atoms(GRXN, GDUMMY_KEY_DCT)
# # assert ZRXN == ZRXN_NEW
# #
# # GBND_KEYS = automol.reac.rotational_bond_keys(GRXN)
# #
# # AXES = sorted(map(sorted, GBND_KEYS))
# # for axis in AXES:
# #     print('axis:', axis)
# #     GROUPS = automol.reac.rotational_groups(GRXN, *axis)
# #     print('\tgroup 1:', GROUPS[0])
# #     print('\tgroup 2:', GROUPS[1])
# #     SYM_NUM = automol.reac.rotational_symmetry_number(GRXN, *axis)
# #     print('\tsymmetry number:', SYM_NUM)
# #
# #
# # ROTORS = automol.rotor.from_zmatrix(ZMA, zrxn=ZRXN)
# # print(ROTORS)
