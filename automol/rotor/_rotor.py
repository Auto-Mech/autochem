"""
 Handle rotor objects

 Rotors: (Rotor1, Rotor2, ..., RotorN)
 Rotor: (tors_obj_1, tors_obj_2, ..., tors_obj_N)
"""

from itertools import chain
import yaml
import numpy
from phydat import phycon
import automol.zmat
import automol.pot
import automol.reac
from automol.rotor import _tors as tors
from automol.rotor._name import group_torsions_into_rotors
from automol.rotor._util import graph_with_keys
from automol.rotor._util import sort_tors_names


# constructors
def from_zmatrix(zma, zrxn=None, tors_names=None, multi=False):
    """ Construct a list-of-lists of torsion objects
    """

    if zrxn is None:
        # Build a graph that is used to get torsion object info
        gra, lin_keys = graph_with_keys(zma, zrxn=zrxn)

        # Build the torsion objects
        tors_lst = tors.torsion_lst(zma, gra, lin_keys)
    else:
        tors_lst = tors.reaction_torsion_lst(zma, zrxn)

    # Place the torsions into order based on rotors
    rotors = group_torsions_into_rotors(
        tors_lst, name_grps=tors_names, multi=multi)

    return rotors


def from_data(zma, tors_inf_dct, tors_names=None, multi=False):
    """ Build the rotors objects from existing data
    """

    tors_lst = ()
    for name, dct in tors_inf_dct.items():
        assert set(dct.keys()) >= {'symmetry', 'axis', 'groups'}, (
            'must have symmetry, axis, and groups in dct tp build torsions'
        )
        tors_lst += (tors.Torsion(zma, name, **dct),)

    rotors = group_torsions_into_rotors(
        tors_lst, name_grps=tors_names, multi=multi)

    return rotors


# Getters
def dimensions(rotor_lst):
    """ Get the dimensions of each of the rotors
    """
    return tuple(len(rotor) for rotor in rotor_lst)


def names(rotor_lst, flat=False):
    """ Get a flat list of names for list of rotors
    """
    _names = tuple(tuple(torsion.name for torsion in rotor)
                   for rotor in rotor_lst)
    if flat:
        _names = tuple(chain(*_names))
    return _names


def axes(rotor_lst, flat=False):
    """ Build a list of list of axes(
    """
    _axes = tuple(tuple(torsion.axis for torsion in rotor)
                  for rotor in rotor_lst)
    if flat:
        _axes = tuple(chain(*_axes))
    return _axes


def groups(rotor_lst, flat=False):
    """ Build a list of list of axes(
    """
    grps = tuple(tuple(torsion.groups for torsion in rotor)
                 for rotor in rotor_lst)
    if flat:
        grps = tuple(chain(*grps))
    return grps


def symmetries(rotor_lst, flat=False):
    """ Build a list of list of axes(
    """
    symms = tuple(tuple(torsion.symmetry for torsion in rotor)
                  for rotor in rotor_lst)
    if flat:
        symms = tuple(chain(*symms))
    return symms


def grids(rotor_lst,
          span=2.0*numpy.pi, increment=30.0*phycon.DEG2RAD, flat=False):
    """ Build a list of list of grids
    """

    rotor_lst_grids = ()
    for rotor in rotor_lst:
        rotor_grids = ()
        for torsion in rotor:
            rotor_grids += (
                automol.pot.grid(
                    torsion.zma, torsion.name,
                    span, torsion.symmetry, increment, from_equilibrium=True),
            )
        rotor_lst_grids += (rotor_grids,)
    if flat:
        rotor_lst_grids = tuple(chain(*rotor_lst_grids))

    return rotor_lst_grids


def zmatrix(rotor_lst):
    """ Get the Z-Matrix for the rotors
    """
    return rotor_lst[0][0].zma


# Manipulate the torsion objects
def relabel_for_geometry(rotor_lst):
    """ relabel the torsion objec tto correspond with a geometry converted
        from a z-matrix
    """
    geo = automol.zmat.geometry(rotor_lst[0][0].zma)
    geo_rotor_lst = tuple(
        tuple(tors.relabel_for_geometry(torsion) for torsion in rotor)
        for rotor in rotor_lst)

    return geo, geo_rotor_lst


# I/O
def string(rotor_lst):
    """ Write a list torsions to a string
    """

    def _encode_idxs(idxs):
        if len(idxs) == 1:
            idx_str = idxs[0]+1
        else:
            idx_str = '-'.join(str(val+1) for val in idxs)
        return idx_str

    tors_dct = {}
    for rotor in rotor_lst:
        for torsion in rotor:
            _axis = torsion.axis
            _grps = torsion.groups
            tors_dct[torsion.name] = {
                'axis1': _axis[0]+1,
                'group1': _encode_idxs(_grps[0]),
                'axis2': _axis[1]+1,
                'group2': _encode_idxs(_grps[1]),
                'symmetry': torsion.symmetry,
            }

    sort_tors_dct = {}
    tors_names = sort_tors_names(list(tors_dct.keys()))
    for name in tors_names:
        sort_tors_dct[name] = tors_dct[name]

    tors_str = yaml.dump(sort_tors_dct, sort_keys=False)

    return tors_str


def from_string(tors_str):
    """ read the transformation from a string
    """

    def _decode_idxs(idxs_str):
        if isinstance(idxs_str, int):
            idxs = (int(idxs_str)-1,)
        else:
            idxs = tuple(map(int, idxs_str.split('-')))
            idxs = tuple(val-1 for val in idxs)
        return idxs

    inf_dct = {}

    tors_dct = yaml.load(tors_str, Loader=yaml.FullLoader)
    for name, dct in tors_dct.items():
        _axis = (dct['axis1']-1, dct['axis2']-1)
        _grps = (_decode_idxs(dct['group1']), _decode_idxs(dct['group2']))
        symm = dct['symmetry']

        inf_dct[name] = {'axis': _axis, 'groups': _grps, 'symmetry': symm}
        # torsions += (Torsion('', name, _axis, _grps, symm),)

    return inf_dct


if __name__ == '__main__':
    # RCT_ICHS = list(map(automol.smiles.inchi, ['[O]O', 'CCC=C[CH]CCCCC']))
    # PRD_ICHS = list(map(automol.smiles.inchi, ['O=O', 'CCCC=CCCCCC']))

    # RCT_GEOS = list(map(automol.inchi.geometry, RCT_ICHS))
    # PRD_GEOS = list(map(automol.inchi.geometry, PRD_ICHS))

    # RCT_GRAS = list(map(automol.geom.graph, RCT_GEOS))
    # PRD_GRAS = list(map(automol.geom.graph, PRD_GEOS))

    # RCT_GRAS = list(map(automol.graph.without_stereo_parities, RCT_GRAS))
    # PRD_GRAS = list(map(automol.graph.without_stereo_parities, PRD_GRAS))

    # RCT_GRAS, _ = automol.graph.standard_keys_for_sequence(RCT_GRAS)
    # PRD_GRAS, _ = automol.graph.standard_keys_for_sequence(PRD_GRAS)

    # RXNS = automol.reac.find(RCT_GRAS, PRD_GRAS)
    # RXN = RXNS[0]

    # RXN, RCT_GEOS, PRD_GEOS = (
    #     automol.reac.standard_keys_with_sorted_geometries(
    #         RXN, RCT_GEOS, PRD_GEOS))

    # GEO = automol.reac.ts_geometry(RXN, RCT_GEOS, log=False)

    # print(automol.geom.string(GEO))

    # with open('zmat.xyz', 'r') as f:
    #     GEO_STR = f.read()
    # with open('zmat.r.yaml', 'r') as f:
    #     RXN_STR = f.read()
    # GEO = autofile.data_types.sread.geometry(GEO_STR)
    # RXN = autofile.data_types.sread.reaction(RXN_STR)

    # ZMA, ZMA_KEYS, DUMMY_KEY_DCT = automol.reac.ts_zmatrix(RXN, GEO)

    # ZRXN = automol.reac.relabel_for_zmatrix(RXN, ZMA_KEYS, DUMMY_KEY_DCT)
    # ZTSG = ZRXN.forward_ts_graph

    # print('zma:\n', automol.zmat.string(ZMA))

    # LIN_KEYS = sorted(
    #     automol.graph.dummy_atoms_neighbor_atom_key(ZTSG).values())
    # BND_KEYS = automol.graph.rotational_bond_keys(ZTSG, lin_keys=LIN_KEYS)
    # print(LIN_KEYS)
    # print(BND_KEYS)

    # NAMES = [automol.zmat.torsion_coordinate_name(ZMA, *k) for k in BND_KEYS]
    # print(NAMES)

    # ROTORS = from_zmatrix(ZMA, zrxn=ZRXN)
    # print(ROTORS)

    # ZMA_STR = autofile.io_.read_file('zmat.zmat')
    # ZMA = autofile.data_types.sread.zmatrix(ZMA_STR)
    # GEO = automol.zmat.geometry(ZMA)
    # R1GEO = GEO[:3]
    # R2GEO = GEO[3:]
    # R1GRA = automol.geom.graph(R1GEO)
    # R2GRA = automol.geom.graph(R2GEO)

    # RXN_STR = autofile.io_.read_file('zmat.r.yaml')

    # RXN = autofile.data_types.sread.reaction(RXN_STR)
    # RCT_GRAS = automol.reac.reactant_graphs(RXN)
    # ISO1 = automol.graph.full_isomorphism(RCT_GRAS[0], R1GRA)
    # ISO2 = automol.graph.full_isomorphism(RCT_GRAS[1], R2GRA)
    # print(automol.graph.inchi(RCT_GRAS[1]))
    # print(automol.graph.inchi(R2GRA))
    # ISO = {**ISO1, **ISO2}
    # print(ISO)

    # # from_zmatrix(ZMA, zrxn=RXN)
