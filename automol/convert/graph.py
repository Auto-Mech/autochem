""" graph conversions
"""
import autoparse.pattern as app
import autoparse.find as apf
from automol.util import dict_
from automol.convert.inchi import standard_form
from automol.convert.inchi import object_to_hardcoded_inchi_by_key
from automol.convert import _molfile
from automol.convert import _rdkit
from automol.convert import _util
from automol.graph._graph_dep import atom_keys
from automol.graph._graph_dep import bond_keys
from automol.graph._graph_dep import atom_symbols
from automol.graph._graph_dep import bond_orders
# getters
from automol.graph._graph_dep import without_dummy_atoms
from automol.graph._graph_dep import atom_bond_valences
from automol.graph._graph_dep import explicit
from automol.graph._graph_dep import atom_unsaturated_valences
# stereo
from automol.graph._graph_dep import has_stereo
from automol.graph._graph_dep import dominant_resonance
# dep
from automol.graph._embed_dep import fake_stereo_geometry
from automol.graph._embed_dep import geometry as embed_geometry
from automol.graph._embed_dep import backbone_isomorphic
from automol.graph.geom import coordinates


# graph => inchi
def inchi(gra, stereo=True):
    """ Generate an InChI string from a molecular graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: str
    """

    ich = object_to_hardcoded_inchi_by_key(
        'graph', gra, comp=_compare)

    if ich is None:
        if not stereo or not has_stereo(gra):
            ich, _ = inchi_with_sort_from_geometry(gra)
            ich = standard_form(ich, stereo=stereo)
        else:
            gra = explicit(gra)
            geo, geo_idx_dct = fake_stereo_geometry(gra)
            ich, _ = inchi_with_sort_from_geometry(
                gra, geo=geo, geo_idx_dct=geo_idx_dct)

    return ich


def _compare(gra1, gra2):
    """ Compare the backbone structure of two moleculare graphs.

        :param gra1: molecular graph 1
        :type gra1: automol graph data structure
        :param gra2: molecular graph 2
        :type gra2: automol graph data structure
        :rtype: bool
    """

    gra1 = without_dummy_atoms(gra1)
    gra2 = without_dummy_atoms(gra2)

    return backbone_isomorphic(gra1, gra2)


def inchi_with_sort_from_geometry(gra, geo=None, geo_idx_dct=None):
    """ Generate an InChI string from a molecular graph.
        If coordinates are passed in, they are used to determine stereo.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct:
        :type geo_idx_dct: dict[:]
        :rtype: (str, tuple(int))
    """
    gra = without_dummy_atoms(gra)
    gra = dominant_resonance(gra)
    atm_keys = sorted(atom_keys(gra))
    bnd_keys = list(bond_keys(gra))
    atm_syms = dict_.values_by_key(atom_symbols(gra), atm_keys)
    atm_bnd_vlcs = dict_.values_by_key(
        atom_bond_valences(gra), atm_keys)
    atm_rad_vlcs = dict_.values_by_key(
        atom_unsaturated_valences(gra), atm_keys)
    bnd_ords = dict_.values_by_key(bond_orders(gra), bnd_keys)

    if geo is not None:
        assert geo_idx_dct is not None
        atm_xyzs = coordinates(geo)
        atm_xyzs = [atm_xyzs[geo_idx_dct[atm_key]] if atm_key in geo_idx_dct
                    else (0., 0., 0.) for atm_key in atm_keys]
    else:
        atm_xyzs = None

    mlf, key_map_inv = _molfile.from_data(
        atm_keys, bnd_keys, atm_syms, atm_bnd_vlcs, atm_rad_vlcs, bnd_ords,
        atm_xyzs=atm_xyzs)
    rdm = _rdkit.from_molfile(mlf)
    ich, aux_info = _rdkit.to_inchi(rdm, with_aux_info=True)
    nums = _parse_sort_order_from_aux_info(aux_info)
    nums = tuple(map(key_map_inv.__getitem__, nums))

    return ich, nums


def _parse_sort_order_from_aux_info(aux_info):
    ptt = app.escape('/N:') + app.capturing(
        app.series(app.UNSIGNED_INTEGER, ','))
    num_str = apf.first_capture(ptt, aux_info)
    nums = tuple(map(int, num_str.split(',')))
    return nums


def geometry(gra):
    """ Convert a molecular graph to a molecular geometry.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :rtype: automol molecular geometry data structure
    """

    symbs = atom_symbols(gra)
    if len(symbs) != 1:
        gra = explicit(gra)
        geo = embed_geometry(gra)
    else:
        symb = list(symbs.values())[0]
        # symb = list(symbs.keys())[0]
        geo = ((symb, (0.00, 0.00, 0.00)),)

    return geo


def formula(gra):
    """ Generate a stoichiometric formula dictionary from a molecular graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :type: dict[str: int]
    """

    gra = explicit(gra)
    syms = atom_symbols(gra).values()
    fml = _util.formula(syms)

    return fml
#
#
# if __name__ == '__main__':
#     import automol
#
#     for ICH in [
#             'InChI=1S/C4H7O2/c1-3-4(2)6-5/h3-5H,1-2H2/t4-/m0/s1',
#             'InChI=1S/C4H7O/c1-4(2)3-5/h3-4H,1H2,2H3/t4-/m1/s1',
#             'InChI=1S/C5H7/c1-3-5-4-2/h1,5H,4H2,2H3',
#             'InChI=1S/C5H9/c1-4-5(2)3/h4-5H,1-2H2,3H3/t5-/m1/s1',
#             'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b4-3+',
#             'InChI=1S/C5H7O/c1-5-3-2-4-6-5/h2-5H,1H3/t5-/m0/s1',
#             'InChI=1S/C6H11/c1-5(2)6(3)4/h5H,1,3H2,2,4H3/t5-/m1/s1',
#             'InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-2,4-5H,3H2',
#             'InChI=1S/C8H15O2/c1-7(2)5-8(3,4)6-10-9/h5,9H,3,6H2,1-2,'
#             '4H3/t8-/m0/s1', ]:
#         GEO = automol.inchi.geometry(ICH)
#         print(automol.geom.string(GEO))
#         print()
#         GRA = automol.geom.graph(GEO)
#         RAD_GRP_DCT = automol.graph.radical_group_dct(GRA)
#         for ATM, GRPS in RAD_GRP_DCT.items():
#             print(len(GRPS))
#             for GRP in GRPS:
#                 print(len(GRP))
#                 print('atom', ATM)
#                 print('group', GRP)
#                 ATM_KEYS = automol.graph.atom_keys(GRP)
#                 GRP_GEO = automol.geom.from_subset(GEO, ATM_KEYS)
#                 print(automol.geom.string(GRP_GEO))
#                 print()
#
#                 print(automol.graph.string(GRP, one_indexed=False))
#                 STE_ATM_KEYS = automol.graph.stereogenic_atom_keys(GRP)
#                 print(STE_ATM_KEYS)
#                 # GRP_ICH = automol.graph.inchi(GRP, stereo=False)
#                 # GRP_ICH = automol.inchi.add_stereo(GRP_ICH)
#                 # print(GRP_ICH)
#                 GRP_ICH = automol.graph.inchi(GRP, stereo=True)
