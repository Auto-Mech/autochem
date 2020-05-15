""" graph conversions
"""
from automol import dict_
import automol.graph
import automol.inchi
import automol.convert.inchi
from automol.convert import _molfile as _molf
from automol.convert import _rdkit
from automol.convert import _util


def geometry(gra):
    """ graph => geometry
    """
    gra = automol.graph.explicit(gra)
    geo, _ = automol.graph.heuristic_geometry(gra)
    return geo


def formula(gra):
    """ graph  => formula
    """
    gra = automol.graph.explicit(gra)
    syms = automol.graph.atom_symbols(gra).values()
    fml = _util.formula(syms)
    return fml


# graph => inchi
def inchi(gra):
    """ graph => inchi
    """
    gras = automol.graph.connected_components(gra)
    ichs = list(map(_connected_inchi, gras))
    ich = automol.inchi.join(ichs)
    ich = automol.inchi.standard_form(ich)
    return ich


def inchi_from_geometry(gra, geo=None, geo_idx_dct=None):
    """ graph => inchi with geometry for stereo information
    """


# hardcoded inchis which neither RDKit nor Pybel can handle
# (any molecule with more than three unpaired electrons on an atom will fail)
HARDCODED_INCHI_GRAPHS_DCT = {
    'InChI=1S/C': ({0: ('C', 0, None)}, {}),
    'InChI=1S/B': ({0: ('B', 0, None)}, {}),
    'InChI=1S/N': ({0: ('N', 0, None)}, {}),
    'InChI=1S/CH/h1H': ({0: ('C', 1, None)}, {}),
    'InChI=1S/CF/c1-2': ({0: ('C', 0, None), 1: ('F', 0, None)},
                         {frozenset({0, 1}): (1, None)}),
    'InChI=1S/CCl/c1-2': ({0: ('C', 0, None), 1: ('Cl', 0, None)},
                          {frozenset({0, 1}): (1, None)}),
    'InChI=1S/CBr/c1-2': ({0: ('C', 0, None), 1: ('Br', 0, None)},
                          {frozenset({0, 1}): (1, None)}),
    'InChI=1S/CI/c1-2': ({0: ('C', 0, None), 1: ('I', 0, None)},
                         {frozenset({0, 1}): (1, None)}),
}


def _connected_inchi(gra):
    if not automol.graph.has_stereo(gra):
        ich = connected_inchi_from_geometry(gra)
        ich = automol.inchi.standard_form(ich, remove_stereo=True)
    else:
        gra = automol.graph.explicit(gra)
        geo, geo_idx_dct = automol.graph.heuristic_geometry(gra)
        ich = connected_inchi_from_geometry(
            gra, geo=geo, geo_idx_dct=geo_idx_dct)

    return ich


def connected_inchi_from_geometry(gra, geo=None, geo_idx_dct=None):
    """ a more general function that takes the stereo geometry explicitly
    """
    ich = None

    # first, see if this is in the list of bad-inchi graphs
    for ich_, gra_ in HARDCODED_INCHI_GRAPHS_DCT.items():
        if _graphs_are_equivalent(gra, gra_):
            ich = ich_

    # otherwise, do it the regular way -- using a molfile
    if ich is None:
        mlf, _ = _molfile(gra, geo=geo, geo_idx_dct=geo_idx_dct)
        rdm = _rdkit.from_molfile(mlf)
        ich = _rdkit.to_inchi(rdm)

    return ich


def _graphs_are_equivalent(gra1, gra2):
    gra1 = automol.graph.without_dummy_atoms(gra1)
    gra2 = automol.graph.without_dummy_atoms(gra2)
    return automol.graph.backbone_isomorphic(gra1, gra2)


def _molfile(gra, geo=None, geo_idx_dct=None):
    """
    """
    gra = automol.graph.without_dummy_atoms(gra)
    gra = automol.graph.dominant_resonance(gra)
    atm_keys = sorted(automol.graph.atom_keys(gra))
    bnd_keys = list(automol.graph.bond_keys(gra))
    atm_syms = dict_.values_by_key(automol.graph.atom_symbols(gra), atm_keys)
    atm_bnd_vlcs = dict_.values_by_key(
        automol.graph.atom_bond_valences(gra), atm_keys)
    atm_rad_vlcs = dict_.values_by_key(
        automol.graph.atom_unsaturated_valences(gra), atm_keys)
    bnd_ords = dict_.values_by_key(automol.graph.bond_orders(gra), bnd_keys)

    if geo is not None:
        assert geo_idx_dct is not None
        atm_xyzs = automol.geom.coordinates(geo)
        atm_xyzs = [atm_xyzs[geo_idx_dct[atm_key]] for atm_key in atm_keys]
    else:
        atm_xyzs = None

    mlf, idx_dct = _molf.from_data(
        atm_keys, bnd_keys, atm_syms, atm_bnd_vlcs, atm_rad_vlcs, bnd_ords,
        atm_xyzs=atm_xyzs)

    return mlf, idx_dct


# *** If we ever need to restore the InChI sort, these are needed ***
# import autoparse.pattern as app
# import autoparse.find as apf
#
# def _parse_sort_order_from_aux_info(aux_info):
#     ptt = app.escape('/N:') + app.capturing(
#         app.series(app.UNSIGNED_INTEGER, ','))
#     num_str = apf.first_capture(ptt, aux_info)
#     nums = tuple(map(int, num_str.split(',')))
#     return nums
#
# def inchi_with_sort_from_geometry(gra, geo=None, geo_idx_dct=None):
#     """ (if coordinates are passed in, they are used to determine stereo)
#     """
#     mlf, idx_dct = _molfile(gra, geo=geo, geo_idx_dct=geo_idx_dct)
#     rdm = _rdkit.from_molfile(mlf)
#     ich, aux_info = _rdkit.to_inchi(rdm, with_aux_info=True)
#     nums = _parse_sort_order_from_aux_info(aux_info)
#     nums = tuple(map(idx_dct.__getitem__, nums))
#     return ich, nums
