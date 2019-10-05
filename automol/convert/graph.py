""" graph conversions
"""
import autoparse.pattern as app
import autoparse.find as apf
from automol import dict_
import automol.graph
import automol.inchi
import automol.convert.inchi
from automol.convert import _molfile
from automol.convert import _rdkit
from automol.convert import _util


# graph => inchi
def inchi(gra):
    """ graph => inchi
    """
    ich = automol.convert.inchi.object_to_hardcoded_inchi_by_key(
        'graph', gra, comp=_compare)

    if ich is None:
        if not automol.graph.has_stereo(gra):
            ich, _ = inchi_with_sort_from_coordinates(gra, atm_xyz_dct=None)
            ich = automol.inchi.standard_form(ich, remove_stereo=True)
        else:
            gra = automol.graph.explicit(gra)
            atm_xyz_dct = automol.graph.atom_stereo_coordinates(gra)
            ich, _ = inchi_with_sort_from_coordinates(
                gra, atm_xyz_dct=atm_xyz_dct)

    return ich


def _compare(gra1, gra2):
    gra1 = automol.graph.without_dummy_atoms(gra1)
    gra2 = automol.graph.without_dummy_atoms(gra2)
    return automol.graph.backbone_isomorphic(gra1, gra2)


def inchi_with_sort_from_coordinates(gra, atm_xyz_dct=None):
    """ connectivity graph => inchi conversion

    (if coordinates are passed in, they are used to determine stereo)
    """
    gra = automol.graph.without_dummy_atoms(gra)
    gra = automol.graph.dominant_resonance(gra)
    atm_keys = list(automol.graph.atom_keys(gra))
    bnd_keys = list(automol.graph.bond_keys(gra))
    atm_syms = dict_.values_by_key(automol.graph.atom_symbols(gra), atm_keys)
    atm_bnd_vlcs = dict_.values_by_key(
        automol.graph.atom_bond_valences(gra), atm_keys)
    atm_rad_vlcs = dict_.values_by_key(
        automol.graph.atom_unsaturated_valences(gra), atm_keys)
    bnd_ords = dict_.values_by_key(automol.graph.bond_orders(gra), bnd_keys)
    atm_xyzs = (None if atm_xyz_dct is None else
                dict_.values_by_key(atm_xyz_dct, atm_keys))
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
    """ graph => geometry
    """
    gra = automol.graph.explicit(gra)
    atm_keys = sorted(automol.graph.atom_keys(gra))
    atm_syms = dict_.values_by_key(automol.graph.atom_symbols(gra), atm_keys)
    atm_xyzs = dict_.values_by_key(
        automol.graph.atom_stereo_coordinates(gra), atm_keys)
    geo = automol.create.geom.from_data(atm_syms, atm_xyzs)
    return geo


def formula(gra):
    """ graph  => formula
    """
    gra = automol.graph.explicit(gra)
    syms = automol.graph.atom_symbols(gra).values()
    fml = _util.formula(syms)
    return fml
