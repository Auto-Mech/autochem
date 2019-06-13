""" graph conversions
"""
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
    ich = hardcoded_inchi(gra)
    if ich is None:
        if not automol.graph.has_stereo(gra):
            ich = inchi_from_coordinates(gra, atm_xyz_dct=None)
            ich = automol.inchi.recalculate(automol.inchi.without_stereo(ich))
        else:
            gra = automol.graph.explicit(gra)
            atm_xyz_dct = automol.graph.atom_stereo_coordinates(gra)
            ich = inchi_from_coordinates(gra, atm_xyz_dct=atm_xyz_dct)

    return ich


def hardcoded_inchi(gra):
    """ (hardcoded) connectivity graph => inchi conversions
    """
    gra = automol.graph.without_ghost_atoms(gra)
    ich = None
    for _ich, _gra in (
            automol.convert.inchi.HARDCODED_INCHI_TO_GRAPH_DCT.items()):
        if automol.graph.backbone_isomorphic(gra, _gra):
            ich = _ich
    return ich


def inchi_from_coordinates(gra, atm_xyz_dct=None):
    """ connectivity graph => inchi conversion

    (if coordinates are passed in, they are used to determine stereo)
    """
    gra = automol.graph.without_ghost_atoms(gra)
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
    mlf, _ = _molfile.from_data(
        atm_keys, bnd_keys, atm_syms, atm_bnd_vlcs, atm_rad_vlcs, bnd_ords,
        atm_xyzs=atm_xyzs)
    rdm = _rdkit.from_molfile(mlf)
    ich = _rdkit.to_inchi(rdm)
    return ich


def formula(gra):
    """ graph  => formula
    """
    gra = automol.graph.explicit(gra)
    syms = automol.graph.atom_symbols(gra).values()
    fml = _util.formula(syms)
    return fml
