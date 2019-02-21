""" graph -> InChI string conversion
"""
from itertools import chain as _chain
from ._molfile import from_data as _mlf_from_data
from ._inchi_aux import sorted_atom_keys as _ich_aux_sorted_atom_keys
from ._rdkit import from_molfile as _rdm_from_molfile
from ._rdkit import to_inchi_with_aux_info as _rdm_to_inchi_with_aux_info
from .._core import atom_keys
from .._core import bond_keys
from .._core import atom_symbols
from .._core import bond_orders
from .._expl import atom_explicit_hydrogen_keys
from .._expl import backbone_keys
from .._expl import backbone_isomorphic
from .._expl import backbone_isomorphism
from .._res import atom_bond_valences
from .._res import atom_radical_valences
from .._res import dominant_resonance
from .._dict import values_by_key as _values_by_key
from .._dict import keys_sorted_by_value as _keys_sorted_by_value


def atom_inchi_numbers(xgr):
    """ InChI numbers, by atom
    """
    _, atm_ich_num_dct = _with_atom_inchi_numbers(xgr)
    return atm_ich_num_dct


def inchi(xgr):
    """ InChI string of this connectivity graph
    """
    ich, _ = _with_atom_inchi_numbers(xgr)
    return ich


# def atom_stereo_inchi_numbers_from_coordinates(xgr, atm_xyz_dct):
#     """ stereo-specific InChI numbers, by atom
#     """
#     _, atm_ich_num_dct = _with_atom_inchi_numbers(xgr, atm_xyz_dct)
#     return atm_ich_num_dct


def stereo_inchi_from_coordinates(xgr, atm_xyz_dct):
    """ stereo-specific InChI of this connectivity graph from coordinates
    """
    ich, _ = _with_atom_inchi_numbers(xgr, atm_xyz_dct)
    return ich


def _with_atom_inchi_numbers(xgr, atm_xyz_dct=None):
    """ InChI string with numbering from a connectivity graph

    For stereo InChIs, pass cartesian coordinates and set the chirality flag.
    """
    ich, bbn_ich_num_dct = _catch_hardcoded(xgr)
    if ich is None:
        mlf, mlf_atm_key_dct = _molfile_with_atom_key_mapping(xgr, atm_xyz_dct)
        rdm = _rdm_from_molfile(mlf)
        ich, ich_aux = _rdm_to_inchi_with_aux_info(rdm)

        # determine the inchi numbering from the AuxInfo string
        ich_srt_mlf_bbn_keys = _ich_aux_sorted_atom_keys(ich_aux)
        ich_srt_bbn_keys = _values_by_key(mlf_atm_key_dct,
                                          ich_srt_mlf_bbn_keys)
        assert set(ich_srt_bbn_keys) == set(backbone_keys(xgr))
        bbn_ich_num_dct = dict(map(reversed, enumerate(ich_srt_bbn_keys)))

    atm_ich_num_dct = _fill_atom_inchi_numbers(xgr, bbn_ich_num_dct)
    return ich, atm_ich_num_dct


def _catch_hardcoded(xgr):
    """ hardcoded molecules with more than 2 unpaired electrons
    """
    ich = bbn_ich_num_dct = None

    graph_dct = {
        'InChI=1S/C': ({0: ('C', 0, None)}, {}),
        'InChI=1S/N': ({0: ('N', 0, None)}, {}),
        'InChI=1S/CH/h1H': ({0: ('C', 1, None)}, {}),
        'InChI=1S/CF/c1-2': ({0: ('C', 0, None), 1: ('F', 0, None)},
                             {frozenset({0, 1}): (1, None)}),
        'InChI=1S/CCl/c1-2': ({0: ('C', 0, None), 1: ('Cl', 0, None)},
                              {frozenset({0, 1}): (1, None)}),
    }
    for ref_ich, ref_xgr in graph_dct.items():
        if backbone_isomorphic(xgr, ref_xgr):
            ich = ref_ich
            bbn_ich_num_dct = backbone_isomorphism(xgr, ref_xgr)

    return ich, bbn_ich_num_dct


def _molfile_with_atom_key_mapping(xgr, atm_xyz_dct=None):
    rgr = dominant_resonance(xgr)
    atm_keys = list(atom_keys(rgr))
    bnd_keys = list(bond_keys(rgr))
    atm_syms = _values_by_key(atom_symbols(rgr), atm_keys)
    atm_bnd_vlcs = _values_by_key(atom_bond_valences(rgr), atm_keys)
    atm_rad_vlcs = _values_by_key(atom_radical_valences(rgr), atm_keys)
    atm_xyzs = (None if atm_xyz_dct is None else
                _values_by_key(atm_xyz_dct, atm_keys))
    bnd_ords = _values_by_key(bond_orders(rgr), bnd_keys)
    mlf, mlf_atm_key_dct = _mlf_from_data(atm_keys, bnd_keys, atm_syms,
                                          atm_bnd_vlcs, atm_rad_vlcs, bnd_ords,
                                          atm_xyzs=atm_xyzs)
    return mlf, mlf_atm_key_dct


def _fill_atom_inchi_numbers(xgr, bbn_ich_num_dct):
    """ atom inchi number dictionary from inchi-sorted backbone keys
    """
    atm_ich_num_dct = bbn_ich_num_dct.copy()

    ich_srt_bbn_keys = _keys_sorted_by_value(bbn_ich_num_dct)
    atm_exp_hyd_keys_dct = atom_explicit_hydrogen_keys(xgr)
    ich_srt_bbn_exp_hyd_keys = _values_by_key(atm_exp_hyd_keys_dct,
                                              ich_srt_bbn_keys)
    ich_srt_exp_hyd_keys = tuple(_chain(*ich_srt_bbn_exp_hyd_keys))
    first_exp_hyd_ich_num = min(bbn_ich_num_dct.values()) - 1
    atm_ich_num_dct.update({
        exp_hyd_key: first_exp_hyd_ich_num - exp_hyd_srt_idx
        for exp_hyd_srt_idx, exp_hyd_key in enumerate(ich_srt_exp_hyd_keys)
    })
    assert set(atm_ich_num_dct.keys()) == atom_keys(xgr)
    return atm_ich_num_dct
