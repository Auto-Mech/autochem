""" InChI string conversions
"""
from ._inchi import key_layer as _key_layer
from ._inchi import core_parent as _core_parent
from ._inchi import atom_stereo_elements as _atom_stereo_elements
from ._inchi import bond_stereo_elements as _bond_stereo_elements
from ._inchi import known_atom_stereo_elements as _known_atom_stereo_elements
from ._inchi import known_bond_stereo_elements as _known_bond_stereo_elements
from ._inchi import has_unknown_stereo_elements as _has_unknown_stereo_elements
from ._rdkit import from_inchi as _rdm_from_inchi
from ._rdkit import to_inchi as _rdm_to_inchi
from ._rdkit import inchi_to_inchi_key as _inchi_to_inchi_key
from ._rdkit import to_smiles as _rdm_to_smiles
from ._rdkit import geometry as _rdm_to_geometry
from ._rdkit import connectivity_graph as _rdm_to_connectivity_graph
from ._pybel import from_inchi as _pbm_from_inchi
from ._pybel import geometry as _pbm_to_geometry
from ..geom import inchi as _inchi_from_geometry
from ..graph import inchi as _inchi_from_graph
from ..graph import stereo_inchi as _stereo_inchi_from_stereo_graph


def inchi_key(ich):
    """ computes InChIKey from an InChI string
    """
    return _inchi_to_inchi_key(ich)


def smiles(ich):
    """ SMILES string from an InChI string
    """
    rdm = _rdm_from_inchi(ich)
    smi = _rdm_to_smiles(rdm)
    return smi


def connectivity_graph(ich):
    """ connectivity graph from an InChI string
    """
    rdm = _rdm_from_inchi(ich)

    # make sure the InChI string was valid
    ich_, _ = _rdm_to_inchi(rdm, with_aux_info=True)
    assert _core_parent(ich) == _core_parent(ich_)

    cgr = _rdm_to_connectivity_graph(rdm)
    cgr_ich = _inchi_from_graph(cgr)
    assert _has_same_connectivity(ich, cgr_ich)
    return cgr


def stereo_graph(ich):
    """ stereo graph from an InChI string
    """
    def _int_minus_one(int_str):
        return int(int_str) - 1

    def _atom_key(ich_atm_key):
        return _int_minus_one(ich_atm_key)

    def _bond_key(ich_bnd_key):
        return frozenset(map(_int_minus_one, str.split(ich_bnd_key, '-')))

    def _value(ich_ste_val):
        assert ich_ste_val in ('-', '+')
        return ich_ste_val == '+'

    atms, cnns = connectivity_graph(ich)
    assert not _has_unknown_stereo_elements(ich)
    atm_ste_dct = {_atom_key(key): _value(val)
                   for key, val in _atom_stereo_elements(ich)}
    bnd_ste_dct = {_bond_key(key): (1, _value(val))
                   for key, val in _bond_stereo_elements(ich)}
    assert set(atm_ste_dct.keys()) <= set(range(len(atms)))
    assert set(bnd_ste_dct.keys()) <= set(cnns.keys())
    atms = {atm_key: ((sym, hcnt, atm_ste_dct[atm_key])
                      if atm_key in atm_ste_dct else (sym, hcnt, None))
            for atm_key, (sym, hcnt, _) in atms.items()}
    bnds = cnns.copy()
    bnds.update(bnd_ste_dct)
    sgr = (atms, bnds)
    sgr_ich = _stereo_inchi_from_stereo_graph(sgr)
    assert _has_same_connectivity(ich, sgr_ich)
    assert _has_compatible_stereo(ich, sgr_ich)
    return sgr


def geometry(ich):
    """ cartesian geometry from an InChI string
    """
    try:
        rdm = _rdm_from_inchi(ich)
        geo = _rdm_to_geometry(rdm)
        geo_ich = _inchi_from_geometry(geo, stereo=True)
        assert _has_same_connectivity(ich, geo_ich)
        assert _has_compatible_stereo(ich, geo_ich)
    except (AssertionError, RuntimeError):
        pbm = _pbm_from_inchi(ich)
        geo = _pbm_to_geometry(pbm)
        geo_ich = _inchi_from_geometry(geo, stereo=True)
        assert _has_same_connectivity(ich, geo_ich)
        assert _has_compatible_stereo(ich, geo_ich)
    return geo


def _has_same_connectivity(ich, other_ich):
    """ do these InChI strings have the same connectivity?
    """
    return (_key_layer(ich, 'c') == _key_layer(other_ich, 'c') and
            _key_layer(ich, 'h') == _key_layer(other_ich, 'h'))


def _has_compatible_stereo(ich, other_ich):
    """ is `other_ich` compatible with `ich`?
    """
    return (set(_known_atom_stereo_elements(ich)) <=
            set(_known_atom_stereo_elements(other_ich)) and
            set(_known_bond_stereo_elements(ich)) <=
            set(_known_bond_stereo_elements(other_ich)))
