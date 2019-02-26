""" core functions defining the graph data structure

xgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_sym, atm_imp_hyd_vlc, atm_ste_par), ...}
bnd_dct: {bnd_key: (bnd_ord, bnd_ste_par), ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
import numpy
from ..constructors.graph import from_data as _from_data
from ._dict import transform_keys as _transform_keys
from ._dict import by_key as _by_key
from ._dict import values_by_key as _values_by_key
from ._dict.multi import by_key_by_position as _by_key_by_position
from ._dict.multi import set_by_key_by_position as _set_by_key_by_position

ATM_SYM_POS = 0
ATM_IMP_HYD_VLC_POS = 1
ATM_STE_PAR_POS = 2

BND_ORD_POS = 0
BND_STE_PAR_POS = 1


# constructors
def from_atoms_and_bonds(atm_dct, bnd_dct):
    """ molecular graph from atoms and bonds
    """
    return (atm_dct, bnd_dct)


def from_dictionaries(atm_sym_dct, bnd_keys, atm_imp_hyd_vlc_dct=None,
                      atm_ste_par_dct=None, bnd_ord_dct=None,
                      bnd_ste_par_dct=None):
    """ molecular graph from dictionaries over atom and bond keys
    """
    atm_keys = sorted(atm_sym_dct.keys())

    def _values(dct, keys, fill_val=None):
        dct = dict() if dct is None else dct
        return _values_by_key(dct, keys, fill_val=fill_val)

    xgr = _from_data(
        atom_symbols=_values(atm_sym_dct, atm_keys),
        bond_keys=[set(map(atm_keys.index, bnd_key)) for bnd_key in bnd_keys],
        atom_implicit_hydrogen_valences=_values(
            atm_imp_hyd_vlc_dct, atm_keys, fill_val=0),
        atom_stereo_parities=_values(atm_ste_par_dct, atm_keys, fill_val=None),
        bond_orders=_values(bnd_ord_dct, bnd_keys, fill_val=1),
        bond_stereo_parities=_values(bnd_ste_par_dct, bnd_keys, fill_val=None)
    )
    xgr = relabel(xgr, dict(enumerate(atm_keys)))
    return xgr


def add_atoms(xgr, sym_dct, imp_hyd_vlc_dct=None, ste_par_dct=None):
    """ add atoms to this molecular graph
    """
    atm_keys = atom_keys(xgr)
    atm_sym_dct = atom_symbols(xgr)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(xgr)
    atm_ste_par_dct = atom_stereo_parities(xgr)
    bnd_keys = bond_keys(xgr)
    bnd_ord_dct = bond_orders(xgr)
    bnd_ste_par_dct = bond_stereo_parities(xgr)

    keys = set(sym_dct.keys())
    imp_hyd_vlc_dct = {} if imp_hyd_vlc_dct is None else imp_hyd_vlc_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct

    assert not keys & atm_keys
    assert set(imp_hyd_vlc_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys

    atm_sym_dct.update(sym_dct)
    atm_imp_hyd_vlc_dct.update(imp_hyd_vlc_dct)
    atm_ste_par_dct.update(ste_par_dct)

    xgr = from_dictionaries(atm_sym_dct, bnd_keys, atm_imp_hyd_vlc_dct,
                            atm_ste_par_dct, bnd_ord_dct, bnd_ste_par_dct)
    return xgr


def add_bonds(xgr, keys, ord_dct=None, ste_par_dct=None):
    """ add bonds to this molecular graph
    """
    atm_sym_dct = atom_symbols(xgr)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(xgr)
    atm_ste_par_dct = atom_stereo_parities(xgr)
    bnd_keys = set(bond_keys(xgr))
    bnd_ord_dct = bond_orders(xgr)
    bnd_ste_par_dct = bond_stereo_parities(xgr)

    keys = set(keys)
    ord_dct = {} if ord_dct is None else ord_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct

    assert not keys & bnd_keys
    assert set(ord_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys

    bnd_keys.update(keys)
    bnd_ord_dct.update(ord_dct)
    bnd_ste_par_dct.update(ste_par_dct)

    xgr = from_dictionaries(atm_sym_dct, bnd_keys, atm_imp_hyd_vlc_dct,
                            atm_ste_par_dct, bnd_ord_dct, bnd_ste_par_dct)
    return xgr


def frozen(xgr):
    """ hashable, sortable, immutable container of graph data
    """
    atm_keys = sorted(atom_keys(xgr))
    bnd_keys = sorted(bond_keys(xgr), key=sorted)

    # make it sortable by replacing Nones with -infinity
    atm_vals = numpy.array(_values_by_key(atoms(xgr), atm_keys))
    bnd_vals = numpy.array(_values_by_key(bonds(xgr), bnd_keys))
    atm_vals[numpy.equal(atm_vals, None)] = -numpy.inf
    bnd_vals[numpy.equal(bnd_vals, None)] = -numpy.inf

    frz_atms = tuple(zip(atm_keys, map(tuple, atm_vals)))
    frz_bnds = tuple(zip(bnd_keys, map(tuple, bnd_vals)))
    return (frz_atms, frz_bnds)


# value getters
def atoms(xgr):
    """ atoms, as a dictionary
    """
    atm_dct, _ = xgr
    return atm_dct


def bonds(xgr):
    """ bonds, as a dictionary
    """
    _, bnd_dct = xgr
    return bnd_dct


def atom_keys(xgr):
    """ sorted atom keys
    """
    return frozenset(atoms(xgr).keys())


def bond_keys(xgr):
    """ sorted bond keys
    """
    return frozenset(bonds(xgr).keys())


def atom_symbols(xgr):
    """ atom symbols, as a dictionary
    """
    return _by_key_by_position(atoms(xgr), atom_keys(xgr), ATM_SYM_POS)


def atom_implicit_hydrogen_valences(xgr):
    """ atom implicit hydrogen valences, as a dictionary
    """
    return _by_key_by_position(atoms(xgr), atom_keys(xgr),
                               ATM_IMP_HYD_VLC_POS)


def atom_stereo_parities(sgr):
    """ atom parities, as a dictionary
    """
    return _by_key_by_position(atoms(sgr), atom_keys(sgr), ATM_STE_PAR_POS)


def bond_orders(rgr):
    """ bond orders, as a dictionary
    """
    return _by_key_by_position(bonds(rgr), bond_keys(rgr), BND_ORD_POS)


def bond_stereo_parities(sgr):
    """ bond parities, as a dictionary
    """
    return _by_key_by_position(bonds(sgr), bond_keys(sgr), BND_STE_PAR_POS)


# value setters
def relabel(xgr, atm_key_dct):
    """ relabel the graph with new atom keys
    """
    orig_atm_keys = atom_keys(xgr)
    assert set(atm_key_dct.keys()) <= orig_atm_keys

    new_atm_key_dct = dict(zip(orig_atm_keys, orig_atm_keys))
    new_atm_key_dct.update(atm_key_dct)

    _relabel_atom_key = new_atm_key_dct.__getitem__

    def _relabel_bond_key(bnd_key):
        return frozenset(map(_relabel_atom_key, bnd_key))

    atm_dct = _transform_keys(atoms(xgr), _relabel_atom_key)
    bnd_dct = _transform_keys(bonds(xgr), _relabel_bond_key)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def set_atom_implicit_hydrogen_valences(xgr, atm_imp_hyd_vlc_dct):
    """ set atom implicit hydrogen valences
    """
    atm_dct = _set_by_key_by_position(atoms(xgr), atm_imp_hyd_vlc_dct,
                                      ATM_IMP_HYD_VLC_POS)
    bnd_dct = bonds(xgr)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def set_atom_stereo_parities(sgr, atm_par_dct):
    """ set atom parities
    """
    atm_dct = _set_by_key_by_position(atoms(sgr), atm_par_dct, ATM_STE_PAR_POS)
    return from_atoms_and_bonds(atm_dct, bonds(sgr))


def set_bond_orders(rgr, bnd_ord_dct):
    """ set bond orders
    """
    bnd_dct = _set_by_key_by_position(bonds(rgr), bnd_ord_dct, BND_ORD_POS)
    return from_atoms_and_bonds(atoms(rgr), bnd_dct)


def set_bond_stereo_parities(sgr, bnd_par_dct):
    """ set bond parities
    """
    bnd_dct = _set_by_key_by_position(bonds(sgr), bnd_par_dct, BND_STE_PAR_POS)
    return from_atoms_and_bonds(atoms(sgr), bnd_dct)


# transformations
def without_bond_orders(xgr):
    """ resonance graph with maximum spin (i.e. no pi bonds)
    """
    bnd_ord_dct = _by_key({}, bond_keys(xgr), fill_val=1)
    return set_bond_orders(xgr, bnd_ord_dct)


def without_stereo_parities(xgr):
    """ graph with stereo assignments wiped out
    """
    atm_ste_par_dct = _by_key({}, atom_keys(xgr), fill_val=None)
    bnd_ste_par_dct = _by_key({}, bond_keys(xgr), fill_val=None)
    xgr = set_atom_stereo_parities(xgr, atm_ste_par_dct)
    xgr = set_bond_stereo_parities(xgr, bnd_ste_par_dct)
    return xgr
