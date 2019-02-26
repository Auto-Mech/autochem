""" graph constructor
"""
import phycon.elements as pce


def from_data(atom_symbols, bond_keys, atom_implicit_hydrogen_valences=None,
              atom_stereo_parities=None, bond_orders=None,
              bond_stereo_parities=None):
    """ construct a molecular graph from data
    """
    atm_dct = _atom_dictionary(atom_symbols, atom_implicit_hydrogen_valences,
                               atom_stereo_parities)
    bnd_dct = _bond_dictionary(atm_dct.keys(), bond_keys, bond_orders,
                               bond_stereo_parities)
    return (atm_dct, bnd_dct)


def _atom_dictionary(syms, vlcs=None, pars=None):
    natms = len(syms)
    vlcs = [0] * natms if vlcs is None else list(vlcs)
    pars = [None] * natms if pars is None else list(pars)

    assert len(vlcs) == natms
    assert len(pars) == natms

    syms = list(map(pce.standard_case, syms))
    vlcs = list(map(int, vlcs))

    assert all(sym in pce.organic_element_keys() for sym in syms)
    assert all(par in (None, False, True) for par in pars)

    atm_dct = dict(enumerate(zip(syms, vlcs, pars)))
    return atm_dct


def _bond_dictionary(atm_keys, keys, ords=None, pars=None):
    nbnds = len(keys)
    atm_keys = set(atm_keys)  # we need atom keys to confirm the bond keys

    ords = [1] * nbnds if ords is None else list(ords)
    pars = [None] * nbnds if pars is None else list(pars)

    assert len(ords) == nbnds
    assert len(pars) == nbnds

    keys = list(map(frozenset, keys))
    ords = list(map(int, ords))

    assert all(len(key) == 2 and key <= atm_keys for key in keys)
    assert all(par in (None, False, True) for par in pars)

    bnd_dct = dict(zip(keys, zip(ords, pars)))
    return bnd_dct
