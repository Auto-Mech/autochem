""" graph constructor
"""
from qcelemental import periodictable as pt
from automol import dict_


def from_data(atom_symbols, bond_keys, atom_implicit_hydrogen_valences=None,
              atom_stereo_parities=None, bond_orders=None,
              bond_stereo_parities=None):
    """ construct a molecular graph from data

    format:
        gra = (atm_dct, bnd_dct)
        atm_dct := {atm_key: (atm_sym, atm_imp_hyd_vlc, atm_ste_par), ...}
        bnd_dct := {bnd_key: (bnd_ord, bnd_ste_par), ...}
        [where bnd_key := frozenset({atm1_key, atm2_key})]

    :param atom_symbols: atomic symbols, by atom key
    :type atom_symbols: dict
    :param bond_keys: bond keys
    :type bond_keys: set
    :param atom_implicit_hydrogen_valences: the number of implicit hydrogens
        associated with each atom, by atom key
    :type atom_implicit_hydrogen_valences: dict
    :param atom_stereo_parities: atom stereo parities, by atom key
    :type atom_stereo_parities: dict
    :param bond_orders: bond orders, by bond key
    :type bond_orders: dict
    :param bond_stereo_parities: bond stereo parities, by bond key
    :type bond_stereo_parities: dict
    """
    atm_dct = atoms_from_data(
        atom_symbols=atom_symbols,
        atom_implicit_hydrogen_valences=atom_implicit_hydrogen_valences,
        atom_stereo_parities=atom_stereo_parities)

    bnd_dct = bonds_from_data(
        bond_keys=bond_keys,
        bond_orders=bond_orders,
        bond_stereo_parities=bond_stereo_parities)

    gra = _from_atoms_and_bonds(atoms=atm_dct, bonds=bnd_dct)
    return gra


def atoms_from_data(atom_symbols, atom_implicit_hydrogen_valences=None,
                    atom_stereo_parities=None):
    """ construct atom dictionary from data

    format:
        atm_dct := {atm_key: (atm_sym, atm_imp_hyd_vlc, atm_ste_par), ...}

    :param atom_symbols: atomic symbols, by atom key
    :type atom_symbols: dict
    :param atom_implicit_hydrogen_valences: the number of implicit hydrogens
        associated with each atom, by atom key
    :type atom_implicit_hydrogen_valences: dict
    :param atom_stereo_parities: atom stereo parities, by atom key
    :type atom_stereo_parities: dict
    """
    keys = sorted(atom_symbols.keys())
    syms = dict_.values_by_key(atom_symbols, keys)
    vlcs = dict_.values_by_key(
        dict_.empty_if_none(atom_implicit_hydrogen_valences), keys, fill_val=0)
    pars = dict_.values_by_key(
        dict_.empty_if_none(atom_stereo_parities), keys, fill_val=None)

    natms = len(syms)
    vlcs = [0] * natms if vlcs is None else list(vlcs)
    pars = [None] * natms if pars is None else list(pars)

    assert len(vlcs) == natms
    assert len(pars) == natms

    syms = list(map(pt.to_E, syms))
    vlcs = list(map(int, vlcs))

    assert all(par in (None, False, True) for par in pars)
    pars = [bool(par) if par is not None else par for par in pars]

    atm_dct = dict(zip(keys, zip(syms, vlcs, pars)))
    return atm_dct


def bonds_from_data(bond_keys, bond_orders=None, bond_stereo_parities=None):
    """ construct bond dictionary graph from data

    format:
        bnd_dct := {bnd_key: (bnd_ord, bnd_ste_par), ...}
        [where bnd_key := frozenset({atm1_key, atm2_key})]

    :param bond_keys: bond keys
    :type bond_keys: set
    :param bond_orders: bond orders, by bond key
    :type bond_orders: dict
    :param bond_stereo_parities: bond stereo parities, by bond key
    :type bond_stereo_parities: dict
    """
    keys = sorted(bond_keys)
    assert all(len(key) == 2 for key in keys)
    ords = dict_.values_by_key(
        dict_.empty_if_none(bond_orders), keys, fill_val=1)
    pars = dict_.values_by_key(
        dict_.empty_if_none(bond_stereo_parities), keys, fill_val=None)

    nbnds = len(keys)

    ords = [1] * nbnds if ords is None else list(ords)
    pars = [None] * nbnds if pars is None else list(pars)

    assert len(ords) == nbnds
    assert len(pars) == nbnds

    keys = list(map(frozenset, keys))
    ords = list(map(int, ords))

    assert all(par in (None, False, True) for par in pars)
    pars = [bool(par) if par is not None else par for par in pars]

    bnd_dct = dict(zip(keys, zip(ords, pars)))
    return bnd_dct


def from_atoms_and_bonds(atoms, bonds):
    """ construct a molecular graph from atom and bond dictionaries

    format:
        gra = (atm_dct, bnd_dct)

    :param atoms: atom dictionary
    :type atoms: dict
    :param bonds: bond dictionary
    :type bonds: dict
    """
    atm_dct = dict(atoms)
    bnd_dct = dict(bonds)

    atm_sym_dct = dict_.transform_values(atm_dct, lambda x: x[0])
    atm_imp_hyd_vlc_dct = dict_.transform_values(atm_dct, lambda x: x[1])
    atm_ste_par_dct = dict_.transform_values(atm_dct, lambda x: x[2])

    bnd_ord_dct = dict_.transform_values(bnd_dct, lambda x: x[0])
    bnd_ste_par_dct = dict_.transform_values(bnd_dct, lambda x: x[1])

    return from_data(
        atm_sym_dct, bnd_dct.keys(),
        atom_implicit_hydrogen_valences=atm_imp_hyd_vlc_dct,
        atom_stereo_parities=atm_ste_par_dct, bond_orders=bnd_ord_dct,
        bond_stereo_parities=bnd_ste_par_dct)


def _from_atoms_and_bonds(atoms, bonds):
    """ construct a molecular graph from atom and bond dictionaries

    format:
        gra = (atm_dct, bnd_dct)

    :param atoms: atom dictionary
    :type atoms: dict
    :param bonds: bond dictionary
    :type bonds: dict
    """
    atm_dct = dict(atoms)
    bnd_dct = dict(bonds)
    atm_keys = set(atm_dct.keys())
    bnd_keys = set(bnd_dct.keys())
    assert all(bnd_key <= atm_keys for bnd_key in bnd_keys)
    return (atm_dct, bnd_dct)
