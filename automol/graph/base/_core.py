""" core graph functions

Data structure:
    gra = (atm_dct, bnd_dct)
    atm_dct := {
        atm_key: (symb, imp_hyd, ste_par),
        ...
    }
    bnd_dct := {
        bnd_key: (ord, ste_par),
        ...
    }
    [where bnd_key := frozenset({atm1_key, atm2_key})]

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import numbers
import itertools
import functools
import numpy
import yaml
from phydat import ptab
from phydat import phycon
import automol.formula
from automol import util
from automol.util import dict_
import automol.util.dict_.multi as mdict

ATM_SYM_POS = 0
ATM_IMP_HYD_POS = 1
ATM_STE_PAR_POS = 2
TS_ATM_PRD_STE_PAR_POS = 3
TS_ATM_FLE_STE_PAR_POS = 4

BND_ORD_POS = 0
BND_STE_PAR_POS = 1
TS_BND_PRD_STE_PAR_POS = 2
TS_BND_FLE_STE_PAR_POS = 3

ATM_PROP_NAMES = ('symbol', 'implicit_hydrogens', 'stereo_parity',
                  'prod_stereo_parity', 'ts_stereo_parity')
BND_PROP_NAMES = ('order', 'stereo_parity',
                  'prod_stereo_parity', 'ts_stereo_parity')


# # constructors
def from_data(atm_symb_dct, bnd_keys, atm_imp_hyd_dct=None,
              atm_ste_par_dct=None,
              atm_prd_ste_par_dct=None, atm_ts_ste_par_dct=None,
              bnd_ord_dct=None, bnd_ste_par_dct=None,
              bnd_prd_ste_par_dct=None, bnd_ts_ste_par_dct=None,
              ts_=None):
    """ Construct a molecular graph from data

    Ordinary graph data structure:
        gra = (atm_dct, bnd_dct)
        atm_dct := {
            atm_key: (symb, imp_hyd, ste_par),
            ...
        }
        bnd_dct := {
            bnd_key: (ord, ste_par),
            ...
        }
        [where bnd_key := frozenset({atm1_key, atm2_key})]

    TS graph data structure:
        gra = (atm_dct, bnd_dct)
        atm_dct := {
            atm_key: (symb, imp_hyd, ste_par, prd_ste_par, ts_ste_par),
            ...
        }
        bnd_dct := {
            bnd_key: (ord, ste_par, prd_ste_par, ts_ste_par),
            ...
        }
        [where bnd_key := frozenset({atm1_key, atm2_key})]

    :param atm_symb_dct: atomic symbols, by atom key
    :type atm_symb_dct: dict
    :param bnd_keys: bond keys
    :type bnd_keys: set
    :param atm_imp_hyd_dct: the number of implicit hydrogens associated
        with each atom, by atom key
    :type atm_imp_hyd_dct: dict
    :param atm_ste_par_dct: stereo parities, by atom key; (in TS graphs,
        reactant stereo parities)
    :type atm_ste_par_dct: dict
    :param atm_prd_ste_par_dct: product stereo parities, by atom key; (TS
        graphs only)
    :type atm_prd_ste_par_dct: dict
    :param atm_ts_ste_par_dct: fleeting TS stereo parities, by atom key; (TS
        graphs only)
    :type atm_ts_ste_par_dct: dict
    :param bnd_ord_dct: bond orders, by bond key
    :type bnd_ord_dct: dict
    :param bnd_ste_par_dct: stereo parities, by bond key; (in TS graphs,
        reactant stereo parities)
    :type bnd_ste_par_dct: dict
    :param bnd_prd_ste_par_dct: product stereo parities, by bond key; (TS
        graphs only)
    :type bnd_prd_ste_par_dct: dict
    :param bnd_ts_ste_par_dct: fleeting TS stereo parities, by bond key; (TS
        graphs only)
    :type bnd_ts_ste_par_dct: dict
    :param ts_: Create a TS graph?
    :type ts_: bool
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    # If `ts_` is `None`, infer whether this is TS or non-TS from the data
    # provided.
    if ts_ is None:
        if (atm_prd_ste_par_dct is not None and
                atm_ts_ste_par_dct is not None and
                bnd_prd_ste_par_dct is not None and
                bnd_ts_ste_par_dct is not None):
            # It has the right format for a TS graph, so set `ts_` to `True`
            ts_ = True
        else:
            # Require that it has the right format for a non-TS graph and set
            # `ts_` to `False`
            assert atm_prd_ste_par_dct is None, (
                f"Cannot add atom product stereo parities to non-TS graph!\n"
                f"atm_prd_ste_par_dct={atm_prd_ste_par_dct}")
            assert atm_ts_ste_par_dct is None, (
                f"Cannot add atom TS stereo parities to non-TS graph!\n"
                f"atm_ts_ste_par_dct={atm_ts_ste_par_dct}")
            assert bnd_prd_ste_par_dct is None, (
                f"Cannot add bond product stereo parities to non-TS graph!\n"
                f"bnd_prd_ste_par_dct={bnd_prd_ste_par_dct}")
            assert bnd_ts_ste_par_dct is None, (
                f"Cannot add bond TS stereo parities to non-TS graph!\n"
                f"bnd_ts_ste_par_dct={bnd_ts_ste_par_dct}")
            ts_ = False

    atm_dct = atoms_from_data(
        atm_symb_dct=atm_symb_dct,
        atm_imp_hyd_dct=atm_imp_hyd_dct,
        atm_ste_par_dct=atm_ste_par_dct,
        atm_prd_ste_par_dct=atm_prd_ste_par_dct,
        atm_ts_ste_par_dct=atm_ts_ste_par_dct,
        ts_=ts_)

    bnd_dct = bonds_from_data(
        bnd_keys=bnd_keys,
        bnd_ord_dct=bnd_ord_dct,
        bnd_ste_par_dct=bnd_ste_par_dct,
        bnd_prd_ste_par_dct=bnd_prd_ste_par_dct,
        bnd_ts_ste_par_dct=bnd_ts_ste_par_dct,
        ts_=ts_)

    atm_keys = set(atm_dct.keys())
    bnd_keys = set(bnd_dct.keys())

    assert all(bnd_key <= atm_keys for bnd_key in bnd_keys), (
        f"\natm_keys: {atm_keys}\nbnd_keys: {bnd_keys}")

    return (atm_dct, bnd_dct)


def atoms_from_data(atm_symb_dct, atm_imp_hyd_dct=None, atm_ste_par_dct=None,
                    atm_prd_ste_par_dct=None, atm_ts_ste_par_dct=None,
                    ts_=False):
    """ Construct an atom dictionary from constituent data.

    Ordinary graph data structure:
        atm_dct := {
            atm_key: (symb, imp_hyd, ste_par),
            ...
        }

    TS graph data structure:
        atm_dct := {
            atm_key: (symb, imp_hyd, ste_par, prd_ste_par, ts_ste_par),
            ...
        }

    :param atm_symb_dct: atomic symbols, by atom key
    :type atm_symb_dct: dict
    :param atm_imp_hyd_dct: the number of implicit hydrogens associated
        with each atom, by atom key
    :type atm_imp_hyd_dct: dict
    :param atm_ste_par_dct: stereo parities, by atom key; (in TS graphs,
        reactant stereo parities)
    :type atm_ste_par_dct: dict
    :param atm_prd_ste_par_dct: product stereo parities, by atom key; (TS
        graphs only)
    :type atm_prd_ste_par_dct: dict
    :param atm_ts_ste_par_dct: fleeting TS stereo parities, by atom key; (TS
        graphs only)
    :type atm_ts_ste_par_dct: dict
    :param ts_: Create a TS graph?
    :type ts_: bool
    :returns: The atoms and their associated properties, as a dictionary
    :rtype: dict[int: tuple]
    """
    keys = sorted(atm_symb_dct.keys())
    symbs = dict_.values_by_key(atm_symb_dct, keys)
    hyds = dict_.values_by_key(
        dict_.empty_if_none(atm_imp_hyd_dct), keys, fill_val=0)
    pars = dict_.values_by_key(
        dict_.empty_if_none(atm_ste_par_dct), keys, fill_val=None)
    if ts_:
        prd_pars = dict_.values_by_key(
            dict_.empty_if_none(atm_prd_ste_par_dct), keys, fill_val=None)
        ts_pars = dict_.values_by_key(
            dict_.empty_if_none(atm_ts_ste_par_dct), keys, fill_val=None)

    natms = len(keys)

    assert len(symbs) == natms
    assert len(hyds) == natms
    assert len(pars) == natms
    assert all(par in (None, False, True) for par in pars)
    if ts_:
        assert len(prd_pars) == natms
        assert all(par in (None, False, True) for par in prd_pars)
        assert len(ts_pars) == natms
        assert all(par in (None, False, True) for par in ts_pars)

    symbs = list(map(ptab.to_symbol, symbs))
    hyds = list(map(int, hyds))
    pars = [bool(par) if par is not None else par for par in pars]
    if ts_:
        prd_pars = [bool(par) if par is not None else par for par in prd_pars]
        ts_pars = [bool(par) if par is not None else par for par in ts_pars]

    if not ts_:
        atm_dct = dict(zip(keys, zip(symbs, hyds, pars)))
    else:
        atm_dct = dict(zip(keys, zip(symbs, hyds, pars, prd_pars, ts_pars)))

    return atm_dct


def bonds_from_data(bnd_keys, bnd_ord_dct=None, bnd_ste_par_dct=None,
                    bnd_prd_ste_par_dct=None, bnd_ts_ste_par_dct=None,
                    ts_=False):
    """ Construct a bond dictionary from constituent data.

    Ordinary graph data structure:
        bnd_dct := {
            bnd_key: (ord, ste_par),
            ...
        }
        [where bnd_key := frozenset({atm1_key, atm2_key})]

    TS graph data structure:
        bnd_dct := {
            bnd_key: (ord, ste_par, prd_ste_par, ts_ste_par),
            ...
        }
        [where bnd_key := frozenset({atm1_key, atm2_key})]

    :param bnd_keys: bond keys
    :type bnd_keys: set
    :param bnd_ord_dct: bond orders, by bond key
    :type bnd_ord_dct: dict
    :param bnd_ste_par_dct: stereo parities, by bond key; (in TS graphs,
        reactant stereo parities)
    :type bnd_ste_par_dct: dict
    :param bnd_prd_ste_par_dct: product stereo parities, by bond key; (TS
        graphs only)
    :type bnd_prd_ste_par_dct: dict
    :param bnd_ts_ste_par_dct: fleeting TS stereo parities, by bond key; (TS
        graphs only)
    :type bnd_ts_ste_par_dct: dict
    :param ts_: Create a TS graph?
    :type ts_: bool
    :returns: The bonds and their associated properties, as a dictionary
    :rtype: dict[frozenset({int, int}): tuple]
    """
    keys = sorted(bnd_keys)
    assert all(len(key) == 2 for key in keys)
    ords = dict_.values_by_key(
        dict_.empty_if_none(bnd_ord_dct), keys, fill_val=1)
    pars = dict_.values_by_key(
        dict_.empty_if_none(bnd_ste_par_dct), keys, fill_val=None)
    if ts_:
        prd_pars = dict_.values_by_key(
            dict_.empty_if_none(bnd_prd_ste_par_dct), keys, fill_val=None)
        ts_pars = dict_.values_by_key(
            dict_.empty_if_none(bnd_ts_ste_par_dct), keys, fill_val=None)

    nbnds = len(keys)

    assert len(ords) == nbnds
    assert len(pars) == nbnds
    assert all(par in (None, False, True) for par in pars)
    if ts_:
        assert len(prd_pars) == nbnds
        assert all(par in (None, False, True) for par in prd_pars)
        assert len(ts_pars) == nbnds
        assert all(par in (None, False, True) for par in ts_pars)

    keys = list(map(frozenset, keys))
    # If ts_ = False, we should assert that round(o) == o
    ords = [int(o) if round(o) == o else float(round(o, 1)) for o in ords]
    pars = [bool(par) if par is not None else par for par in pars]
    if ts_:
        prd_pars = [bool(par) if par is not None else par for par in prd_pars]
        ts_pars = [bool(par) if par is not None else par for par in ts_pars]

    if not ts_:
        bnd_dct = dict(zip(keys, zip(ords, pars)))
    else:
        bnd_dct = dict(zip(keys, zip(ords, pars, prd_pars, ts_pars)))

    return bnd_dct


def from_atoms_and_bonds(atm_dct, bnd_dct, ts_=None):
    """ Construct a molecular graph from atom and bond dictionaries.

    Data structure:
        gra = (atm_dct, bnd_dct)

    :param atm_dct: atom dictionary
    :type atm_dct: dict[int: tuple]
    :param bnd_dct: bond dictionary
    :type bnd_dct: dict[frozenset({int, int}): tuple]
    :param ts_: Create a TS graph? If `None`, the answer will be inferred
        from the atom and bond dictionaries.
    :type ts_: bool
    :returns: A molecular graph
    :rtype: automol graph data structure
    """

    atm_dct = dict(atm_dct)
    bnd_dct = dict(bnd_dct)

    atm_nprops = dict_.transform_values(atm_dct, len).values()
    bnd_nprops = dict_.transform_values(bnd_dct, len).values()

    if (all(n == 5 for n in atm_nprops) and
            all(n == 4 for n in bnd_nprops)):
        # It has the right format for a TS graph, so set `ts_` to `True`
        data_has_ts_format = True
        ts_ = True if ts_ is None else ts_
    else:
        # Require that it has the right format for a non-TS graph and set
        # `ts_` to `False`
        assert (
            all(n == 3 for n in atm_nprops) and
            all(n == 2 for n in bnd_nprops)
        ), (
            f"Atom or bond dictionary has improper format for non-TS"
            f"graph:\natm_dct:\n{atm_dct}\nbnd_dct:\n{bnd_dct}\n"
        )
        data_has_ts_format = False
        ts_ = False if ts_ is None else ts_

    atm_nprops = dict_.transform_values(atm_dct, len)
    bnd_nprops = dict_.transform_values(bnd_dct, len)

    atm_keys = sorted(atm_dct.keys())
    atm_symb_dct = (
        mdict.by_key_by_position(atm_dct, atm_keys, ATM_SYM_POS))
    atm_imp_hyd_dct = (
        mdict.by_key_by_position(atm_dct, atm_keys, ATM_IMP_HYD_POS))
    atm_ste_par_dct = (
        mdict.by_key_by_position(atm_dct, atm_keys, ATM_STE_PAR_POS))

    atm_prd_ste_par_dct = (
        mdict.by_key_by_position(atm_dct, atm_keys, TS_ATM_PRD_STE_PAR_POS)
        if ts_ and data_has_ts_format else None)
    atm_ts_ste_par_dct = (
        mdict.by_key_by_position(atm_dct, atm_keys, TS_ATM_FLE_STE_PAR_POS)
        if ts_ and data_has_ts_format else None)

    bnd_keys = sorted(bnd_dct.keys())
    bnd_ord_dct = (
        mdict.by_key_by_position(bnd_dct, bnd_keys, BND_ORD_POS))
    bnd_ste_par_dct = (
        mdict.by_key_by_position(bnd_dct, bnd_keys, BND_STE_PAR_POS))
    bnd_prd_ste_par_dct = (
        mdict.by_key_by_position(bnd_dct, bnd_keys, TS_BND_PRD_STE_PAR_POS)
        if ts_ and data_has_ts_format else None)
    bnd_ts_ste_par_dct = (
        mdict.by_key_by_position(bnd_dct, bnd_keys, TS_BND_FLE_STE_PAR_POS)
        if ts_ and data_has_ts_format else None)

    return from_data(
        atm_symb_dct, bnd_dct.keys(),
        atm_imp_hyd_dct=atm_imp_hyd_dct,
        atm_ste_par_dct=atm_ste_par_dct,
        atm_prd_ste_par_dct=atm_prd_ste_par_dct,
        atm_ts_ste_par_dct=atm_ts_ste_par_dct,
        bnd_ord_dct=bnd_ord_dct,
        bnd_ste_par_dct=bnd_ste_par_dct,
        bnd_prd_ste_par_dct=bnd_prd_ste_par_dct,
        bnd_ts_ste_par_dct=bnd_ts_ste_par_dct, ts_=ts_)


# # getters
def atoms(gra):
    """ Get the atoms of this graph, along with their associated properties

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atoms and their associated properties, as a dictionary
    :rtype: dict[int: tuple]
    """
    atm_dct, _ = gra
    return atm_dct


def bonds(gra):
    """ Get the bonds of this graph, along with their associated properties

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The bonds and their associated properties, as a dictionary
    :rtype: dict[frozenset({int, int}): tuple]
    """
    _, bnd_dct = gra
    return bnd_dct


def atom_keys(gra, symb=None, excl_symbs=()):
    """ Get the atom keys of this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symb: Optionally, restrict this to atoms with a particular
        atomic symbol (e.g., 'H' for hydrogens).
    :type symb: str
    :param excl_symbs: Optionally, exclude atoms with particular atomic
        symbols.
    :type excl_symbs: tuple[str]
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    atm_keys = frozenset(atoms(gra).keys())
    if symb is not None:
        atm_sym_dct = atom_symbols(gra)
        atm_keys = frozenset(k for k in atm_keys
                             if atm_sym_dct[k] == symb and
                             atm_sym_dct[k] not in excl_symbs)
    return atm_keys


def bond_keys(gra, ts_=True):
    """ Get the bond keys of this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: The bond keys
    :rtype: frozenset[{int, int}]
    """
    gra = gra if ts_ else ts_reactants_graph(gra)
    return frozenset(bonds(gra).keys())


def atom_symbols(gra):
    """ Get the atom symbols of this molecular graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of atomic symbols, by atom key
    :rtype: dict[int: str]
    """
    return mdict.by_key_by_position(atoms(gra), atom_keys(gra), ATM_SYM_POS)


def bond_orders(gra, ts_=True):
    """ Get the bond orders of this molecular graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: A dictionary of bond orders, by bond key
    :rtype: dict[frozenset: int or float]
    """
    gra = gra if ts_ else ts_reactants_graph(gra)
    return mdict.by_key_by_position(bonds(gra), bond_keys(gra), BND_ORD_POS)


def atom_implicit_hydrogens(gra):
    """ Get the implicit hydrogen valences of atoms in this molecular graph, as
        a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of implicit hydrogen valences, by atom key
    :rtype: dict[int: int]
    """
    return mdict.by_key_by_position(atoms(gra), atom_keys(gra),
                                    ATM_IMP_HYD_POS)


def atom_stereo_parities(gra):
    """ Get the atom stereo parities of this molecular graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of atom stereo parities, by atom key
    :rtype: dict[int: bool or NoneType]
    """
    return mdict.by_key_by_position(atoms(gra), atom_keys(gra),
                                    ATM_STE_PAR_POS)


def bond_stereo_parities(gra, ts_=True):
    """ Get the bond stereo parities of this molecular graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: A dictionary of bond stereo parities, by bond key
    :rtype: dict[frozenset: bool or NoneType]
    """
    gra = gra if ts_ else ts_reactants_graph(gra)
    return mdict.by_key_by_position(bonds(gra), bond_keys(gra),
                                    BND_STE_PAR_POS)


def stereo_parities(gra, ts_=True):
    """ Get the atom and bond stereo parities of this molecular graph, as a
        single dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: A dictionary of stereo parities, by atom/bond key
    :rtype: dict[int or frozenset: bool or NoneType]
    """
    par_dct = atom_stereo_parities(gra)
    par_dct.update(bond_stereo_parities(gra, ts_=ts_))
    return par_dct


# # TS graph constructor
def ts_graph(gra, frm_bnd_keys, brk_bnd_keys,
             atm_prd_ste_par_dct=None, atm_ts_ste_par_dct=None,
             bnd_prd_ste_par_dct=None, bnd_ts_ste_par_dct=None):
    """ Construct a TS graph from a molecular graph

    :param gra: molecular graph, representing the reactants
    :type gra: automol graph data structure
    :param frm_bnd_keys: Keys to bonds which are forming in the TS
    :type frm_bnd_keys: tuple[frozenset]
    :param brk_bnd_keys: Keys to bonds which are breaking in the TS
    :type brk_bnd_keys: tuple[frozenset]
    :param atm_prd_ste_par_dct: product stereo parities, by atom key
    :type atm_prd_ste_par_dct: dict
    :param atm_ts_ste_par_dct: fleeting TS stereo parities, by atom key
    :type atm_ts_ste_par_dct: dict
    :param bnd_prd_ste_par_dct: product stereo parities, by bond key
    :type bnd_prd_ste_par_dct: dict
    :param bnd_ts_ste_par_dct: fleeting TS stereo parities, by bond key
    :type bnd_ts_ste_par_dct: dict
    :returns: TS graph
    :rtype: automol TS graph data structure
    """
    assert not is_ts_graph(gra), (
        f"Attempting to construct a new TS graph from a TS graph:\n{gra}")

    # Construct a TS graph object from the reactants graph
    bnd_ord_dct = bond_orders(gra)
    tsg = from_data(
        atm_symb_dct=atom_symbols(gra),
        bnd_keys=bond_keys(gra),
        atm_imp_hyd_dct=atom_implicit_hydrogens(gra),
        atm_ste_par_dct=atom_stereo_parities(gra),
        atm_prd_ste_par_dct=atm_prd_ste_par_dct,
        atm_ts_ste_par_dct=atm_ts_ste_par_dct,
        bnd_ord_dct=bnd_ord_dct,
        bnd_ste_par_dct=bond_stereo_parities(gra),
        bnd_prd_ste_par_dct=bnd_prd_ste_par_dct,
        bnd_ts_ste_par_dct=bnd_ts_ste_par_dct,
        ts_=True
    )

    # Encode forming and breaking bonds
    frm_bnd_keys = frozenset(map(frozenset, frm_bnd_keys))
    brk_bnd_keys = frozenset(map(frozenset, brk_bnd_keys))

    frm_ord_dct = {k: 0.1 for k in frm_bnd_keys}
    brk_ord_dct = {k: 0.9 for k in brk_bnd_keys}

    tsg = add_bonds(tsg, frm_bnd_keys, ord_dct=frm_ord_dct, check=False)
    tsg = add_bonds(tsg, brk_bnd_keys, ord_dct=brk_ord_dct, check=False)
    return tsg


# # TS graph getters
def ts_atom_product_stereo_parities(tsg):
    """ Get the product atom stereo parities of this TS graph, as a dictionary

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :returns: A dictionary of atom stereo parities, by atom key
    :rtype: dict[int: bool or NoneType]
    """
    if is_ts_graph(tsg):
        ret = mdict.by_key_by_position(atoms(tsg), atom_keys(tsg),
                                       TS_ATM_PRD_STE_PAR_POS)
    else:
        ret = None
    return ret


def ts_atom_fleeting_stereo_parities(tsg):
    """ Get the fleeting atom stereo parities of this TS graph, as a dictionary

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :returns: A dictionary of atom stereo parities, by atom key
    :rtype: dict[int: bool or NoneType]
    """
    if is_ts_graph(tsg):
        ret = mdict.by_key_by_position(atoms(tsg), atom_keys(tsg),
                                       TS_ATM_FLE_STE_PAR_POS)
    else:
        ret = None
    return ret


def ts_bond_product_stereo_parities(tsg):
    """ Get the product bond stereo parities of this TS graph, as a dictionary

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :returns: A dictionary of bond stereo parities, by bond key
    :rtype: dict[int: bool or NoneType]
    """
    if is_ts_graph(tsg):
        ret = mdict.by_key_by_position(bonds(tsg), bond_keys(tsg),
                                       TS_BND_PRD_STE_PAR_POS)
    else:
        ret = None
    return ret


def ts_bond_fleeting_stereo_parities(tsg):
    """ Get the fleeting bond stereo parities of this TS graph, as a dictionary

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :returns: A dictionary of bond stereo parities, by bond key
    :rtype: dict[int: bool or NoneType]
    """
    if is_ts_graph(tsg):
        ret = mdict.by_key_by_position(bonds(tsg), bond_keys(tsg),
                                       TS_BND_FLE_STE_PAR_POS)
    else:
        ret = None
    return ret


def ts_forming_bond_keys(tsg):
    """ Get the forming bonds from a TS graph

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :returns: The keys to forming bonds
    :rtype: frozenset[frozenset[{int, int}]]
    """
    ord_dct = bond_orders(tsg)
    frm_bnd_keys = [k for k, o in ord_dct.items() if round(o, 1) == 0.1]
    return frozenset(map(frozenset, frm_bnd_keys))


def ts_breaking_bond_keys(tsg):
    """ Get the breaking bonds from a TS graph

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :returns: The keys to breaking bonds
    :rtype: frozenset[frozenset[{int, int}]]
    """
    ord_dct = bond_orders(tsg)
    brk_bnd_keys = [k for k, o in ord_dct.items() if round(o % 1, 1) == 0.9]
    return frozenset(map(frozenset, brk_bnd_keys))


def ts_reacting_bonds(tsg):
    """ Get all of the bonds involved in the reaction

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :returns: The keys to bonds involved in the reaction
    :rtype: frozenset[frozenset[{int, int}]]
    """
    bnd_keys = ts_forming_bond_keys(tsg) | ts_breaking_bond_keys(tsg)
    return bnd_keys


def ts_reacting_atoms(tsg):
    """ Get all of the atoms involved in the reaction

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :returns: The keys to atoms involved in the reaction
    :rtype: frozenset[int]
    """
    bnd_keys = ts_reacting_bonds(tsg)
    atm_keys = frozenset(itertools.chain(*bnd_keys))
    return atm_keys


def ts_without_reacting_bond_orders(tsg, keep_zeros=False):
    """ Remove reacting bonds from a TS graph, replacing them with their values
    for the reactants

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :param keep_zeros: Keep the bonds with a resulting bond order of 0?
    :type keep_zeros: bool
    :returns: The TS graph, without reacting bond orders
    :rtype: automol TS graph data structure
    """
    rxn_bnd_keys = ts_reacting_bonds(tsg)
    rxn_ord_dct = dict_.by_key(bond_orders(tsg), rxn_bnd_keys)
    # Round the bond orders for forming bonds, and remove forming bonds
    tsg = set_bond_orders(tsg, dict_.transform_values(rxn_ord_dct, round))
    if not keep_zeros:
        tsg = without_null_bonds(tsg, except_dummies=True)
    return tsg


def ts_reactants_graph(tsg):
    """ Generate a graph representing the reactants of a TS graph

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :returns: the reactants graph
    :rtype: automol graph data structure
    """
    # Round the bond orders for forming bonds, and remove forming bonds
    tsg = ts_without_reacting_bond_orders(tsg)
    # Remove the extra stereo columns
    gra = from_atoms_and_bonds(atoms(tsg), bonds(tsg), ts_=False)
    return gra


# # setters
def set_atom_symbols(gra, atm_symb_dct):
    """ Set the atom symbols of this molecular graph with a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_symb_dct: A dictionary of atomic symbols, by atom key.
    :type atm_symb_dct: dict[int: str]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_dct = mdict.set_by_key_by_position(atoms(gra), atm_symb_dct,
                                           ATM_SYM_POS)
    return from_atoms_and_bonds(atm_dct, bonds(gra))


def set_bond_orders(gra, bnd_ord_dct):
    """ Set the bond orders of this molecular graph with a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bnd_ord_dct: A dictionary of bond orders, by bond key.
    :type bnd_ord_dct: dict[frozenset: int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    bnd_dct = mdict.set_by_key_by_position(bonds(gra), bnd_ord_dct,
                                           BND_ORD_POS)
    return from_atoms_and_bonds(atoms(gra), bnd_dct)


def set_atom_implicit_hydrogens(gra, atm_imp_hyd_dct):
    """ Set the implicit hydrogen valences of atoms in this molecular graph, as
        a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_imp_hyd_dct: A dictionary of implicit hydrogen valences, by
        atom key
    :type atm_imp_hyd_dct: dict[int: int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_dct = mdict.set_by_key_by_position(atoms(gra), atm_imp_hyd_dct,
                                           ATM_IMP_HYD_POS)
    bnd_dct = bonds(gra)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def set_atom_stereo_parities(gra, atm_par_dct):
    """ Set the atom stereo parities of this molecular graph with a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_par_dct: A dictionary of atom stereo parities, by atom key
    :type atm_par_dct: dict[int: bool or NoneType]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_dct = mdict.set_by_key_by_position(atoms(gra), atm_par_dct,
                                           ATM_STE_PAR_POS)
    return from_atoms_and_bonds(atm_dct, bonds(gra))


def set_bond_stereo_parities(gra, bnd_par_dct):
    """ Set the bond stereo parities of this molecular graph with a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bnd_par_dct: A dictionary of bond stereo parities, by bond key
    :type bnd_par_dct: dict[frozenset: bool or NoneType]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    bnd_dct = mdict.set_by_key_by_position(bonds(gra), bnd_par_dct,
                                           BND_STE_PAR_POS)
    return from_atoms_and_bonds(atoms(gra), bnd_dct)


def set_stereo_parities(gra, par_dct):
    """ Set the atom and bond stereo parities of this molecular graph with a
        single dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param par_dct: A dictionary of stereo parities, by atom/bond key
    :type par_dct: dict[int or frozenset: bool or NoneType]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_par_dct = dict_.filter_by_key(
        par_dct, lambda x: isinstance(x, numbers.Number))
    bnd_par_dct = dict_.filter_by_key(
        par_dct, lambda x: not isinstance(x, numbers.Number))
    atm_dct = mdict.set_by_key_by_position(atoms(gra), atm_par_dct,
                                           ATM_STE_PAR_POS)
    bnd_dct = mdict.set_by_key_by_position(bonds(gra), bnd_par_dct,
                                           BND_STE_PAR_POS)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


# # TS graph getters
def ts_set_atom_product_stereo_parities(tsg, atm_par_dct):
    """ Set the product atom stereo parities of this TS graph with a dictionary

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :param atm_par_dct: A dictionary of atom stereo parities, by atom key
    :type atm_par_dct: dict[int: bool or NoneType]
    :returns: A TS graph
    :rtype: automol graph data structure
    """
    assert is_ts_graph(tsg), (
        f"Attempting to set TS properties on a non-TS graph:\n{string(tsg)}")

    atm_dct = mdict.set_by_key_by_position(atoms(tsg), atm_par_dct,
                                           TS_ATM_PRD_STE_PAR_POS)
    return from_atoms_and_bonds(atm_dct, bonds(tsg))


def ts_set_atom_fleeting_stereo_parities(tsg, atm_par_dct):
    """ Set the fleeting atom stereo parities of this TS graph with a
        dictionary

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :param atm_par_dct: A dictionary of atom stereo parities, by atom key
    :type atm_par_dct: dict[int: bool or NoneType]
    :returns: A TS graph
    :rtype: automol graph data structure
    """
    assert is_ts_graph(tsg), (
        f"Attempting to set TS properties on a non-TS graph:\n{string(tsg)}")

    atm_dct = mdict.set_by_key_by_position(atoms(tsg), atm_par_dct,
                                           TS_ATM_FLE_STE_PAR_POS)
    return from_atoms_and_bonds(atm_dct, bonds(tsg))


def ts_set_bond_product_stereo_parities(tsg, bnd_par_dct):
    """ Set the product bond stereo parities of this TS graph with a dictionary

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :param bnd_par_dct: A dictionary of bond stereo parities, by bond key
    :type bnd_par_dct: dict[frozenset: bool or NoneType]
    :returns: A TS graph
    :rtype: automol graph data structure
    """
    assert is_ts_graph(tsg), (
        f"Attempting to set TS properties on a non-TS graph:\n{string(tsg)}")

    bnd_dct = mdict.set_by_key_by_position(bonds(tsg), bnd_par_dct,
                                           TS_BND_PRD_STE_PAR_POS)
    return from_atoms_and_bonds(atoms(tsg), bnd_dct)


def ts_set_bond_fleeting_stereo_parities(tsg, bnd_par_dct):
    """ Set the fleeting bond stereo parities of this TS graph with a
        dictionary

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :param bnd_par_dct: A dictionary of bond stereo parities, by bond key
    :type bnd_par_dct: dict[frozenset: bool or NoneType]
    :returns: A TS graph
    :rtype: automol graph data structure
    """
    assert is_ts_graph(tsg), (
        f"Attempting to set TS properties on a non-TS graph:\n{string(tsg)}")

    bnd_dct = mdict.set_by_key_by_position(bonds(tsg), bnd_par_dct,
                                           TS_BND_FLE_STE_PAR_POS)
    return from_atoms_and_bonds(atoms(tsg), bnd_dct)


# # I/O
def string(gra, one_indexed=True):
    """ Generate a string representation of the graph, in YAML format

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param one_indexed: Switch to one-indexing for keys?
    :type one_indexed: bool
    :returns: A string representation of the molecular graph
    :rtype: str
    """
    yaml_gra_dct = yaml_dictionary(gra, one_indexed=one_indexed)
    gra_str = yaml.dump(yaml_gra_dct, default_flow_style=None, sort_keys=False)
    return gra_str


def yaml_dictionary(gra, one_indexed=True):
    """ Generate a YAML dictionary representation of the graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param one_indexed: Switch to one-indexing for keys?
    :type one_indexed: bool
    :returns: A YAML-friendly dictionary representation of the graph
    :rtype: dict
    """
    if one_indexed:
        # shift to one-indexing when we print
        atm_key_dct = {atm_key: atm_key+1 for atm_key in atom_keys(gra)}
        gra = relabel(gra, atm_key_dct)

    yaml_atm_dct = atoms(gra)
    yaml_bnd_dct = bonds(gra)

    ts_ = is_ts_graph(gra)
    atm_nprops = 5 if ts_ else 3
    bnd_nprops = 4 if ts_ else 2

    # prepare the atom dictionary
    yaml_atm_dct = dict(sorted(yaml_atm_dct.items()))
    yaml_atm_dct = dict_.transform_values(
        yaml_atm_dct, lambda x: dict(zip(ATM_PROP_NAMES[:atm_nprops], x)))

    # perpare the bond dictionary
    yaml_bnd_dct = dict_.transform_keys(
        yaml_bnd_dct, lambda x: tuple(sorted(x)))
    yaml_bnd_dct = dict(sorted(yaml_bnd_dct.items()))
    yaml_bnd_dct = dict_.transform_keys(
        yaml_bnd_dct, lambda x: '-'.join(map(str, x)))
    yaml_bnd_dct = dict_.transform_values(
        yaml_bnd_dct, lambda x: dict(zip(BND_PROP_NAMES[:bnd_nprops], x)))

    yaml_gra_dct = {'atoms': yaml_atm_dct, 'bonds': yaml_bnd_dct}
    return yaml_gra_dct


def from_string(gra_str, one_indexed=True):
    """ Generate a graph from a string representation, in YAML format

    :param gra_str: A string representation of the molecular graph
    :type gra_str: str
    :param one_indexed: Assume one-indexing for string keys?
    :type one_indexed: bool
    :returns: molecular graph
    :rtype: automol graph data structure
    """
    yaml_gra_dct = yaml.load(gra_str, Loader=yaml.FullLoader)
    gra = from_yaml_dictionary(yaml_gra_dct, one_indexed=one_indexed)
    return gra


def from_yaml_dictionary(yaml_gra_dct, one_indexed=True):
    """ Generate a graph from a YAML dictionary representation

    :param yaml_gra_dct: A YAML-friendly dictionary representation of the graph
    :type yaml_gra_dct: dict
    :param one_indexed: Assume one-indexing for YAML dict keys?
    :type one_indexed: bool
    :returns: molecular graph
    :rtype: automol graph data structure
    """
    atm_dct = yaml_gra_dct['atoms']
    bnd_dct = yaml_gra_dct['bonds']

    atm_dct = dict_.transform_values(
        atm_dct, lambda x: tuple(map(x.__getitem__, ATM_PROP_NAMES[:len(x)])))

    bnd_dct = dict_.transform_keys(
        bnd_dct, lambda x: frozenset(map(int, x.split('-'))))

    bnd_dct = dict_.transform_values(
        bnd_dct, lambda x: tuple(map(x.__getitem__, BND_PROP_NAMES[:len(x)])))

    gra = from_atoms_and_bonds(atm_dct, bnd_dct)

    if one_indexed:
        # revert one-indexing if the input is one-indexed
        atm_key_dct = {atm_key: atm_key-1 for atm_key in atom_keys(gra)}
        gra = relabel(gra, atm_key_dct)

    return gra


# # conversions
def frozen(gra):
    """ Generate a hashable, sortable, immutable representation of the graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A hashable, sortable, immutable representation of the graph
    :rtype: (tuple, tuple)
    """
    atm_keys = sorted(atom_keys(gra))
    bnd_keys = sorted(bond_keys(gra), key=sorted)

    # make it sortable by replacing Nones with -infinity
    atm_vals = numpy.array(dict_.values_by_key(atoms(gra), atm_keys),
                           dtype=object)
    bnd_vals = numpy.array(dict_.values_by_key(bonds(gra), bnd_keys),
                           dtype=object)
    atm_vals[numpy.equal(atm_vals, None)] = -numpy.inf
    bnd_vals[numpy.equal(bnd_vals, None)] = -numpy.inf

    frz_atms = tuple(zip(atm_keys, map(tuple, atm_vals)))
    frz_bnds = tuple(zip(bnd_keys, map(tuple, bnd_vals)))
    return (frz_atms, frz_bnds)


def formula(gra):
    """ Generate a stoichiometric formula dictionary from a molecular graph.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :rtype: dict[str: int]
    """

    gra = explicit(gra)
    syms = atom_symbols(gra).values()
    fml = util.formula_from_symbols(syms)

    return fml


# # properties
def atom_count(gra, symb=None, heavy_only=False, dummy=False,
               with_implicit=True, keys=None):
    """ Count the number of atoms in the molecule

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symb: Atomic symbol to count, defaults to None
    :type symb: str or NoneType
    :param heavy_only: Restrict the count to heavy atoms?
    :type heavy_only: bool
    :param dummy: Include dummy atoms?, defaults to False
    :type dummy: bool, optional
    :param with_implicit: Include implicit hydrogens?, defaults to True
    :type with_implicit: bool, optional
    :param keys: optionally, restrict the count to a subset of atoms
    :type keys: tuple[int]
    :return: The number of atoms
    :rtype: int
    """
    if not dummy and symb != 'X':
        gra = without_dummy_atoms(gra)

    # Get the keys
    keys = atom_keys(gra) if keys is None else keys
    symb_dct = atom_symbols(gra)
    if heavy_only:
        # If only including heavy atoms, filter non-heavy atoms out
        symb_dct = dict_.filter_by_value(
            symb_dct, lambda s: ptab.to_number(s) != 1)
    # Restrict count to the requested subset of atoms
    symbs = [symb_dct[k] for k in keys if k in symb_dct]
    natms = len(symbs) if symb is None else symbs.count(symb)

    if with_implicit and symb in ('H', None) and not heavy_only:
        atm_imp_hyd_dct = atom_implicit_hydrogens(gra)
        natms += sum(atm_imp_hyd_dct.values())

    return natms


def electron_count(gra, charge=0):
    """ Count the number of electrons in the molecule, based on its charge

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param charge: The molecular charge, defaults to 0
    :type charge: int, optional
    :return: The number of electrons
    :rtype: int
    """
    atm_symb_dct = atom_symbols(explicit(gra))
    nelec = sum(map(ptab.to_number, atm_symb_dct.values())) - charge
    return nelec


def atom_stereo_keys(gra, symb=None, excl_symbs=(), ts_all=False):
    """ Get the keys of atom stereo-centers this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symb: Optionally, restrict this to atoms with a particular
        atomic symbol (e.g., 'H' for hydrogens).
    :type symb: str
    :param excl_symbs: Optionally, exclude atoms with particular atomic
        symbols.
    :type excl_symbs: tuple[str]
    :param ts_all: If this is a TS graph, return *all* stereo keys?
        Otherwise, only reactant stereo keys will be returned.
    :type ts_all: bool
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    keys = atom_keys(gra, symb=symb, excl_symbs=excl_symbs)
    par_dct = dict_.by_key(atom_stereo_parities(gra), keys)
    ste_keys = dict_.keys_by_value(par_dct, lambda x: x is not None)
    if ts_all and is_ts_graph(gra):
        # Add product stereo atoms
        prd_par_dct = dict_.by_key(ts_atom_product_stereo_parities(gra), keys)
        ste_keys |= dict_.keys_by_value(prd_par_dct, lambda x: x is not None)

        # Add fleeting stereo atoms
        ts_par_dct = dict_.by_key(ts_atom_fleeting_stereo_parities(gra), keys)
        ste_keys |= dict_.keys_by_value(ts_par_dct, lambda x: x is not None)
    return ste_keys


def bond_stereo_keys(gra, ts_all=False):
    """ Get the keys of bond stereo-centers this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_all: If this is a TS graph, return *all* stereo keys?
        Otherwise, only reactant stereo keys will be returned.
    :type ts_all: bool
    :returns: The bond keys
    :rtype: frozenset[{int, int}]
    """
    ste_keys = dict_.keys_by_value(bond_stereo_parities(gra),
                                   lambda x: x is not None)
    if ts_all and is_ts_graph(gra):
        ste_keys |= dict_.keys_by_value(ts_bond_product_stereo_parities(gra),
                                        lambda x: x is not None)
        ste_keys |= dict_.keys_by_value(ts_bond_fleeting_stereo_parities(gra),
                                        lambda x: x is not None)
    return ste_keys


def stereo_keys(gra, ts_all=False):
    """ Get the keys of atom and bond stereo-centers this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_all: If this is a TS graph, return *all* stereo keys?
        Otherwise, only reactant stereo keys will be returned.
    :type ts_all: bool
    :returns: The atom and bond keys
    :rtype: frozenset
    """
    ste_keys = (atom_stereo_keys(gra, ts_all=ts_all) |
                bond_stereo_keys(gra, ts_all=ts_all))
    return ste_keys


def has_atom_stereo(gra, symb=None, excl_symbs=(), ts_all=False):
    """ Does this graph have atom stereochemistry?

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symb: Optionally, restrict this to atoms with a particular
        atomic symbol (e.g., 'H' for hydrogens).
    :type symb: str
    :param excl_symbs: Optionally, exclude atoms with particular atomic
        symbols.
    :type excl_symbs: tuple[str]
    :param ts_all: If this is a TS graph, return *all* stereo keys?
        Otherwise, only reactant stereo keys will be returned.
    :type ts_all: bool
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return bool(
        atom_stereo_keys(gra, symb=symb, excl_symbs=excl_symbs, ts_all=ts_all))


def has_bond_stereo(gra, ts_all=False):
    """ Does this graph have bond stereochemistry?

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_all: If this is a TS graph, return *all* stereo keys?
        Otherwise, only reactant stereo keys will be returned.
    :type ts_all: bool
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return bool(bond_stereo_keys(gra, ts_all=ts_all))


def has_stereo(gra, ts_all=False):
    """ Does this graph have stereochemistry of any kind?

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_all: If this is a TS graph, return *all* stereo keys?
        Otherwise, only reactant stereo keys will be returned.
    :type ts_all: bool
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return bool(stereo_keys(gra, ts_all=ts_all))


def has_pi_bonds(gra):
    """ Does this graph have pi bonds?

    Returns `True` if it has bond orders other than 0, 1, or, for TS graphs,
    0.1 or 0.9

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    bnd_ords = bond_orders(gra).values()
    return any(round(o, 1) not in (0, 0.1, 0.9, 1) for o in bnd_ords)


def is_ts_graph(gra):
    """ Is this a TS graph?

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: `True` if it is, `False` if it isn't
    :rtype: bool
    """
    atm_nprops = dict_.transform_values(atoms(gra), len).values()
    bnd_nprops = dict_.transform_values(bonds(gra), len).values()
    if all(n == 5 for n in atm_nprops) and all(n == 4 for n in bnd_nprops):
        ret = True
    elif all(n == 3 for n in atm_nprops) and all(n == 2 for n in bnd_nprops):
        ret = False
    else:
        raise ValueError(f"Not a valid graph data structure:\n{gra}\n")
    return ret


def atomic_numbers(gra):
    """ Get atomic numbers for atoms in this graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of atomic numbers, by atom key
    :rtype: dict[int: int]
    """
    symb_dct = atom_symbols(gra)
    anum_dct = dict_.transform_values(symb_dct, ptab.to_number)
    return anum_dct


def mass_numbers(gra):
    """ Get mass numbers for atoms in this graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of mass numbers, by atom key
    :rtype: dict[int: float]
    """
    symb_dct = atom_symbols(gra)
    mnum_dct = dict_.transform_values(symb_dct, ptab.to_mass_number)
    return mnum_dct


def atomic_valences(gra):
    """ Get atomic valences for atoms in this graph, as a dictionary

    Valences here refers to the number of bonds that each element is
    intrinsically capable of forming. Ex: 'C' => 4, 'O' => 2, etc.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of atomic valences, by atom key
    :rtype: dict[int: int]
    """
    atm_symb_dct = atom_symbols(gra)
    atm_vlc_dct = dict_.transform_values(atm_symb_dct, ptab.valence)
    return atm_vlc_dct


def atom_lone_pairs(gra):
    """ Get the number of lone pairs for atoms in this graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of lone pair counts, by atom key
    :rtype: dict[int: int]
    """
    atm_symb_dct = atom_symbols(gra)
    atm_lpc_dct = dict_.transform_values(atm_symb_dct, ptab.lone_pair_count)
    return atm_lpc_dct


def lone_pair_atom_keys(gra):
    """ Get the keys of atoms in this graph that have lone pairs

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    lpc_dct = atom_lone_pairs(gra)
    return frozenset(key for key, lpc in lpc_dct.items() if lpc > 0)


def atom_van_der_waals_radius(gra, key, angstrom=True):
    """ Get the van der Waals radius of an atom in the graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param key: the atom key
    :type key: int
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    symb_dct = atom_symbols(gra)
    symb = symb_dct[key]
    rad = ptab.van_der_waals_radius(symb)
    if angstrom:
        rad *= phycon.BOHR2ANG
    return rad


def atom_bond_counts(gra, bond_order=True, with_implicit=True):
    """ Get the bond counts for the atoms in this graph, as a dictionary

    The bond count will be based on the graph's bond orders, so one must
    convert to a Kekule graph in order to include pi bonds in the count.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bond_order: Use the bond orders in the graph?, defaults to True
    :type bond_order: bool, optional
    :param with_implicit: Include implicit hydrogens?, defaults to True
    :type with_implicit: bool, optional
    :returns: A dictionary of bond counts, by atom key
    :rtype: dict[int: int]
    """
    atm_keys = list(atom_keys(gra))
    if with_implicit:
        gra = explicit(gra)
    gra = ts_reactants_graph(gra)
    if not bond_order:
        gra = without_pi_bonds(gra)

    atm_nbhs = dict_.values_by_key(atom_neighborhoods(gra), atm_keys)
    atm_nbnds = [sum(bond_orders(nbh).values()) for nbh in atm_nbhs]
    atm_nbnd_dct = dict_.transform_values(dict(zip(atm_keys, atm_nbnds)), int)
    return atm_nbnd_dct


def atom_unpaired_electrons(gra, bond_order=True):
    """ The number of unpaired electrons on each atom, calculated as the atomic
    valences minus the bond count

    In a non-Kekule graph, these will be *either* radical *or* pi-bonding
    electrons. In a Kekule graph, they will be only radical electrons.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bond_order: Use the bond orders in the graph?, defaults to True
    :type bond_order: bool, optional
    :returns: the atom unsaturated valences, by atom key
    :rtype: dict
    """
    atm_keys = list(atom_keys(gra))
    gra = ts_reactants_graph(gra)
    if not bond_order:
        gra = without_pi_bonds(gra)

    atm_bnd_vlcs = dict_.values_by_key(atom_bond_counts(gra), atm_keys)
    atm_tot_vlcs = dict_.values_by_key(atomic_valences(gra), atm_keys)
    atm_rad_vlcs = numpy.subtract(atm_tot_vlcs, atm_bnd_vlcs)
    atm_unp_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_rad_vlcs)), int)
    return atm_unp_dct


def bond_unpaired_electrons(gra, bond_order=True):
    """ The number of adjacent pairs of unpaired electrons across each bond,
    which are available for pi-bonding

    Calculated as the minimum number of unpaired electrons for the atoms on
    either side of the bond.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bond_order: Use the bond orders in the graph?, defaults to True
    :type bond_order: bool, optional
    :returns: the bond unsaturated valences, by bond key
    :rtype: dict
    """
    bnd_keys = list(bond_keys(gra))
    gra = ts_reactants_graph(gra)
    if not bond_order:
        gra = without_pi_bonds(gra)

    # determine unsaturated valences for each atom
    atm_unp_dct = atom_unpaired_electrons(gra)
    # a bond's unsaturated is set by its limiting atom
    bnd_unsats = [min(map(atm_unp_dct.__getitem__, k)) for k in bnd_keys]
    bnd_unsat_dct = dict(zip(bnd_keys, bnd_unsats))
    return bnd_unsat_dct


def tetrahedral_atom_keys(gra):
    """ Keys to tetrahedral atoms (possible stereo centers).

    Atoms will be considered tetrahedral if either:
        a. They are bonded to 4 other atoms.
        b. They are bonded to 3 other atoms and have one lone pair.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    gra = ts_reactants_graph(gra)
    bnd_vlc_dct = atom_bond_counts(gra)
    lpc_dct = atom_lone_pairs(gra)

    def _is_tetrahedral(key):
        bnd_vlc = bnd_vlc_dct[key]
        lpc = lpc_dct[key]

        ret = (lpc == 0) and (bnd_vlc == 4) or (lpc == 1) and (bnd_vlc == 3)
        return ret

    tet_atm_keys = frozenset(filter(_is_tetrahedral, atom_keys(gra)))
    return tet_atm_keys


def maximum_spin_multiplicity(gra, bond_order=True):
    """ The highest possible spin multiplicity for this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bond_order: Use the bond orders in the graph?, defaults to True
    :type bond_order: bool, optional
    :returns: The maximum spin multiplicity, as a number
    :rtype: int
    """
    atm_rad_dct = atom_unpaired_electrons(gra, bond_order=bond_order)
    return sum(atm_rad_dct.values()) + 1


def possible_spin_multiplicities(gra, bond_order=True):
    """ Possible spin multiplicities for this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bond_order: Use the bond orders in the graph?, defaults to True
    :type bond_order: bool, optional
    :returns: The list of possible spin multiplicities
    :rtype: tuple
    """
    mult_max = maximum_spin_multiplicity(gra, bond_order=bond_order)
    mult_min = 2 if mult_max % 2 == 0 else 1
    mults = tuple(range(mult_min, mult_max+1, 2))
    return mults


def atom_symbol_keys(gra):
    """ Group the atom keys by atomic symbol

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary whose keys are the atomic symbols and whose values
        are the atom keys associated with each
    :rtype: dict[str: list]
    """

    idx_symb_dct = atom_symbols(gra)

    symb_idx_dct = {}
    for idx, symb in idx_symb_dct.items():
        if symb not in symb_idx_dct:
            symb_idx_dct[symb] = [idx]
        else:
            symb_idx_dct[symb].append(idx)

    return symb_idx_dct


def backbone_hydrogen_keys(gra):
    """ Get the backbone hydrogen keys of this molecular graph

    These are hydrogen keys which cannot be made implicit, because they are
    part of the backbone. There are two cases: (1.) one of the hydrogens in
    H2 must be considered a backbone hydrogen, and (2.) any multivalent
    hydrogen must be treated as part of the backbone.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    hyd_keys = atom_keys(gra, symb='H')
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    def _is_backbone(hyd_key):
        is_h2 = all(ngb_key in hyd_keys and hyd_key < ngb_key
                    for ngb_key in atm_ngb_keys_dct[hyd_key])
        is_multivalent = len(atm_ngb_keys_dct[hyd_key]) > 1
        return is_h2 or is_multivalent

    bbn_hyd_keys = frozenset(filter(_is_backbone, hyd_keys))
    return bbn_hyd_keys


def nonbackbone_hydrogen_keys(gra):
    """ Get the hydrogen keys of this graph, with or without backbone hydrogens

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    return atom_keys(gra, symb='H') - backbone_hydrogen_keys(gra)


def atom_nonbackbone_hydrogen_keys(gra):
    """ Get the hydrogen keys of each atom, with or without backbone hydrogens,
    as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary giving the hydrogen keys for each atom, by key.
    :rtypoe: dict[int: frozenset]
    """
    hyd_keys = nonbackbone_hydrogen_keys(gra)
    atm_hyd_keys_dct = dict_.transform_values(
        atoms_neighbor_atom_keys(gra), lambda x: x & hyd_keys)
    return atm_hyd_keys_dct


def backbone_keys(gra):
    """ Get the backbone atom keys of this graph, including backbone hydrogens

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    return atom_keys(gra) - nonbackbone_hydrogen_keys(gra)


def backbone_bond_keys(gra, terminal=False):
    """ Get the backbone bond keys of this graph (bonds between backbone atoms)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param terminal: Include bond keys with terminal atoms?
    :type terminal: bool, optional
    :return: The bond keys
    :rtype: frozenset(frozenset[int])
    """
    bbn_keys = backbone_keys(gra)

    # If requested, remove terminal atom keys
    if not terminal:
        bbn_keys -= terminal_atom_keys(gra)

    bnd_keys = bond_keys(gra)
    return frozenset(filter(lambda x: x == (x & bbn_keys), bnd_keys))


def atom_backbone_hydrogen_keys(gra):
    """ Get backbone hydrogen keys for each atom, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary giving the hydrogen keys for each atom, by key
    :rtype: dict[int: frozenset]
    """
    hyd_keys = backbone_hydrogen_keys(gra)
    atm_hyd_keys_dct = dict_.transform_values(
        atoms_neighbor_atom_keys(gra), lambda x: x & hyd_keys)
    return atm_hyd_keys_dct


def terminal_atom_keys(gra, backbone=True):
    """ Get the backbone atom keys of this graph, including backbone hydrogens

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param backbone: Restrict this to backbone atoms?
    :type backbone: bool
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    if backbone:
        gra = implicit(gra)

    atm_nkeys_dct = atoms_neighbor_atom_keys(gra)
    atm_keys = [key for key, nkeys in atm_nkeys_dct.items() if len(nkeys) <= 1]
    return frozenset(atm_keys)


def unsaturated_atom_keys(gra):
    """ Get the keys of unsaturated (radical or pi-bonded) atoms

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    atm_unp_dct = atom_unpaired_electrons(gra, bond_order=False)
    unsat_atm_keys = frozenset(dict_.keys_by_value(atm_unp_dct, bool))
    return unsat_atm_keys


def unsaturated_bond_keys(gra):
    """ Get the keys of unsaturated (radical or pi-bonded) bonds

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The bond keys
    :rtype: frozenset[frozenset[int]]
    """
    bnd_unp_dct = bond_unpaired_electrons(gra, bond_order=False)
    unsat_bnd_keys = frozenset(dict_.keys_by_value(bnd_unp_dct, bool))
    return unsat_bnd_keys


def angle_keys(gra):
    """ Get triples of keys for pairs of adjacent bonds, with the central atom
    in the middle

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The angle keys
    :rtype: frozenset[tuple[int]]
    """
    bnd_keys = bond_keys(gra)

    ang_keys = []
    for bnd_key1, bnd_key2 in itertools.combinations(bnd_keys, r=2):
        if bnd_key1 != bnd_key2 and bnd_key1 & bnd_key2:
            atm2_key, = bnd_key1 & bnd_key2
            atm1_key, = bnd_key1 - {atm2_key}
            atm3_key, = bnd_key2 - {atm2_key}
            ang_keys.append((atm1_key, atm2_key, atm3_key))
            ang_keys.append((atm3_key, atm2_key, atm1_key))

    return frozenset(ang_keys)


# # relabeling and changing keys
def relabel(gra, atm_key_dct):
    """ Relabel the atoms in the graph with new keys, using a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_key_dct: New keys for a subset of the atoms, by current atom key
    :type atm_key_dct: dict[int: int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    orig_atm_keys = atom_keys(gra)
    assert set(atm_key_dct.keys()) <= orig_atm_keys, (
        f'{set(atm_key_dct.keys())}\n{orig_atm_keys}')

    new_atm_key_dct = dict(zip(orig_atm_keys, orig_atm_keys))
    new_atm_key_dct.update(atm_key_dct)

    _relabel_atom_key = new_atm_key_dct.__getitem__

    def _relabel_bond_key(bnd_key):
        return frozenset(map(_relabel_atom_key, bnd_key))

    atm_dct = dict_.transform_keys(atoms(gra), _relabel_atom_key)
    bnd_dct = dict_.transform_keys(bonds(gra), _relabel_bond_key)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def standard_keys(gra):
    """ Relabel the atoms in the graph with standard zero-indexed keys

    The new keys will follow the same sort order as the original ones.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_key_dct = dict(map(reversed, enumerate(sorted(atom_keys(gra)))))
    return relabel(gra, atm_key_dct)


def standard_keys_for_sequence(gras):
    """ Give standard, non-overlapping keys to a sequence of graphs and return
    the relabelling dictionary along with each graph

    The new keys will start from zero in the first graph and count up
    sequentially, with each next graph starting from the next possible integer
    value

    :param gras: A sequence of molecular graphs
    :returns: A sequence of relabelled molecular graphs, and a sequence of
        relabelling dictionaries for each
    """
    atm_key_dcts = []

    shift = 0
    for gra in gras:
        natms = atom_count(gra, dummy=True, with_implicit=False)

        atm_key_dct = {atm_key: idx+shift
                       for idx, atm_key in enumerate(sorted(atom_keys(gra)))}
        atm_key_dcts.append(atm_key_dct)

        shift += natms

    gras = tuple(relabel(gra, atm_key_dct)
                 for gra, atm_key_dct in zip(gras, atm_key_dcts))
    atm_key_dcts = tuple(atm_key_dcts)

    return gras, atm_key_dcts


def relabel_for_zmatrix(gra, zma_keys, dummy_key_dct):
    """ Relabel a graph to line up keys with a geometry => z-matrix conversion

    The original graph keys should correspond to the geometry that was used to
    generate the z-matrix. This function will insert the dummy atoms and
    sort/relabel the graph to match the z-matrix indices.

    Note: This assumes that the conversion was performed using
    `auotmol.geom.zmatrix`, which returns `zma_keys` and `dummy_key_dct`.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param zma_keys: graph keys in the order they appear in the z-matrix
    :type zma_keys: list[int]
    :param dummy_key_dct: dummy keys introduced on z-matrix conversion, by atom
        they are attached to
    :type dummy_key_dct: dict[int: int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    gra = add_dummy_atoms(gra, dummy_key_dct)
    key_dct = dict(map(reversed, enumerate(zma_keys)))
    gra = relabel(gra, key_dct)
    return gra


def relabel_for_geometry(gra):
    """ Relabel a graph to line up keys with a z-matrix => geometry conversion

    The original graph keys should correspond to the z-matrix that was used to
    generate the geometry. This function will remove dummy atoms and relabel
    the keys to match.

    Note: This assumes that the conversion was performed using
    `auotmol.zmat.geometry`.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    def _shift_remove_dummy_atom(gra, dummy_key):
        keys = sorted(atom_keys(gra))
        idx = keys.index(dummy_key)
        key_dct = {}
        key_dct.update({k: k for k in keys[:idx]})
        key_dct.update({k: k-1 for k in keys[(idx+1):]})
        gra = remove_atoms(gra, [dummy_key], stereo=True)
        gra = relabel(gra, key_dct)
        return gra

    dummy_keys = sorted(atom_keys(gra, symb='X'))

    for dummy_key in reversed(dummy_keys):
        gra = _shift_remove_dummy_atom(gra, dummy_key)
    return gra


def negate_nonbackbone_hydrogen_keys(gra):
    """ Flip the signs of hydrogen keys

    :param gra: molecular graph
    :type gra: automol graph data structure
    :return: molecular graph with hydrogen keys negated
    :rtype: automol graph data structure
    """
    gra = ts_reactants_graph(gra)
    hyd_keys = nonbackbone_hydrogen_keys(gra)
    atm_key_dct = {k: -abs(k) for k in hyd_keys}
    return relabel(gra, atm_key_dct)


# # add/remove/insert/without
def add_atoms(gra, symb_dct, imp_hyd_dct=None, ste_par_dct=None,
              prd_ste_par_dct=None, ts_ste_par_dct=None, check=True):
    """ Add atoms to this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symb_dct: atomic symbols, by atom key
    :type symb_dct: dict
    :param imp_hyd_dct: the number of implicit hydrogens associated with
        each atom, by atom key
    :type imp_hyd_dct: dict
    :param ste_par_dct: stereo parities, by atom key; (in TS graphs, reactant
        stereo parities)
    :type ste_par_dct: dict
    :param prd_ste_par_dct: product stereo parities, by atom key; (TS graphs
        only)
    :type prd_ste_par_dct: dict
    :param ts_ste_par_dct: fleeting TS stereo parities, by atom key; (TS graphs
        only)
    :type ts_ste_par_dct: dict
    :param check: Check that we aren't trying to add atoms with duplicate keys?
    :type check: bool
    :return: molecular graph (TS or non-TS)
    :rtype: automol graph data structure
    """
    ts_ = is_ts_graph(gra)
    atm_keys = atom_keys(gra)
    atm_symb_dct = atom_symbols(gra)
    atm_imp_hyd_dct = atom_implicit_hydrogens(gra)
    atm_ste_par_dct = atom_stereo_parities(gra)
    atm_prd_ste_par_dct = ts_atom_product_stereo_parities(gra)
    atm_ts_ste_par_dct = ts_atom_fleeting_stereo_parities(gra)

    keys = set(symb_dct.keys())
    imp_hyd_dct = {} if imp_hyd_dct is None else imp_hyd_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct
    if ts_:
        prd_ste_par_dct = {} if prd_ste_par_dct is None else prd_ste_par_dct
        ts_ste_par_dct = {} if ts_ste_par_dct is None else ts_ste_par_dct

    if check:
        assert not keys & atm_keys, (
            f'{keys} and {atm_keys} have a non-empty intersection')

    assert not keys & atm_keys
    assert set(imp_hyd_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys
    if ts_:
        assert set(prd_ste_par_dct.keys()) <= keys
        assert set(ts_ste_par_dct.keys()) <= keys

    atm_symb_dct.update(symb_dct)
    atm_imp_hyd_dct.update(imp_hyd_dct)
    atm_ste_par_dct.update(ste_par_dct)
    if ts_:
        atm_prd_ste_par_dct.update(prd_ste_par_dct)
        atm_ts_ste_par_dct.update(ts_ste_par_dct)

    atm_dct = atoms_from_data(
        atm_symb_dct=atm_symb_dct, atm_imp_hyd_dct=atm_imp_hyd_dct,
        atm_ste_par_dct=atm_ste_par_dct,
        atm_prd_ste_par_dct=atm_prd_ste_par_dct,
        atm_ts_ste_par_dct=atm_ts_ste_par_dct, ts_=ts_)
    bnd_dct = bonds(gra)
    gra = from_atoms_and_bonds(atm_dct=atm_dct, bnd_dct=bnd_dct)
    return gra


def add_bonds(gra, keys, ord_dct=None, ste_par_dct=None, prd_ste_par_dct=None,
              ts_ste_par_dct=None, check=True):
    """ Add bonds to this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param keys: bond keys
    :type keys: set
    :param ord_dct: bond orders, by bond key
    :type ord_dct: dict
    :param ste_par_dct: stereo parities, by bond key; (in TS graphs, reactant
        stereo parities)
    :type ste_par_dct: dict
    :param prd_ste_par_dct: product stereo parities, by bond key; (TS graphs
        only)
    :type prd_ste_par_dct: dict
    :param ts_ste_par_dct: fleeting TS stereo parities, by bond key; (TS graphs
        only)
    :type ts_ste_par_dct: dict
    :param check: Check that we aren't trying to add bonds with duplicate keys?
    :type check: bool
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    ts_ = is_ts_graph(gra)
    bnd_keys = set(bond_keys(gra))
    bnd_ord_dct = bond_orders(gra)
    bnd_ste_par_dct = bond_stereo_parities(gra)
    bnd_prd_ste_par_dct = ts_bond_product_stereo_parities(gra)
    bnd_ts_ste_par_dct = ts_bond_fleeting_stereo_parities(gra)

    keys = set(map(frozenset, keys))
    ord_dct = {} if ord_dct is None else ord_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct
    if ts_:
        prd_ste_par_dct = {} if prd_ste_par_dct is None else prd_ste_par_dct
        ts_ste_par_dct = {} if ts_ste_par_dct is None else ts_ste_par_dct

    ord_dct = dict_.transform_keys(ord_dct, frozenset)
    ste_par_dct = dict_.transform_keys(ste_par_dct, frozenset)
    if ts_:
        prd_ste_par_dct = dict_.transform_keys(prd_ste_par_dct, frozenset)
        ts_ste_par_dct = dict_.transform_keys(ts_ste_par_dct, frozenset)

    if check:
        assert not keys & bnd_keys, (
            f'{keys} and {bnd_keys} have a non-empty intersection')

    assert set(ord_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys
    if ts_:
        assert set(prd_ste_par_dct.keys()) <= keys
        assert set(ts_ste_par_dct.keys()) <= keys

    bnd_keys.update(keys)
    bnd_ord_dct.update(ord_dct)
    bnd_ste_par_dct.update(ste_par_dct)
    if ts_:
        bnd_prd_ste_par_dct.update(prd_ste_par_dct)
        bnd_ts_ste_par_dct.update(ts_ste_par_dct)

    atm_dct = atoms(gra)
    bnd_dct = bonds_from_data(
        bnd_keys=bnd_keys, bnd_ord_dct=bnd_ord_dct,
        bnd_ste_par_dct=bnd_ste_par_dct,
        bnd_prd_ste_par_dct=bnd_prd_ste_par_dct,
        bnd_ts_ste_par_dct=bnd_ts_ste_par_dct, ts_=ts_)

    gra = from_atoms_and_bonds(atm_dct=atm_dct, bnd_dct=bnd_dct)
    return gra


def remove_atoms(gra, atm_keys, check=True, stereo=True):
    """ Remove atoms from this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_keys: The keys of atoms to be removed
    :type atm_keys: list[int]
    :param check: Check that these atoms actually exist in the graph?
    :type check: bool
    :param stereo: Keep stereo information?
    :type stereo: bool
    :return: molecular graph (TS or non-TS)
    :rtype: automol graph data structure
    """
    all_atm_keys = atom_keys(gra)
    atm_keys = set(atm_keys)

    if check:
        assert atm_keys <= all_atm_keys

    atm_keys_left = all_atm_keys - atm_keys
    return subgraph(gra, atm_keys_left, stereo=stereo)


def remove_bonds(gra, bnd_keys, check=True, stereo=True):
    """ Remove bonds from this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bnd_keys: The keys of bonds to be removed
    :type bnd_keys: list[int]
    :param check: Check that these bonds actually exist in the graph?
    :type check: bool
    :param stereo: Keep stereo information?
    :type stereo: bool
    :return: molecular graph (TS or non-TS)
    :rtype: automol graph data structure
    """
    if not stereo:
        gra = without_stereo(gra)

    all_bnd_keys = bond_keys(gra)
    bnd_keys = set(map(frozenset, bnd_keys))

    if check:
        assert bnd_keys <= all_bnd_keys

    bnd_keys = all_bnd_keys - bnd_keys
    atm_dct = atoms(gra)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def change_implicit_hydrogens(gra, imp_hyd_change_dct):
    """ Change the implicit hydrogen count for atoms in this graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param imp_hyd_change_dct: A dictionary telling how many implicit
        hydrogens to add (positive integer) or remove (negative integer)
        for each atom.
    :type imp_hyd_change_dct: dict[int: int]
    :return: molecular graph
    :rtype: automol graph data structure
    """
    atm_keys = list(imp_hyd_change_dct.keys())
    atm_imp_hyds = numpy.add(
        dict_.values_by_key(atom_implicit_hydrogens(gra), atm_keys),
        dict_.values_by_key(imp_hyd_change_dct, atm_keys))
    assert all(atm_imp_hyd >= 0 for atm_imp_hyd in atm_imp_hyds)
    atm_imp_hyd_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_imp_hyds)), int)
    return set_atom_implicit_hydrogens(gra, atm_imp_hyd_dct)


def add_atom_explicit_hydrogens(gra, exp_hyd_keys_dct):
    """ Add explicit hydrogens by atom

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param exp_hyd_keys_dct: A dictionary of new keys for explicit hydrogens to
        be added. The keys of this dictionary are the parent atoms (already in
        `gra`) that these new hydrogens will be connected to.
    :type exp_hyd_keys_dct: dict[int: frozenset]
    :return: molecular graph
    :rtype: automol graph data structure
    """
    assert set(exp_hyd_keys_dct.keys()) <= atom_keys(gra), (
        f'{set(exp_hyd_keys_dct.keys())}'
        ' !<= '
        f'{atom_keys(gra)}'
    )
    for atm_key, atm_exp_hyd_keys in exp_hyd_keys_dct.items():
        assert not set(atm_exp_hyd_keys) & atom_keys(gra)
        atm_exp_hyd_bnd_keys = {frozenset({atm_key, atm_exp_hyd_key})
                                for atm_exp_hyd_key in atm_exp_hyd_keys}
        atm_exp_hyd_symb_dct = dict_.by_key({}, atm_exp_hyd_keys, fill_val='H')
        gra = add_atoms(gra, atm_exp_hyd_symb_dct)
        gra = add_bonds(gra, atm_exp_hyd_bnd_keys)
    return gra


def add_bonded_atom(gra, symb, atm_key, bnd_atm_key=None, imp_hyd=None,
                    atm_ste_par=None, bnd_ord=None, bnd_ste_par=None):
    """ Add a single atom and connect it to an atom already in the graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symb: The atomic symbol of the new atom
    :type symb: str
    :param atm_key: The key of the atom it will be connected to
    :type atm_key: int
    :param bnd_atm_key: The key of the new atom, defaults to None
    :type bnd_atm_key: int, optional
    :param imp_hyd: Implicit hydrogen count of the new atom, defaults to None
    :type imp_hyd: int, optional
    :param atm_ste_par: Stereo parity of the new atom, defaults to None
    :type atm_ste_par: bool or NoneType, optional
    :param bnd_ord: Order of the bond, defaults to None
    :type bnd_ord: int, optional
    :param bnd_ste_par: Stereo parity of the bond, defaults to None
    :type bnd_ste_par: bool or NoneType, optional
    :return: molecular graph
    :rtype: automol graph data structure
    """
    atm_keys = atom_keys(gra)

    bnd_atm_key = max(atm_keys) + 1 if bnd_atm_key is None else bnd_atm_key
    assert bnd_atm_key not in atm_keys, (
        f'Cannot add atom with key {bnd_atm_key}.\n'
        f'It is already in the graph:\n{gra}')

    symb_dct = {bnd_atm_key: symb}
    imp_hyd_dct = ({bnd_atm_key: imp_hyd} if imp_hyd is not None else None)
    atm_ste_par_dct = ({bnd_atm_key: atm_ste_par}
                       if atm_ste_par is not None else None)

    gra = add_atoms(gra, symb_dct, imp_hyd_dct=imp_hyd_dct,
                    ste_par_dct=atm_ste_par_dct)

    bnd_key = frozenset({bnd_atm_key, atm_key})
    bnd_ord_dct = {bnd_key: bnd_ord} if bnd_ord is not None else None
    bnd_ste_par_dct = ({bnd_key: bnd_ste_par}
                       if bnd_ste_par is not None else None)

    gra = add_bonds(gra, [bnd_key], ord_dct=bnd_ord_dct,
                    ste_par_dct=bnd_ste_par_dct)

    return gra


def shift_insert_bonded_atom(gra, sym, atm_key, bnd_atm_key=None, imp_hyd=None,
                             atm_ste_par=None, bnd_ord=None, bnd_ste_par=None):
    """ Insert a single atom and connect it to an atom already in the graph,
    shifting the keys after the new atom to make room for it

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symb: The atomic symbol of the new atom
    :type symb: str
    :param atm_key: The key of the atom it will be connected to
    :type atm_key: int
    :param bnd_atm_key: The key of the new atom, defaults to None
    :type bnd_atm_key: int, optional
    :param imp_hyd: Implicit hydrogen count of the new atom, defaults to None
    :type imp_hyd: int, optional
    :param atm_ste_par: Stereo parity of the new atom, defaults to None
    :type atm_ste_par: bool or NoneType, optional
    :param bnd_ord: Order of the bond, defaults to None
    :type bnd_ord: int, optional
    :param bnd_ste_par: Stereo parity of the bond, defaults to None
    :type bnd_ste_par: bool or NoneType, optional
    :return: molecular graph
    :rtype: automol graph data structure
    """
    keys = sorted(atom_keys(gra))
    bnd_atm_key_ = max(keys) + 1

    gra = add_bonded_atom(gra, sym, atm_key, bnd_atm_key=bnd_atm_key_,
                          imp_hyd=imp_hyd, atm_ste_par=atm_ste_par,
                          bnd_ord=bnd_ord, bnd_ste_par=bnd_ste_par)
    if bnd_atm_key != bnd_atm_key_:
        assert bnd_atm_key in keys
        idx = keys.index(bnd_atm_key)
        key_dct = {}
        key_dct.update({k: k for k in keys[:idx]})
        key_dct[bnd_atm_key_] = bnd_atm_key
        key_dct.update({k: k+1 for k in keys[idx:]})
        gra = relabel(gra, key_dct)

    return gra


def add_dummy_atoms(gra, dummy_key_dct):
    """ Insert dummy atoms to the graph, connecting them to existing atoms by
    null (order 0) bonds

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param dummy_key_dct: Keys of new dummy atoms, by atom that they connect to
    :type dummy_key_dct: dict[int: int]
    :return: molecular graph
    :rtype: automol graph data structure
    """
    atm_keys = atom_keys(gra)
    assert set(dummy_key_dct.keys()) <= atm_keys, (
        "Keys must be existing atoms in the graph.")
    assert not set(dummy_key_dct.values()) & atm_keys, (
        "Dummy atom keys cannot overlap with existing atoms.")

    for key, dummy_key in sorted(dummy_key_dct.items()):
        gra = add_bonded_atom(gra, 'X', key, bnd_atm_key=dummy_key, bnd_ord=0)

    return gra


def shift_insert_dummy_atoms(gra, dummy_key_dct):
    """ Insert dummy atoms to the graph, connecting them to existing atoms by
    null (order 0) bonds and shifting the keys after each new atom to make room
    for it

    Works with `dummy_key_dct` returned by `shift_remove_dummy_atoms`, undoing
    its operation.

    Note: The key values in `dummy_key_dct` are unavoidably confusing. All of
    the atom keys in `dummy_key_dct` refer to what they keys *will be* after
    insertion. For example, {0: 1, 2: 3} describes insertion of a dummy atom
    with key 1 after atom 0. The next key, 2, then refers to this atom's key
    *after* the dummy atoms are inserted. Originally, this atom had key 1.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param dummy_key_dct: Keys of new dummy atoms, by atom that they connect
        to, obtained from `shift_remove_dummy_atoms`
    :type dummy_key_dct: dict[int: int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    # These lists are used to track how the keys change:
    old_keys = sorted(atom_keys(gra))
    new_keys = sorted(atom_keys(gra))

    for key, dummy_key in reversed(sorted(dummy_key_dct.items(),
                                          key=lambda x: x[::-1])):
        key_dct = dict(zip(old_keys, new_keys))
        key_ = key_dct[key]

        gra = shift_insert_bonded_atom(
            gra, 'X', key_, bnd_atm_key=dummy_key, bnd_ord=0)

        idx = new_keys.index(dummy_key)
        new_keys[idx:] = numpy.add(new_keys[idx:], 1)
        new_keys.insert(idx, dummy_key)
        old_keys.insert(idx, None)

    return gra


def shift_remove_dummy_atoms(gra):
    """ Remove dummy atoms and standardize keys, returning the dummy key
    dictionary for converting back

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A molecular graph and a dummy key dictionary for converting back
        using `shift_insert_dummy_atoms()`
    :rtype: automol graph data structure
    """
    dummy_ngb_key_dct = dummy_atoms_neighbor_atom_key(gra)

    # These lists are used to track how the keys change:
    old_keys = sorted(atom_keys(gra))
    new_keys = sorted(atom_keys(gra))

    dummy_keys_dct = {}
    for dummy_key, key in sorted(dummy_ngb_key_dct.items()):
        key_dct = dict(zip(old_keys, new_keys))
        dummy_key_ = key_dct[dummy_key]
        key_ = key_dct[key]

        dummy_keys_dct[key_] = dummy_key_
        gra = remove_atoms(gra, [dummy_key], stereo=True)

        idx = old_keys.index(dummy_key)
        old_keys.pop(idx)
        new_keys.pop(idx)
        new_keys[idx:] = numpy.subtract(new_keys[idx:], 1)

    gra = standard_keys(gra)
    return gra, dummy_keys_dct


def without_pi_bonds(gra):
    """ Get a version of this graph without any pi-bonds

    All bond orders will be set to 1, except for dummy bonds (order 0) and
    reacting bonds (order 0.1 or 0.9).

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    bnd_keys = list(bond_keys(gra))
    # don't set dummy bonds to one!
    bnd_ord_dct = bond_orders(gra)
    bnd_ords = [0 if round(v, 1) == 0
                else round(v % 1, 1) if round(v % 1, 1) in (0.1, 0.9)
                else 1
                for v in map(bnd_ord_dct.__getitem__, bnd_keys)]
    bnd_ord_dct = dict(zip(bnd_keys, bnd_ords))
    return set_bond_orders(gra, bnd_ord_dct)


def without_dummy_atoms(gra):
    """ Remove dummy atoms from this graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    symb_dct = atom_symbols(gra)
    keys = [key for key, sym in symb_dct.items() if ptab.to_number(sym)]
    return subgraph(gra, keys, stereo=True)


def without_null_bonds(gra, except_dummies=True):
    """ Remove 0-order bonds from this graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param except_dummies: Keep 0-order bonds to dummy atoms?
    :type except_dummies: bool
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    dummy_atm_keys = atom_keys(gra, symb='X')
    ord_dct = dict_.filter_by_value(
        bond_orders(gra), func=lambda x: round(x, 1) == 0)
    null_bnd_keys = ord_dct.keys()
    if except_dummies:
        null_bnd_keys = [bk for bk in null_bnd_keys
                         if not any(k in dummy_atm_keys for k in bk)]
    gra = remove_bonds(gra, null_bnd_keys)
    return gra


def without_stereo(gra, atm_keys=None, bnd_keys=None):
    """ Remove stereo information (atom and bond parities) from this graph

    For TS graphs, fleeting and product stereochemistry will be removed as well

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_keys: Optionally, restrict this operation to a subset of bonds
    :type atm_keys: list[int]
    :param bnd_keys: Optionally, restrict this operation to a subset of bonds
    :type bnd_keys: list[frozenset[int]]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    if atm_keys is None and bnd_keys is not None:
        atm_keys = ()
    if bnd_keys is None and atm_keys is not None:
        bnd_keys = ()
    atm_keys = atom_keys(gra) if atm_keys is None else atm_keys
    bnd_keys = bond_keys(gra) if bnd_keys is None else bnd_keys
    atm_ste_par_dct = dict_.by_key({}, atm_keys, fill_val=None)
    bnd_ste_par_dct = dict_.by_key({}, bnd_keys, fill_val=None)
    gra = set_atom_stereo_parities(gra, atm_ste_par_dct)
    gra = set_bond_stereo_parities(gra, bnd_ste_par_dct)
    if is_ts_graph(gra):
        gra = ts_set_atom_product_stereo_parities(gra, atm_ste_par_dct)
        gra = ts_set_atom_fleeting_stereo_parities(gra, atm_ste_par_dct)
        gra = ts_set_bond_product_stereo_parities(gra, bnd_ste_par_dct)
        gra = ts_set_bond_fleeting_stereo_parities(gra, bnd_ste_par_dct)
    return gra


def explicit(gra, atm_keys=None):
    """ Make implicit hydrogens in this graph explicit

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_keys: Optionally, restrict this operation to a subset of atoms
    :type atm_keys: list[int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_keys = backbone_keys(gra) if atm_keys is None else atm_keys
    atm_keys = sorted(atm_keys)
    atm_imp_hyd_dct = dict_.by_key(atom_implicit_hydrogens(gra), atm_keys)

    atm_exp_hyd_keys_dct = {}
    next_atm_key = max(atom_keys(gra)) + 1
    for atm_key in atm_keys:
        imp_hyd = atm_imp_hyd_dct[atm_key]
        atm_exp_hyd_keys_dct[atm_key] = set(
            range(next_atm_key, next_atm_key+imp_hyd))
        next_atm_key += imp_hyd

    gra = set_atom_implicit_hydrogens(
        gra, dict_.by_key({}, atm_keys, fill_val=0))
    gra = add_atom_explicit_hydrogens(gra, atm_exp_hyd_keys_dct)
    return gra


def implicit(gra, atm_keys=None):
    """ Make explicit hydrogens in this graph implicit

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_keys: Optionally, restrict this operation to a subset of atoms
    :type atm_keys: list[int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_keys = backbone_keys(gra) if atm_keys is None else atm_keys

    atm_exp_hyd_keys_dct = dict_.by_key(
        atom_nonbackbone_hydrogen_keys(gra), atm_keys)

    inc_imp_hyd_keys_dct = dict_.transform_values(atm_exp_hyd_keys_dct, len)
    gra = change_implicit_hydrogens(gra, inc_imp_hyd_keys_dct)

    exp_hyd_keys = set(itertools.chain(*atm_exp_hyd_keys_dct.values()))
    gra = remove_atoms(gra, exp_hyd_keys, stereo=True)
    return gra


# # unions
def union(gra1, gra2, check=True):
    """ Get the union of two molecular graphs

    This is a disconnected graph consisting of the union of atom sets from each
    graph and the union of bond sets from each graph

    :param gra1: molecular graph
    :type gra1: automol graph data structure
    :param gra2: molecular graph
    :type gra2: automol graph data structure
    :param check: check that no keys overlap?
    :type check: bool
    :returns: molecular graph
    :rtype: automol graph data structure
    """
    if check:
        assert not atom_keys(gra1) & atom_keys(gra2)
    atm_dct = {}
    atm_dct.update(atoms(gra1))
    atm_dct.update(atoms(gra2))

    bnd_dct = {}
    bnd_dct.update(bonds(gra1))
    bnd_dct.update(bonds(gra2))
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def union_from_sequence(gras, check=True, shift_keys=False):
    """ Get the union of a sequence of graphs

    This is a disconnected graph consisting of the union of atom sets from each
    graph and the union of bond sets from each graph

    :param gras: A sequence of molecular graphs
    :param shift_keys: shift keys to prevent key overlap?
    :type shift_keys: bool
    :param check: check that no keys overlap?
    :type check: bool
    :returns: molecular graph
    :rtype: automol graph data structure
    """
    def _union(gra1, gra2):
        return union(gra1, gra2, check=check)

    if shift_keys:
        gras, _ = standard_keys_for_sequence(gras)

    return tuple(functools.reduce(_union, gras))


# # subgraphs and neighborhoods
def subgraph(gra, atm_keys, stereo=False):
    """ Get the subgraph of a subset of the atoms in this graph

    All bonds between the specified atoms will be included.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_keys: The atom keys to be included in the subgraph
    :type atm_keys: list[int]
    :param stereo: Keep stereo information in the subgraph?
    :type stereo: bool
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_keys = set(atm_keys)
    assert atm_keys <= atom_keys(gra)
    bnd_keys = set(filter(lambda x: x <= atm_keys, bond_keys(gra)))
    atm_dct = dict_.by_key(atoms(gra), atm_keys)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    sub = from_atoms_and_bonds(atm_dct, bnd_dct)
    if not stereo:
        sub = without_stereo(sub)
    return sub


def bond_induced_subgraph(gra, bnd_keys, stereo=False):
    """ Get the subgraph of a subset of bonds in this graph

    All atoms contained in the specified bonds will be included.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bnd_keys: The bond keys to be included in the subgraph
    :type bnd_keys: list[int]
    :param stereo: Keep stereo information in the subgraph?
    :type stereo: bool
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_keys = set(itertools.chain(*bnd_keys))
    bnd_keys = set(bnd_keys)
    assert atm_keys <= atom_keys(gra)
    atm_dct = dict_.by_key(atoms(gra), atm_keys)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    sub = from_atoms_and_bonds(atm_dct, bnd_dct)
    if not stereo:
        sub = without_stereo(sub)
    return sub


def atom_neighborhood(gra, atm_key, bnd_keys=None, stereo=False, ts_=True):
    """ Get the neighborhood subgraph of a specific atom

    The neighborhood subgraph contains the atom, its neighbors, and the bonds
    between them.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_key: the atom key
    :type atm_key: int
    :param bnd_keys: optionally, restrict this to a subset of the bond keys
    :type bnd_keys: tuple[frozenset[int]]
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    gra = gra if ts_ else ts_reactants_graph(gra)
    bnd_keys = (bond_keys(gra, ts_=ts_)
                if bnd_keys is None else bnd_keys)
    nbh_bnd_keys = set(k for k in bnd_keys if atm_key in k)
    nbh = bond_induced_subgraph(gra, nbh_bnd_keys, stereo=stereo)
    return nbh


def atom_neighborhoods(gra, bnd_keys=None, stereo=False, ts_=True):
    """ Get the neighborhood subgraphs of each atom, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bnd_keys: optionally, restrict this to a subset of the bond keys
    :type bnd_keys: tuple[frozenset[int]]
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: neighborhood subgraphs, by atom key
    :rtype: dict[int: automol graph data structure]
    """
    bnd_keys = (bond_keys(gra, ts_=ts_)
                if bnd_keys is None else bnd_keys)

    def _neighborhood(atm_key):
        return atom_neighborhood(gra, atm_key, bnd_keys=bnd_keys,
                                 stereo=stereo, ts_=ts_)

    atm_keys = list(atom_keys(gra))
    atm_nbh_dct = dict(zip(atm_keys, map(_neighborhood, atm_keys)))
    return atm_nbh_dct


def bond_neighborhood(gra, bnd_key, bnd_keys=None, stereo=False, ts_=True):
    """ Get the neighborhood subgraph of a specific bond

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bnd_key: the bond key
    :type bnd_key: frozenset[int]
    :param bnd_keys: optionally, restrict this to a subset of the bond keys
    :type bnd_keys: tuple[frozenset[int]]
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    bnd_keys = (bond_keys(gra, ts_=ts_)
                if bnd_keys is None else bnd_keys)
    nbh_bnd_keys = set(filter(lambda x: bnd_key & x, bnd_keys))
    nbh = bond_induced_subgraph(gra, nbh_bnd_keys, stereo=stereo)
    return nbh


def bond_neighborhoods(gra, bnd_keys=None, stereo=False, ts_=True):
    """ Get the neighborhood subgraphs of each bond, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bnd_keys: optionally, restrict this to a subset of the bond keys
    :type bnd_keys: tuple[frozenset[int]]
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: neighborhood subgraphs, by bond key
    :rtype: dict[frozenset: automol graph data structure]
    """
    bnd_keys = list(bond_keys(gra, ts_=ts_) if bnd_keys is None else bnd_keys)

    def _neighborhood(bnd_key):
        return bond_neighborhood(gra, bnd_key, bnd_keys=bnd_keys,
                                 stereo=stereo, ts_=ts_)

    bnd_nbh_dct = dict(zip(bnd_keys, map(_neighborhood, bnd_keys)))
    return bnd_nbh_dct


def atom_neighbor_atom_keys(gra, atm_key, bnd_keys=None, symb=None,
                            excl_symbs=(), ts_=True):
    """ Get keys for an atom's neighbors

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_key: the atom key
    :type atm_key: int
    :param bnd_keys: optionally, restrict this to a subset of the bond keys
    :type bnd_keys: tuple[frozenset[int]]
    :param symb: Optionally, restrict this to atoms with a particular
        atomic symbol (e.g., 'H' for hydrogens).
    :type symb: str
    :param excl_symbs: Optionally, exclude atoms with particular atomic
        symbols.
    :type excl_symbs: tuple[str]
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: The keys of neighboring atoms
    :rtype: frozenset[int]
    """
    atm_nbh = atom_neighborhood(
        gra, atm_key, bnd_keys=bnd_keys, ts_=ts_)
    atm_nbh_keys = atom_keys(atm_nbh, symb=symb, excl_symbs=excl_symbs)
    atm_ngb_keys = frozenset(atm_nbh_keys - {atm_key})
    return atm_ngb_keys


def atoms_neighbor_atom_keys(gra, ts_=True):
    """ Get the keys for each atom's neighboring atoms, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: Neighboring atom keys by atom, as a dictionary
    :rtype: dict[int: frozenset]
    """
    def _neighbor_keys(atm_key, atm_nbh):
        return frozenset(atom_keys(atm_nbh) - {atm_key})

    atm_ngb_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(gra, ts_=ts_), _neighbor_keys)
    return atm_ngb_keys_dct


def atoms_sorted_neighbor_atom_keys(gra, symbs_first=('C',), symbs_last=('H',),
                                    ords_last=(0.1,), prioritize_keys=(),
                                    ts_=True):
    """ Get keys for each atom's neighbors, sorted in a particular order

    :param gra: the graph
    :param symbs_first: Atom types to put first, defaults to ('C',)
    :type symbs_first: tuple, optional
    :param symbs_last: Atom types to put last, defaults to ('H',)
    :type symbs_last: tuple, optional
    :param ords_last: Put neighbors bonded with these bond orders last,
        defaults to (0.1,). (Mainly for internal use.)
    :type ords_last: tuple, optional
    :param prioritize_keys: Keys to put first no matter what, defaults to ()
    :type prioritize_keys: tuple, optional
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: Sorted neighboring atom keys by atom, as a dictionary
    :rtype: dict[int: tuple]
    """
    atm_symb_dct = atom_symbols(gra)
    bnd_ord_dct = bond_orders(gra)

    def _neighbor_keys(atm_key, atm_nbh):
        keys = sorted(atom_keys(atm_nbh) - {atm_key})
        bnd_keys = [frozenset({atm_key, k}) for k in keys]
        ords = list(map(bnd_ord_dct.__getitem__, bnd_keys))
        ords = [-1 if o not in ords_last else ords_last.index(o)
                for o in ords]
        symbs = list(map(atm_symb_dct.__getitem__, keys))
        pris = [0 if k in prioritize_keys else 1 for k in keys]
        srt_vals = list(zip(ords, pris, symbs))
        srt = automol.formula.argsort_symbols(
            srt_vals, symbs_first, symbs_last, idx=2)
        keys = tuple(map(keys.__getitem__, srt))
        return keys

    atm_ngb_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(gra, ts_=ts_), _neighbor_keys)
    return atm_ngb_keys_dct


def atom_sorted_neighbor_atom_keys(gra, atm_key, excl_atm_keys=(),
                                   incl_atm_keys=None, symbs_first=('C',),
                                   symbs_last=('H',), ts_=True):
    """ Get keys for this atom's neighbors, sorted in a particular order

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_key: The atom key
    :type atm_key: int
    :param excl_atm_keys: Atom keys to exclude, defaults to ()
    :type excl_atm_keys: tuple, optional
    :param incl_atm_keys: Restrict the search to a subset of atom keys,
        defaults to None
    :type incl_atm_keys: tuple, optional
    :param symbs_first: Atom types to put first, defaults to ('C',)
    :type symbs_first: tuple, optional
    :param symbs_last: Atom types to put last, defaults to ('H',)
    :type symbs_last: tuple, optional
    :return: The neighboring atom keys, in the requested order
    :rtype: tuple[int]
    """
    atm_symb_dct = atom_symbols(gra)
    incl_atm_keys = atom_keys(gra) if incl_atm_keys is None else incl_atm_keys

    atm_nbh = atom_neighborhood(gra, atm_key, ts_=ts_)
    atm_keys = sorted(atom_keys(atm_nbh) - {atm_key} - set(excl_atm_keys))
    atm_keys = [k for k in atm_keys if k in incl_atm_keys]

    symbs = list(map(atm_symb_dct.__getitem__, atm_keys))
    srt = automol.formula.argsort_symbols(symbs, symbs_first, symbs_last)
    atm_keys = tuple(map(atm_keys.__getitem__, srt))
    return atm_keys


def atom_neighbor_atom_key(gra, atm_key, excl_atm_keys=(), incl_atm_keys=None,
                           symbs_first=('C',), symbs_last=('H',)):
    """ Get a key for one of an atom's neighbors

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_key: The atom key
    :type atm_key: int
    :param excl_atm_keys: Atom keys to exclude, defaults to ()
    :type excl_atm_keys: tuple, optional
    :param incl_atm_keys: Restrict the search to a subset of atom keys,
        defaults to None
    :type incl_atm_keys: tuple, optional
    :param symbs_first: Atom symbols to search for first, defaults to ('C',)
    :type symbs_first: tuple, optional
    :param symbs_last: Atom symbols to search for last, defaults to ('H',)
    :type symbs_last: tuple, optional
    :return: The neighboring atom key
    :rtype: int
    """
    atm_keys = atom_sorted_neighbor_atom_keys(
        gra, atm_key, excl_atm_keys=excl_atm_keys, incl_atm_keys=incl_atm_keys,
        symbs_first=symbs_first, symbs_last=symbs_last)
    return atm_keys[0] if atm_keys else None


def atom_bond_keys(gra, atm_key, ts_=True):
    """ Get the bond keys of a specific atom

    This is equivalent to getting the bond keys of this atom's neighborhood

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: The bond keys containing this atom
    :rtype: frozenset[frozenset[int]]
    """
    return bond_keys(atom_neighborhood(gra, atm_key, ts_=ts_))


def atoms_bond_keys(gra, ts_=True):
    """ Get the bond keys of each atom, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: The bond keys of each atom, as a dictionary
    :rtype: dict[int: frozenset]
    """
    atm_nbhs = atom_neighborhoods(gra, ts_=ts_)
    return dict_.transform_values(atm_nbhs, bond_keys)


def dummy_atoms_neighbor_atom_key(gra, ts_=True):
    """ Get the atoms that are connected to dummy atoms, by dummy atom key
    (Requires that each dummy atom only be connected to one neighbor)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: The atoms that are connected to dummy atoms, by dummy atom key
    :rtype: dict[int: int]
    """
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra, ts_=ts_)
    dummy_atm_keys = atom_keys(gra, symb='X')

    dummy_ngb_key_dct = {}
    for key in dummy_atm_keys:
        ngb_keys = atm_ngb_keys_dct[key]
        assert len(ngb_keys) == 1, (
            "Dummy atoms should only be connected to one atom!")
        ngb_key, = ngb_keys
        dummy_ngb_key_dct[key] = ngb_key

    return dummy_ngb_key_dct


def bonds_neighbor_atom_keys(gra, ts_=True):
    """ Get the keys of each bond's neighboring atoms, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: Neighboring atom keys by bond, as a dictionary
    :rtype: dict[frozenset: frozenset]
    """
    def _neighbor_keys(bnd_key, bnd_nbh):
        return frozenset(atom_keys(bnd_nbh) - bnd_key)

    bnd_ngb_keys_dct = dict_.transform_items_to_values(
        bond_neighborhoods(gra, ts_=ts_), _neighbor_keys)
    return bnd_ngb_keys_dct


def bonds_neighbor_bond_keys(gra, ts_=True):
    """ Get the keys of each bond's neighboring bonds, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: Neighboring bond keys by bond, as a dictionary
    :rtype: dict[frozenset: frozenset]
    """
    def _neighbor_keys(bnd_key, bnd_nbh):
        bnd_keys = bond_keys(bnd_nbh)
        bnd_keys -= {bnd_key}
        bnd_keys = frozenset(key for key in bnd_keys if key & bnd_key)
        return bnd_keys

    bnd_ngb_keys_dct = dict_.transform_items_to_values(
        bond_neighborhoods(gra, ts_=ts_), _neighbor_keys)
    return bnd_ngb_keys_dct
