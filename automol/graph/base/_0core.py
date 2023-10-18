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
import functools
import itertools
import numbers
from typing import Dict, List, Tuple

import numpy
import yaml
from phydat import phycon, ptab

import automol.formula
import automol.util.dict_.multi as mdict
from automol import util
from automol.util import ZmatConv, dict_, zmat_conv

ATM_SYM_POS = 0
ATM_IMP_HYD_POS = 1
ATM_STE_PAR_POS = 2

BND_ORD_POS = 0
BND_STE_PAR_POS = 1

ATM_PROP_NAMES = ("symbol", "implicit_hydrogens", "stereo_parity")
BND_PROP_NAMES = ("order", "stereo_parity")


# # constructors
def from_data(
    atm_symb_dct,
    bnd_keys,
    atm_imp_hyd_dct=None,
    atm_ste_par_dct=None,
    bnd_ord_dct=None,
    bnd_ste_par_dct=None,
):
    """Construct a molecular graph from data

    Graph data structure:
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

    Bonds to dummy atoms will automatically be given an order of 0.

    :param atm_symb_dct: atomic symbols, by atom key
    :type atm_symb_dct: dict
    :param bnd_keys: bond keys
    :type bnd_keys: set
    :param atm_imp_hyd_dct: the number of implicit hydrogens associated
        with each atom, by atom key
    :type atm_imp_hyd_dct: dict
    :param atm_ste_par_dct: stereo parities, by atom key;
    :type atm_ste_par_dct: dict
    :param bnd_ord_dct: bond orders, by bond key
    :type bnd_ord_dct: dict
    :param bnd_ste_par_dct: stereo parities, by bond key;
    :type bnd_ste_par_dct: dict
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    dummy_keys = dict_.keys_by_value(atm_symb_dct, lambda s: s.upper() == "X")

    atm_dct = atoms_from_data(
        atm_symb_dct=atm_symb_dct,
        atm_imp_hyd_dct=atm_imp_hyd_dct,
        atm_ste_par_dct=atm_ste_par_dct,
    )

    bnd_dct = bonds_from_data(
        bnd_keys=bnd_keys,
        bnd_ord_dct=bnd_ord_dct,
        bnd_ste_par_dct=bnd_ste_par_dct,
        dummy_atm_keys=dummy_keys,
    )

    atm_keys = set(atm_dct.keys())
    bnd_keys = set(bnd_dct.keys())

    assert all(
        bnd_key <= atm_keys for bnd_key in bnd_keys
    ), f"\natm_keys: {atm_keys}\nbnd_keys: {bnd_keys}"

    return (atm_dct, bnd_dct)


def atoms_from_data(atm_symb_dct, atm_imp_hyd_dct=None, atm_ste_par_dct=None):
    """Construct an atom dictionary from constituent data.

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
    :param atm_ste_par_dct: stereo parities, by atom key;
    :type atm_ste_par_dct: dict
    :returns: The atoms and their associated properties, as a dictionary
    :rtype: dict[int: tuple]
    """
    keys = sorted(atm_symb_dct.keys())
    symbs = dict_.values_by_key(atm_symb_dct, keys)
    hyds = dict_.values_by_key(dict_.empty_if_none(atm_imp_hyd_dct), keys, fill_val=0)
    pars = dict_.values_by_key(
        dict_.empty_if_none(atm_ste_par_dct), keys, fill_val=None
    )

    natms = len(keys)

    assert len(symbs) == natms
    assert len(hyds) == natms
    assert len(pars) == natms
    assert all(par in (None, False, True) for par in pars)

    symbs = list(map(ptab.to_symbol, symbs))
    hyds = list(map(int, hyds))
    pars = [bool(par) if par is not None else par for par in pars]

    atm_dct = dict(zip(keys, zip(symbs, hyds, pars)))

    return atm_dct


def bonds_from_data(
    bnd_keys, bnd_ord_dct=None, bnd_ste_par_dct=None, dummy_atm_keys=()
):
    """Construct a bond dictionary from constituent data.

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
    :param bnd_ste_par_dct: stereo parities, by bond key;
    :type bnd_ste_par_dct: dict
    :param dummy_atm_keys: keys to dummy atoms, resulting in 0-order bonds for these
    :type dummy_atm_keys: list
    :returns: The bonds and their associated properties, as a dictionary
    :rtype: dict[frozenset({int, int}): tuple]
    """
    dummy_atm_keys = set(dummy_atm_keys)
    keys = sorted(bnd_keys)
    assert all(len(key) == 2 for key in keys)
    ords = dict_.values_by_key(dict_.empty_if_none(bnd_ord_dct), keys, fill_val=1)
    pars = dict_.values_by_key(
        dict_.empty_if_none(bnd_ste_par_dct), keys, fill_val=None
    )

    nbnds = len(keys)

    assert len(ords) == nbnds
    assert len(pars) == nbnds
    assert all(par in (None, False, True) for par in pars)

    keys = list(map(frozenset, keys))
    # If ts_ = False, we should assert that round(o) == o
    dummy_keys = [k for k in keys if k & dummy_atm_keys]
    ords = [
        0 if k in dummy_keys else int(o) if round(o) == o else float(round(o, 1))
        for k, o in zip(keys, ords)
    ]
    pars = [bool(par) if par is not None else par for par in pars]

    bnd_dct = dict(zip(keys, zip(ords, pars)))

    return bnd_dct


def from_atoms_and_bonds(atm_dct, bnd_dct):
    """Construct a molecular graph from atom and bond dictionaries.

    Data structure:
        gra = (atm_dct, bnd_dct)

    :param atm_dct: atom dictionary
    :type atm_dct: dict[int: tuple]
    :param bnd_dct: bond dictionary
    :type bnd_dct: dict[frozenset({int, int}): tuple]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """

    atm_dct = dict(atm_dct)
    bnd_dct = dict(bnd_dct)

    atm_nprops = dict_.transform_values(atm_dct, len).values()
    bnd_nprops = dict_.transform_values(bnd_dct, len).values()

    # Require that it has the right format for a non-TS graph and set
    # `ts_` to `False`
    assert all(n == 3 for n in atm_nprops) and all(n == 2 for n in bnd_nprops), (
        f"Atom or bond dictionary has improper format:\n"
        f"atm_dct:\n{atm_dct}\nbnd_dct:\n{bnd_dct}\n"
    )

    atm_nprops = dict_.transform_values(atm_dct, len)
    bnd_nprops = dict_.transform_values(bnd_dct, len)

    atm_keys = sorted(atm_dct.keys())
    atm_symb_dct = mdict.by_key_by_position(atm_dct, atm_keys, ATM_SYM_POS)
    atm_imp_hyd_dct = mdict.by_key_by_position(atm_dct, atm_keys, ATM_IMP_HYD_POS)
    atm_ste_par_dct = mdict.by_key_by_position(atm_dct, atm_keys, ATM_STE_PAR_POS)

    bnd_keys = sorted(bnd_dct.keys())
    bnd_ord_dct = mdict.by_key_by_position(bnd_dct, bnd_keys, BND_ORD_POS)
    bnd_ste_par_dct = mdict.by_key_by_position(bnd_dct, bnd_keys, BND_STE_PAR_POS)

    return from_data(
        atm_symb_dct,
        bnd_dct.keys(),
        atm_imp_hyd_dct=atm_imp_hyd_dct,
        atm_ste_par_dct=atm_ste_par_dct,
        bnd_ord_dct=bnd_ord_dct,
        bnd_ste_par_dct=bnd_ste_par_dct,
    )


# # getters
def atoms(gra):
    """Get the atoms of this graph, along with their associated properties

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atoms and their associated properties, as a dictionary
    :rtype: dict[int: tuple]
    """
    atm_dct, _ = gra
    return atm_dct


def bonds(gra):
    """Get the bonds of this graph, along with their associated properties

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The bonds and their associated properties, as a dictionary
    :rtype: dict[frozenset({int, int}): tuple]
    """
    _, bnd_dct = gra
    return bnd_dct


def atom_keys(gra, symb=None, excl_symbs=()):
    """Get the atom keys of this molecular graph

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
    symb_dct = atom_symbols(gra)
    atm_keys = frozenset(k for k in atoms(gra).keys() if symb_dct[k] not in excl_symbs)
    if symb is not None:
        atm_keys = frozenset(k for k in atm_keys if symb_dct[k] == symb)
    return atm_keys


def bond_keys(gra, ts_=True):
    """Get the bond keys of this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: The bond keys
    :rtype: frozenset[{int, int}]
    """
    gra = gra if ts_ else ts_reagents_graph_without_stereo(gra)
    return frozenset(bonds(gra).keys())


def atom_symbols(gra, dummy_symbol=None):
    """Get the atom symbols of this molecular graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param dummy_symbol: Use this symbol for dummy atoms, defaults to None
        (If `None`, the usual symbol 'X' will be used)
    :type dummy_symbol: str
    :returns: A dictionary of atomic symbols, by atom key
    :rtype: dict[int: str]
    """
    symb_dct = mdict.by_key_by_position(atoms(gra), atoms(gra).keys(), ATM_SYM_POS)
    if dummy_symbol is not None:
        symb_dct = dict_.transform_values(
            symb_dct, lambda s: dummy_symbol if s == "X" else s
        )
    return symb_dct


def bond_orders(gra, ts_=True):
    """Get the bond orders of this molecular graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: A dictionary of bond orders, by bond key
    :rtype: dict[frozenset: int or float]
    """
    gra = gra if ts_ else ts_reagents_graph_without_stereo(gra)
    return mdict.by_key_by_position(bonds(gra), bonds(gra).keys(), BND_ORD_POS)


def atom_implicit_hydrogens(gra):
    """Get the implicit hydrogen valences of atoms in this molecular graph, as
        a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of implicit hydrogen valences, by atom key
    :rtype: dict[int: int]
    """
    return mdict.by_key_by_position(atoms(gra), atoms(gra).keys(), ATM_IMP_HYD_POS)


def atom_stereo_parities(gra):
    """Get the atom stereo parities of this molecular graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of atom stereo parities, by atom key
    :rtype: dict[int: bool or NoneType]
    """
    ret = mdict.by_key_by_position(atoms(gra), atoms(gra).keys(), ATM_STE_PAR_POS)
    return ret


def bond_stereo_parities(gra):
    """Get the bond stereo parities of this molecular graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of bond stereo parities, by bond key
    :rtype: dict[frozenset: bool or NoneType]
    """
    ret = mdict.by_key_by_position(bonds(gra), bonds(gra).keys(), BND_STE_PAR_POS)
    return ret


def stereo_parities(gra):
    """Get the atom and bond stereo parities of this molecular graph, as a
        single dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of stereo parities, by atom/bond key
    :rtype: dict[int or frozenset: bool or NoneType]
    """
    par_dct = atom_stereo_parities(gra)
    par_dct.update(bond_stereo_parities(gra))
    return par_dct


# # TS graph constructor
def ts_graph(gra, frm_bnd_keys, brk_bnd_keys):
    """Construct a TS graph from a molecular graph

    :param gra: molecular graph, representing the reactants
    :type gra: automol graph data structure
    :param frm_bnd_keys: Keys to bonds which are forming in the TS
    :type frm_bnd_keys: tuple[frozenset]
    :param brk_bnd_keys: Keys to bonds which are breaking in the TS
    :type brk_bnd_keys: tuple[frozenset]
    :returns: TS graph
    :rtype: automol graph data structure
    """
    frm_bnd_keys = frozenset(map(frozenset, frm_bnd_keys))
    brk_bnd_keys = frozenset(map(frozenset, brk_bnd_keys))

    frm_ord_dct = {k: 0.1 for k in frm_bnd_keys}
    brk_ord_dct = {k: 0.9 for k in brk_bnd_keys}

    gra = add_bonds(gra, frm_bnd_keys, ord_dct=frm_ord_dct, check=False)
    gra = add_bonds(gra, brk_bnd_keys, ord_dct=brk_ord_dct, check=False)
    return gra


# # TS graph getters
def ts_forming_bond_keys(tsg):
    """Get the forming bonds from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: The keys to forming bonds
    :rtype: frozenset[frozenset[{int, int}]]
    """
    ord_dct = bond_orders(tsg)
    frm_bnd_keys = [k for k, o in ord_dct.items() if round(o, 1) == 0.1]
    return frozenset(map(frozenset, frm_bnd_keys))


def ts_breaking_bond_keys(tsg):
    """Get the breaking bonds from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: The keys to breaking bonds
    :rtype: frozenset[frozenset[{int, int}]]
    """
    ord_dct = bond_orders(tsg)
    brk_bnd_keys = [k for k, o in ord_dct.items() if round(o % 1, 1) == 0.9]
    return frozenset(map(frozenset, brk_bnd_keys))


def ts_reacting_bond_keys(tsg):
    """Get all of the bonds involved in the reaction

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: The keys to bonds involved in the reaction
    :rtype: frozenset[frozenset[{int, int}]]
    """
    bnd_keys = ts_forming_bond_keys(tsg) | ts_breaking_bond_keys(tsg)
    return bnd_keys


def ts_reacting_atom_keys(tsg):
    """Get all of the atoms involved in the reaction

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: The keys to atoms involved in the reaction
    :rtype: frozenset[int]
    """
    bnd_keys = ts_reacting_bond_keys(tsg)
    atm_keys = frozenset(itertools.chain(*bnd_keys))
    return atm_keys


def ts_transferring_atoms(tsg) -> Dict[int, Tuple[int, int]]:
    """Get get transferring atoms, along with their donors and acceptors

    An atom is "transferring" if it is breaking and forming exactly one bond:

        D---T...A
         (b) (f)

    :param tsg: _description_
    :type tsg: _type_
    :returns: A dictionary mapping each transferring atom onto its donor (what it is
        transferring from) and its acceptor (what it is transferring to)
    :rtype: Dict[int, Tuple[int, int]]
    """
    bkeys_dct = atoms_bond_keys(tsg)
    rxn_akeys = ts_reacting_atom_keys(tsg)
    frm_bkeys = ts_forming_bond_keys(tsg)
    brk_bkeys = ts_breaking_bond_keys(tsg)

    tra_dct = {}
    for akey in rxn_akeys:
        afrm_bkeys = bkeys_dct[akey] & frm_bkeys
        abrk_bkeys = bkeys_dct[akey] & brk_bkeys
        if len(afrm_bkeys) == len(abrk_bkeys) == 1:
            (frm_bkey,) = afrm_bkeys
            (brk_bkey,) = abrk_bkeys
            (tra_key,) = brk_bkey & frm_bkey
            (don_key,) = brk_bkey - frm_bkey
            (acc_key,) = frm_bkey - brk_bkey
            tra_dct[tra_key] = (don_key, acc_key)

    return tra_dct


def ts_reverse(tsg):
    """Reverse the direction of a TS graph

    Since the parities are relative to the canonical keys *in the canonical
    direction*, no change to the parities occurs upon reversal.

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A TS graph for the reverse reaction
    :rtype: automol graph data structure
    """
    assert without_pi_bonds(tsg) == tsg, f"TS graphs shouldn't have pi bonds:\n{tsg}"
    ord_dct = bond_orders(tsg)
    rev_ord_dct = {
        k: (0.1 if round(o, 1) == 0.9 else 0.9 if round(o, 1) == 0.1 else o)
        for k, o in ord_dct.items()
    }
    rev_tsg = set_bond_orders(tsg, rev_ord_dct)

    return rev_tsg


def ts_reagents_graph_without_stereo(
    tsg, prod=False, dummy=True, keep_zeros=False, keep_stereo=False
):
    """Get the reactants or products from a TS graph, without stereo

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param prod: Replace reacting bond orders with product values instead?
    :type prod: bool
    :param dummy: Keep dummy atoms? default True
    :type dummy: bool, optional
    :param keep_zeros: Keep the bonds with a resulting bond order of 0?
    :type keep_zeros: bool
    :param keep_stereo: Keep stereo, even though it is invalid?
    :type keep_stereo: bool
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    bnd_keys = ts_reacting_bond_keys(tsg)
    ord_dct = dict_.by_key(bond_orders(tsg), bnd_keys)
    if prod:
        ord_dct = dict_.transform_values(ord_dct, lambda o: 1 - round(o, 1))
    ord_dct = dict_.transform_values(ord_dct, round)

    # Remove dummy atoms, unless requested to keep
    if not dummy:
        tsg = without_dummy_atoms(tsg)

    # Round the bond orders for forming bonds, and remove forming bonds
    tsg = set_bond_orders(tsg, ord_dct)
    if not keep_zeros:
        tsg = without_bonds_by_orders(tsg, ords=0, skip_dummies=True)

    # Remove invalid stereo, unless requested otherwise
    if not keep_stereo:
        tsg = without_stereo(tsg)

    return tsg


# # setters
def set_atom_symbols(gra, atm_symb_dct):
    """Set the atom symbols of this molecular graph with a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_symb_dct: A dictionary of atomic symbols, by atom key.
    :type atm_symb_dct: dict[int: str]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_dct = mdict.set_by_key_by_position(atoms(gra), atm_symb_dct, ATM_SYM_POS)
    return from_atoms_and_bonds(atm_dct, bonds(gra))


def set_bond_orders(gra, bnd_ord_dct):
    """Set the bond orders of this molecular graph with a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bnd_ord_dct: A dictionary of bond orders, by bond key.
    :type bnd_ord_dct: dict[frozenset: int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    bnd_dct = mdict.set_by_key_by_position(bonds(gra), bnd_ord_dct, BND_ORD_POS)
    return from_atoms_and_bonds(atoms(gra), bnd_dct)


def set_atom_implicit_hydrogens(gra, atm_imp_hyd_dct):
    """Set the implicit hydrogen valences of atoms in this molecular graph, as
        a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_imp_hyd_dct: A dictionary of implicit hydrogen valences, by
        atom key
    :type atm_imp_hyd_dct: dict[int: int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_dct = mdict.set_by_key_by_position(atoms(gra), atm_imp_hyd_dct, ATM_IMP_HYD_POS)
    bnd_dct = bonds(gra)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def set_atom_stereo_parities(gra, atm_par_dct):
    """Set the atom stereo parities of this molecular graph with a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_par_dct: A dictionary of atom stereo parities, by atom key
    :type atm_par_dct: dict[int: bool or NoneType]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_dct = mdict.set_by_key_by_position(atoms(gra), atm_par_dct, ATM_STE_PAR_POS)
    return from_atoms_and_bonds(atm_dct, bonds(gra))


def set_bond_stereo_parities(gra, bnd_par_dct):
    """Set the bond stereo parities of this molecular graph with a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bnd_par_dct: A dictionary of bond stereo parities, by bond key
    :type bnd_par_dct: dict[frozenset: bool or NoneType]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    bnd_dct = mdict.set_by_key_by_position(bonds(gra), bnd_par_dct, BND_STE_PAR_POS)
    return from_atoms_and_bonds(atoms(gra), bnd_dct)


def set_stereo_parities(gra, par_dct):
    """Set the atom and bond stereo parities of this molecular graph with a
        single dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param par_dct: A dictionary of stereo parities, by atom/bond key
    :type par_dct: dict[int or frozenset: bool or NoneType]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_par_dct = dict_.filter_by_key(par_dct, lambda x: isinstance(x, numbers.Number))
    bnd_par_dct = dict_.filter_by_key(
        par_dct, lambda x: not isinstance(x, numbers.Number)
    )
    gra = set_atom_stereo_parities(gra, atm_par_dct)
    gra = set_bond_stereo_parities(gra, bnd_par_dct)
    return gra


# # I/O
def string(gra, one_indexed=True):
    """Generate a string representation of the graph, in YAML format

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param one_indexed: Switch to one-indexing for keys?
    :type one_indexed: bool
    :returns: A string representation of the molecular graph
    :rtype: str
    """
    yaml_gra_dct = yaml_data(gra, one_indexed=one_indexed)
    gra_str = yaml.dump(yaml_gra_dct, default_flow_style=None, sort_keys=False)
    return gra_str


def yaml_data(gra, one_indexed=True):
    """Generate a YAML dictionary representation of the graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param one_indexed: Switch to one-indexing for keys?
    :type one_indexed: bool
    :returns: A YAML-friendly dictionary representation of the graph
    :rtype: dict
    """
    if one_indexed:
        # shift to one-indexing when we print
        atm_key_dct = {atm_key: atm_key + 1 for atm_key in atom_keys(gra)}
        gra = relabel(gra, atm_key_dct)

    yaml_atm_dct = atoms(gra)
    yaml_bnd_dct = bonds(gra)

    ts_ = is_ts_graph(gra)
    atm_nprops = 5 if ts_ else 3
    bnd_nprops = 4 if ts_ else 2

    # prepare the atom dictionary
    yaml_atm_dct = dict(sorted(yaml_atm_dct.items()))
    yaml_atm_dct = dict_.transform_values(
        yaml_atm_dct, lambda x: dict(zip(ATM_PROP_NAMES[:atm_nprops], x))
    )

    # perpare the bond dictionary
    yaml_bnd_dct = dict_.transform_keys(yaml_bnd_dct, lambda x: tuple(sorted(x)))
    yaml_bnd_dct = dict(sorted(yaml_bnd_dct.items()))
    yaml_bnd_dct = dict_.transform_keys(yaml_bnd_dct, lambda x: "-".join(map(str, x)))
    yaml_bnd_dct = dict_.transform_values(
        yaml_bnd_dct, lambda x: dict(zip(BND_PROP_NAMES[:bnd_nprops], x))
    )

    yaml_gra_dct = {"atoms": yaml_atm_dct, "bonds": yaml_bnd_dct}
    return yaml_gra_dct


def from_string(gra_str, one_indexed=True):
    """Generate a graph from a string representation, in YAML format

    :param gra_str: A string representation of the molecular graph
    :type gra_str: str
    :param one_indexed: Assume one-indexing for string keys?
    :type one_indexed: bool
    :returns: molecular graph
    :rtype: automol graph data structure
    """
    yaml_gra_dct = yaml.load(gra_str, Loader=yaml.FullLoader)
    gra = from_yaml_data(yaml_gra_dct, one_indexed=one_indexed)
    return gra


def from_yaml_data(yaml_gra_dct, one_indexed=True):
    """Generate a graph from a YAML dictionary representation

    :param yaml_gra_dct: A YAML-friendly dictionary representation of the graph
    :type yaml_gra_dct: dict
    :param one_indexed: Assume one-indexing for YAML dict keys?
    :type one_indexed: bool
    :returns: molecular graph
    :rtype: automol graph data structure
    """
    atm_dct = yaml_gra_dct["atoms"]
    bnd_dct = yaml_gra_dct["bonds"]

    atm_dct = dict_.transform_values(
        atm_dct, lambda x: tuple(map(x.__getitem__, ATM_PROP_NAMES[: len(x)]))
    )

    bnd_dct = dict_.transform_keys(bnd_dct, lambda x: frozenset(map(int, x.split("-"))))

    bnd_dct = dict_.transform_values(
        bnd_dct, lambda x: tuple(map(x.__getitem__, BND_PROP_NAMES[: len(x)]))
    )

    gra = from_atoms_and_bonds(atm_dct, bnd_dct)

    if one_indexed:
        # revert one-indexing if the input is one-indexed
        atm_key_dct = {atm_key: atm_key - 1 for atm_key in atom_keys(gra)}
        gra = relabel(gra, atm_key_dct)

    return gra


def from_old_yaml_data(yaml_gra_dct, one_indexed=True):
    """Generate a graph from a YAML dictionary representation

    :param yaml_gra_dct: A YAML-friendly dictionary representation of the graph
    :type yaml_gra_dct: dict
    :param one_indexed: Assume one-indexing for YAML dict keys?
    :type one_indexed: bool
    :returns: molecular graph
    :rtype: automol graph data structure
    """
    old_atm_prop_names = ("symbol", "implicit_hydrogen_valence", "stereo_parity")

    atm_dct = yaml_gra_dct["atoms"]
    bnd_dct = yaml_gra_dct["bonds"]

    if "implicit_hydrogens" in list(atm_dct.values())[0]:
        atm_dct = dict_.transform_values(
            atm_dct, lambda x: tuple(map(x.__getitem__, ATM_PROP_NAMES[: len(x)]))
        )
    else:
        atm_dct = dict_.transform_values(
            atm_dct, lambda x: tuple(map(x.__getitem__, old_atm_prop_names[: len(x)]))
        )

    bnd_dct = dict_.transform_keys(bnd_dct, lambda x: frozenset(map(int, x.split("-"))))

    bnd_dct = dict_.transform_values(
        bnd_dct, lambda x: tuple(map(x.__getitem__, BND_PROP_NAMES[: len(x)]))
    )

    gra = from_atoms_and_bonds(atm_dct, bnd_dct)

    if one_indexed:
        # revert one-indexing if the input is one-indexed
        atm_key_dct = {atm_key: atm_key - 1 for atm_key in atom_keys(gra)}
        gra = relabel(gra, atm_key_dct)

    return gra


# # conversions
def frozen(gra):
    """Generate a hashable, sortable, immutable representation of the graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A hashable, sortable, immutable representation of the graph
    :rtype: (tuple, tuple)
    """
    atm_keys = sorted(atom_keys(gra))
    bnd_keys = sorted(bond_keys(gra), key=sorted)

    # make it sortable by replacing Nones with -infinity
    atm_vals = numpy.array(dict_.values_by_key(atoms(gra), atm_keys), dtype=object)
    bnd_vals = numpy.array(dict_.values_by_key(bonds(gra), bnd_keys), dtype=object)
    atm_vals[numpy.equal(atm_vals, None)] = -numpy.inf
    bnd_vals[numpy.equal(bnd_vals, None)] = -numpy.inf

    frz_atms = tuple(zip(atm_keys, map(tuple, atm_vals)))
    frz_bnds = tuple(zip(bnd_keys, map(tuple, bnd_vals)))
    return (frz_atms, frz_bnds)


def formula(gra):
    """Generate a stoichiometric formula dictionary from a molecular graph.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :rtype: dict[str: int]
    """

    gra = explicit(gra)
    syms = atom_symbols(gra).values()
    fml = util.formula_from_symbols(syms)

    return fml


def symbols(gra) -> tuple:
    """Get the atomic symbols for this graph, sorted by atom key

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The symbols as a list
    :rtype: tuple
    """
    keys = sorted(atom_keys(gra))
    symbs = tuple(dict_.values_by_key(atom_symbols(gra), keys))
    return symbs


# # sorting
def argsort_by_size(gras, descending: bool = True) -> List[int]:
    """Get indices to sort graphs by size

    Currently sorts on three criteria, in order:
        1. Heavy atom count
        2. Total atom count
        3. Electron count

    :param gras: The list of molecule graphs to sort
    :param descending: Sort in descening order, largest to smallest?, defaults to True
    :type descending: bool, optional
    :return: The list indices in sorted order
    :rtype: List[int]
    """

    def _sort_value(arg):
        _, gra = arg
        return (
            atom_count(gra, heavy_only=True),  # 1. heavy atom count
            atom_count(gra, heavy_only=False),  # 2. total atom count
            electron_count(gra),  # 3. electron count
        )

    idxs, _ = zip(*sorted(enumerate(gras), key=_sort_value, reverse=descending))
    return idxs


def sort_by_size(gras, descending: bool = True):
    """Sort graphs by size

    Currently sorts on three criteria, in order:
        1. Heavy atom count
        2. Total atom count
        3. Electron count

    :param gras: The list of molecule graphs to sort
    :param descending: Sort in descening order, largest to smallest?, defaults to True
    :type descending: bool, optional
    :return: The sorted list of molecular grpahs
    """
    gras = tuple(gras)
    idxs = argsort_by_size(gras, descending=descending)
    return tuple(map(gras.__getitem__, idxs))


# # properties
def atom_count(
    gra, symb=None, heavy_only=False, dummy=False, with_implicit=True, keys=None
):
    """Count the number of atoms in the molecule

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
    if not dummy and symb != "X":
        gra = without_dummy_atoms(gra)

    # Get the keys
    keys = atom_keys(gra) if keys is None else keys
    symb_dct = atom_symbols(gra)
    if heavy_only:
        # If only including heavy atoms, filter non-heavy atoms out
        symb_dct = dict_.filter_by_value(symb_dct, lambda s: ptab.to_number(s) > 1)
    # Restrict count to the requested subset of atoms
    symbs = [symb_dct[k] for k in keys if k in symb_dct]
    natms = len(symbs) if symb is None else symbs.count(symb)

    if with_implicit and symb in ("H", None) and not heavy_only:
        atm_imp_hyd_dct = atom_implicit_hydrogens(gra)
        natms += sum(atm_imp_hyd_dct.values())

    return natms


def electron_count(gra, charge=0):
    """Count the number of electrons in the molecule, based on its charge

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


def atom_stereo_keys(gra, symb=None, excl_symbs=()):
    """Get the keys of atom stereo-centers this molecular graph

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
    keys = atom_keys(gra, symb=symb, excl_symbs=excl_symbs)
    par_dct = dict_.by_key(atom_stereo_parities(gra), keys)
    ste_keys = dict_.keys_by_value(par_dct, lambda x: x is not None)
    return frozenset(ste_keys)


def bond_stereo_keys(gra):
    """Get the keys of bond stereo-centers this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The bond keys
    :rtype: frozenset[{int, int}]
    """
    par_dct = bond_stereo_parities(gra)
    ste_keys = dict_.keys_by_value(par_dct, lambda x: x is not None)
    return frozenset(ste_keys)


def stereo_keys(gra):
    """Get the keys of atom and bond stereo-centers this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom and bond keys
    :rtype: frozenset
    """
    ste_keys = atom_stereo_keys(gra) | bond_stereo_keys(gra)
    return ste_keys


def has_atom_stereo(gra, symb=None, excl_symbs=()):
    """Does this graph have atom stereochemistry?

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symb: Optionally, restrict this to atoms with a particular
        atomic symbol (e.g., 'H' for hydrogens).
    :type symb: str
    :param excl_symbs: Optionally, exclude atoms with particular atomic
        symbols.
    :type excl_symbs: tuple[str]
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return bool(atom_stereo_keys(gra, symb=symb, excl_symbs=excl_symbs))


def has_bond_stereo(gra):
    """Does this graph have bond stereochemistry?

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return bool(bond_stereo_keys(gra))


def has_stereo(gra):
    """Does this graph have stereochemistry of any kind?

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return bool(stereo_keys(gra))


def has_pi_bonds(gra):
    """Does this graph have pi bonds?

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
    """Is this a TS graph?

    This is based on whether or not it has 0.1 or 0.9-order bonds.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: `True` if it is, `False` if it isn't
    :rtype: bool
    """
    bnd_ords = bond_orders(gra).values()
    return any(round(o, 1) in (0.1, 0.9) for o in bnd_ords)


def atomic_numbers(gra):
    """Get atomic numbers for atoms in this graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of atomic numbers, by atom key
    :rtype: dict[int: int]
    """
    symb_dct = atom_symbols(gra)
    anum_dct = dict_.transform_values(symb_dct, ptab.to_number)
    return anum_dct


def mass_numbers(gra):
    """Get mass numbers for atoms in this graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of mass numbers, by atom key
    :rtype: dict[int: float]
    """
    symb_dct = atom_symbols(gra)
    mnum_dct = dict_.transform_values(symb_dct, ptab.to_mass_number)
    return mnum_dct


def van_der_waals_radii(gra, angstrom=False):
    """Get van der Waals radii for atoms in this graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param angstrom: Return in angstroms, instead of bohr?, default False
    :type angstrom: bool, optional
    :returns: A dictionary of covalent, by atom key
    :rtype: dict[int: float]
    """
    symb_dct = atom_symbols(gra)
    rvdw_dct = dict_.transform_values(
        symb_dct, lambda s: ptab.van_der_waals_radius(s, angstrom=angstrom)
    )
    return rvdw_dct


def covalent_radii(gra, angstrom=False):
    """Get covalent radii for atoms in this graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param angstrom: Return in angstroms, instead of bohr?, default False
    :type angstrom: bool, optional
    :returns: A dictionary of covalent, by atom key
    :rtype: dict[int: float]
    """
    symb_dct = atom_symbols(gra)
    rcov_dct = dict_.transform_values(
        symb_dct, lambda s: ptab.covalent_radius(s, angstrom=angstrom)
    )
    return rcov_dct


def atomic_valences(gra):
    """Get atomic valences for atoms in this graph, as a dictionary

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
    """Get the number of lone pairs for atoms in this graph, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary of lone pair counts, by atom key
    :rtype: dict[int: int]
    """
    atm_symb_dct = atom_symbols(gra)
    atm_lpc_dct = dict_.transform_values(atm_symb_dct, ptab.lone_pair_count)
    return atm_lpc_dct


def lone_pair_atom_keys(gra):
    """Get the keys of atoms in this graph that have lone pairs

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    lpc_dct = atom_lone_pairs(gra)
    return frozenset(key for key, lpc in lpc_dct.items() if lpc > 0)


def atom_van_der_waals_radius(gra, key, angstrom=True):
    """Get the van der Waals radius of an atom in the graph

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


def atom_bond_counts(gra, bond_order=True, with_implicit=True, ts_=True):
    """Get the bond counts for the atoms in this graph, as a dictionary

    The bond count will be based on the graph's bond orders, so one must
    convert to a Kekule graph in order to include pi bonds in the count.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bond_order: Use the bond orders in the graph?, defaults to True
    :type bond_order: bool, optional
    :param with_implicit: Include implicit hydrogens?, defaults to True
    :type with_implicit: bool, optional
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: A dictionary of bond counts, by atom key
    :rtype: dict[int: int]
    """
    atm_keys = list(atom_keys(gra))
    gra = gra if ts_ else ts_reagents_graph_without_stereo(gra)
    # Convert to explicit graph if we watn implicit hydrogens in the count
    gra = explicit(gra) if with_implicit else gra

    if not bond_order:
        gra = without_pi_bonds(gra)

    atm_nbhs = dict_.values_by_key(atom_neighborhoods(gra), atm_keys)
    atm_nbnds = [sum(bond_orders(nbh).values()) for nbh in atm_nbhs]
    atm_nbnd_dct = dict_.transform_values(dict(zip(atm_keys, atm_nbnds)), int)
    return atm_nbnd_dct


def atom_unpaired_electrons(gra, bond_order=True):
    """The number of unpaired electrons on each atom, calculated as the atomic
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
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_keys = list(atom_keys(gra))
    if not bond_order:
        gra = without_pi_bonds(gra)

    atm_bnd_vlcs = dict_.values_by_key(atom_bond_counts(gra), atm_keys)
    atm_tot_vlcs = dict_.values_by_key(atomic_valences(gra), atm_keys)
    atm_rad_vlcs = numpy.subtract(atm_tot_vlcs, atm_bnd_vlcs)
    atm_unp_dct = dict_.transform_values(dict(zip(atm_keys, atm_rad_vlcs)), int)
    return atm_unp_dct


def bond_unpaired_electrons(gra, bond_order=True):
    """The number of adjacent pairs of unpaired electrons across each bond,
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
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    bnd_keys = list(bond_keys(gra))
    if not bond_order:
        gra = without_pi_bonds(gra)

    # determine unsaturated valences for each atom
    atm_unp_dct = atom_unpaired_electrons(gra)
    # a bond's unsaturated is set by its limiting atom
    bnd_unsats = [min(map(atm_unp_dct.__getitem__, k)) for k in bnd_keys]
    bnd_unsat_dct = dict(zip(bnd_keys, bnd_unsats))
    return bnd_unsat_dct


def tetrahedral_atom_keys(gra):
    """Get keys to stereo candidate atoms, which will be stereogenic if their
    groups are all distinct

    Atoms will be considered stereo candidates if they are tetrahderal, which
    happens when:
        a. They are bonded to 4 other atoms and have no lone pairs.
        b. They are bonded to 3 other atoms and have one lone pair.

    For TS graphs, atom keys which meet these criteria for *either* the
    reactants *or* the products will be included, to insure that all
    stereochemistry gets captured.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    if is_ts_graph(gra):
        gras = [
            ts_reagents_graph_without_stereo(gra, prod=False),
            ts_reagents_graph_without_stereo(gra, prod=True),
        ]
    else:
        gras = [gra]

    keys = sorted(atom_keys(gra))
    tet_atm_keys = []
    for gra_ in gras:
        nbnds_lst = dict_.values_by_key(atom_bond_counts(gra_), keys)
        nlps_lst = dict_.values_by_key(atom_lone_pairs(gra_), keys)
        for key, nbnds, nlps in zip(keys, nbnds_lst, nlps_lst):
            # Cases: 4 bonds 0 lone pairs OR 3 bonds and 1 lone pair
            if (nbnds == 4 and nlps == 0) or (nbnds == 3 and nlps == 1):
                tet_atm_keys.append(key)

    return frozenset(tet_atm_keys)


def maximum_spin_multiplicity(gra, bond_order=True):
    """The highest possible spin multiplicity for this molecular graph

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
    """Possible spin multiplicities for this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bond_order: Use the bond orders in the graph?, defaults to True
    :type bond_order: bool, optional
    :returns: The list of possible spin multiplicities
    :rtype: tuple
    """
    mult_max = maximum_spin_multiplicity(gra, bond_order=bond_order)
    mult_min = 2 if mult_max % 2 == 0 else 1
    mults = tuple(range(mult_min, mult_max + 1, 2))
    return mults


def atom_symbol_keys(gra):
    """Group the atom keys by atomic symbol

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
    """Get the backbone hydrogen keys of this molecular graph

    These are hydrogen keys which cannot be made implicit, because they are
    part of the backbone. There are two cases: (1.) one of the hydrogens in
    H2 must be considered a backbone hydrogen, and (2.) any multivalent
    hydrogen must be treated as part of the backbone.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    hyd_keys = atom_keys(gra, symb="H")
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    atm_imp_hyd_dct = atom_implicit_hydrogens(gra)
    atm_rxn_keys = ts_reacting_atom_keys(gra)

    def _is_backbone(hyd_key):
        is_h2 = all(
            ngb_key in hyd_keys and hyd_key < ngb_key
            for ngb_key in atm_ngb_keys_dct[hyd_key]
        )
        is_multivalent = len(atm_ngb_keys_dct[hyd_key]) > 1
        has_implicit_hydrogens = atm_imp_hyd_dct[hyd_key] > 0
        is_reacting = hyd_key in atm_rxn_keys
        return is_h2 or is_multivalent or has_implicit_hydrogens or is_reacting

    bbn_hyd_keys = frozenset(filter(_is_backbone, hyd_keys))
    return bbn_hyd_keys


def nonbackbone_hydrogen_keys(gra):
    """Get the hydrogen keys of this graph, with or without backbone hydrogens

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    return atom_keys(gra, symb="H") - backbone_hydrogen_keys(gra)


def atom_nonbackbone_hydrogen_keys(gra):
    """Get the hydrogen keys of each atom, with or without backbone hydrogens,
    as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary giving the hydrogen keys for each atom, by key.
    :rtypoe: dict[int: frozenset]
    """
    hyd_keys = nonbackbone_hydrogen_keys(gra)
    atm_hyd_keys_dct = dict_.transform_values(
        atoms_neighbor_atom_keys(gra), lambda x: x & hyd_keys
    )
    return atm_hyd_keys_dct


def backbone_keys(gra, hyd=True):
    """Get the backbone atom keys of this graph, including backbone hydrogens

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param hyd: Include backbone hydrogen keys?
    :type hyd: bool
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    bbn_keys = atom_keys(gra, excl_symbs=("H",))
    if hyd:
        bbn_keys |= backbone_hydrogen_keys(gra)
    return bbn_keys


def backbone_bond_keys(gra, terminal=False):
    """Get the backbone bond keys of this graph (bonds between backbone atoms)

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
    """Get backbone hydrogen keys for each atom, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A dictionary giving the hydrogen keys for each atom, by key
    :rtype: dict[int: frozenset]
    """
    hyd_keys = backbone_hydrogen_keys(gra)
    atm_hyd_keys_dct = dict_.transform_values(
        atoms_neighbor_atom_keys(gra), lambda x: x & hyd_keys
    )
    return atm_hyd_keys_dct


def terminal_atom_keys(gra, backbone=True):
    """Get the backbone atom keys of this graph, including backbone hydrogens

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
    """Get the keys of unsaturated (radical or pi-bonded) atoms

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The atom keys
    :rtype: frozenset[int]
    """
    atm_unp_dct = atom_unpaired_electrons(gra, bond_order=False)
    unsat_atm_keys = frozenset(dict_.keys_by_value(atm_unp_dct, bool))
    return unsat_atm_keys


def unsaturated_bond_keys(gra):
    """Get the keys of unsaturated (radical or pi-bonded) bonds

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The bond keys
    :rtype: frozenset[frozenset[int]]
    """
    bnd_unp_dct = bond_unpaired_electrons(gra, bond_order=False)
    unsat_bnd_keys = frozenset(dict_.keys_by_value(bnd_unp_dct, bool))
    return unsat_bnd_keys


def angle_keys(gra):
    """Get triples of keys for pairs of adjacent bonds, with the central atom
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
            (atm2_key,) = bnd_key1 & bnd_key2
            (atm1_key,) = bnd_key1 - {atm2_key}
            (atm3_key,) = bnd_key2 - {atm2_key}
            ang_keys.append((atm1_key, atm2_key, atm3_key))
            ang_keys.append((atm3_key, atm2_key, atm1_key))

    return frozenset(ang_keys)


# # relabeling and changing keys
def relabel(gra, atm_key_dct, check=True):
    """Relabel the atoms in the graph with new keys, using a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_key_dct: New keys for a subset of the atoms, by current atom key
    :type atm_key_dct: dict[int: int]
    :param check: Check that all keys in `atm_key_dct` are present in the graph?
    :type check: bool
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    orig_atm_keys = atom_keys(gra)
    if check:
        assert (
            set(atm_key_dct.keys()) <= orig_atm_keys
        ), f"{set(atm_key_dct.keys())}\n{orig_atm_keys}"

    new_atm_key_dct = dict(zip(orig_atm_keys, orig_atm_keys))
    new_atm_key_dct.update(atm_key_dct)

    _relabel_atom_key = new_atm_key_dct.__getitem__

    def _relabel_bond_key(bnd_key):
        return frozenset(map(_relabel_atom_key, bnd_key))

    atm_dct = dict_.transform_keys(atoms(gra), _relabel_atom_key)
    bnd_dct = dict_.transform_keys(bonds(gra), _relabel_bond_key)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def standard_keys(gra):
    """Relabel the atoms in the graph with standard zero-indexed keys

    The new keys will follow the same sort order as the original ones.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_key_dct = dict(map(reversed, enumerate(sorted(atom_keys(gra)))))
    return relabel(gra, atm_key_dct)


def standard_keys_for_sequence(gras):
    """Give standard, non-overlapping keys to a sequence of graphs and return
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

        atm_key_dct = {
            atm_key: idx + shift for idx, atm_key in enumerate(sorted(atom_keys(gra)))
        }
        atm_key_dcts.append(atm_key_dct)

        shift += natms

    gras = tuple(
        relabel(gra, atm_key_dct) for gra, atm_key_dct in zip(gras, atm_key_dcts)
    )
    atm_key_dcts = tuple(atm_key_dcts)

    return gras, atm_key_dcts


def apply_zmatrix_conversion(gra, zc_: ZmatConv):
    """Apply a z-matrix conversion to the graph

    This can be used to match a z-matrix after geometry -> z-matrix conversion.

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param zc_: A z-matrix conversion
    :type zc_: ZmatConv
    :returns: The transformed graph
    :rtype: automol graph data structure
    """
    gra = relabel(gra, zmat_conv.relabel_dict(zc_))
    for dummy_key, parent_key in zmat_conv.insert_dict(zc_).items():
        gra = add_bonded_atom(gra, "X", parent_key, dummy_key, bnd_ord=0)
    return gra


def reverse_zmatrix_conversion(gra, zc_: ZmatConv):
    """Reverse a z-matrix conversion, recovering the original graph

    This can be used to match the original geometry after geometry -> z-matrix
    conversion.

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param zc_: A z-matrix conversion
    :type zc_: ZmatConv
    :returns: The transformed graph
    :rtype: automol graph data structure
    """
    gra = remove_atoms(gra, zmat_conv.insert_dict(zc_))
    gra = relabel(gra, zmat_conv.relabel_dict(zc_, rev=True))
    return gra


# # add/remove/insert/without
def add_atoms(gra, symb_dct, imp_hyd_dct=None, ste_par_dct=None, check=True):
    """Add atoms to this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symb_dct: atomic symbols, by atom key
    :type symb_dct: dict
    :param imp_hyd_dct: the number of implicit hydrogens associated with
        each atom, by atom key
    :type imp_hyd_dct: dict
    :param ste_par_dct: stereo parities, by atom key;
    :type ste_par_dct: dict
    :param check: Check that we aren't trying to add atoms with duplicate keys?
    :type check: bool
    :return: molecular graph (TS or non-TS)
    :rtype: automol graph data structure
    """
    atm_keys = atom_keys(gra)
    atm_symb_dct = atom_symbols(gra)
    atm_imp_hyd_dct = atom_implicit_hydrogens(gra)
    atm_ste_par_dct = atom_stereo_parities(gra)

    keys = set(symb_dct.keys())
    imp_hyd_dct = {} if imp_hyd_dct is None else imp_hyd_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct

    if check:
        assert (
            not keys & atm_keys
        ), f"{keys} and {atm_keys} have a non-empty intersection"

    assert not keys & atm_keys
    assert set(imp_hyd_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys

    atm_symb_dct.update(symb_dct)
    atm_imp_hyd_dct.update(imp_hyd_dct)
    atm_ste_par_dct.update(ste_par_dct)

    atm_dct = atoms_from_data(
        atm_symb_dct=atm_symb_dct,
        atm_imp_hyd_dct=atm_imp_hyd_dct,
        atm_ste_par_dct=atm_ste_par_dct,
    )
    bnd_dct = bonds(gra)
    gra = from_atoms_and_bonds(atm_dct=atm_dct, bnd_dct=bnd_dct)
    return gra


def add_bonds(gra, keys, ord_dct=None, ste_par_dct=None, check=True):
    """Add bonds to this molecular graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param keys: bond keys
    :type keys: set
    :param ord_dct: bond orders, by bond key
    :type ord_dct: dict
    :param ste_par_dct: stereo parities, by bond key;
    :type ste_par_dct: dict
    :param check: Check that we aren't trying to add bonds with duplicate keys?
    :type check: bool
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    bnd_keys = set(bond_keys(gra))
    bnd_ord_dct = bond_orders(gra)
    bnd_ste_par_dct = bond_stereo_parities(gra)

    keys = set(map(frozenset, keys))
    ord_dct = {} if ord_dct is None else ord_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct

    ord_dct = dict_.transform_keys(ord_dct, frozenset)
    ste_par_dct = dict_.transform_keys(ste_par_dct, frozenset)

    if check:
        assert (
            not keys & bnd_keys
        ), f"{keys} and {bnd_keys} have a non-empty intersection"

    assert set(ord_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys

    bnd_keys.update(keys)
    bnd_ord_dct.update(ord_dct)
    bnd_ste_par_dct.update(ste_par_dct)

    atm_dct = atoms(gra)
    bnd_dct = bonds_from_data(
        bnd_keys=bnd_keys, bnd_ord_dct=bnd_ord_dct, bnd_ste_par_dct=bnd_ste_par_dct
    )

    gra = from_atoms_and_bonds(atm_dct=atm_dct, bnd_dct=bnd_dct)
    return gra


def remove_atoms(gra, atm_keys, check=True, stereo=True):
    """Remove atoms from this molecular graph

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
    """Remove bonds from this molecular graph

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
    """Change the implicit hydrogen count for atoms in this graph

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
        dict_.values_by_key(imp_hyd_change_dct, atm_keys),
    )
    assert all(atm_imp_hyd >= 0 for atm_imp_hyd in atm_imp_hyds)
    atm_imp_hyd_dct = dict_.transform_values(dict(zip(atm_keys, atm_imp_hyds)), int)
    return set_atom_implicit_hydrogens(gra, atm_imp_hyd_dct)


def add_atom_explicit_hydrogens(gra, exp_hyd_keys_dct):
    """Add explicit hydrogens by atom

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
        f"{set(exp_hyd_keys_dct.keys())}" " !<= " f"{atom_keys(gra)}"
    )
    for atm_key, atm_exp_hyd_keys in exp_hyd_keys_dct.items():
        assert not set(atm_exp_hyd_keys) & atom_keys(gra)
        atm_exp_hyd_bnd_keys = {
            frozenset({atm_key, atm_exp_hyd_key})
            for atm_exp_hyd_key in atm_exp_hyd_keys
        }
        atm_exp_hyd_symb_dct = dict_.by_key({}, atm_exp_hyd_keys, fill_val="H")
        gra = add_atoms(gra, atm_exp_hyd_symb_dct)
        gra = add_bonds(gra, atm_exp_hyd_bnd_keys)
    return gra


def add_bonded_atom(
    gra,
    symb,
    atm_key,
    bnd_atm_key=None,
    imp_hyd=None,
    atm_ste_par=None,
    bnd_ord=None,
    bnd_ste_par=None,
):
    """Add a single atom and connect it to an atom already in the graph

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
        f"Cannot add atom with key {bnd_atm_key}.\n"
        f"It is already in the graph:\n{gra}"
    )

    symb_dct = {bnd_atm_key: symb}
    imp_hyd_dct = {bnd_atm_key: imp_hyd} if imp_hyd is not None else None
    atm_ste_par_dct = {bnd_atm_key: atm_ste_par} if atm_ste_par is not None else None

    gra = add_atoms(gra, symb_dct, imp_hyd_dct=imp_hyd_dct, ste_par_dct=atm_ste_par_dct)

    bnd_key = frozenset({bnd_atm_key, atm_key})
    bnd_ord_dct = {bnd_key: bnd_ord} if bnd_ord is not None else None
    bnd_ste_par_dct = {bnd_key: bnd_ste_par} if bnd_ste_par is not None else None

    gra = add_bonds(gra, [bnd_key], ord_dct=bnd_ord_dct, ste_par_dct=bnd_ste_par_dct)

    return gra


def without_pi_bonds(gra):
    """Get a version of this graph without any pi-bonds

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
    bnd_ords = [
        0
        if round(v, 1) == 0
        else round(v % 1, 1)
        if round(v % 1, 1) in (0.1, 0.9)
        else 1
        for v in map(bnd_ord_dct.__getitem__, bnd_keys)
    ]
    bnd_ord_dct = dict(zip(bnd_keys, bnd_ords))
    return set_bond_orders(gra, bnd_ord_dct)


def without_reacting_bonds(gra):
    """Get a copy of this graph without reacting bonds

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    return without_bonds_by_orders(gra, ords=(0.1, 0.9))


def without_bonds_by_orders(gra, ords=(0,), skip_dummies=True):
    """Remove bonds of certain orders from the graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ords: The orders (or order) of bonds to be removed
    :type ords: list or int
    :param skip_dummies: Skip bonds to dummy atoms in the removal?
    :type skip_dummies: bool
    """
    ords = [ords] if isinstance(ords, numbers.Number) else ords
    ord_dct = dict_.filter_by_value(
        bond_orders(gra), func=lambda x: round(x, 1) in ords
    )
    rem_keys = ord_dct.keys()
    if skip_dummies:
        xkeys = atom_keys(gra, symb="X")
        rem_keys = [bk for bk in rem_keys if not any(k in xkeys for k in bk)]
    gra = remove_bonds(gra, rem_keys)
    return gra


def without_dummy_atoms(gra):
    """Remove dummy atoms from this graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    symb_dct = atom_symbols(gra)
    keys = [key for key, sym in symb_dct.items() if ptab.to_number(sym)]
    return subgraph(gra, keys, stereo=True)


def without_stereo(gra, atm_keys=None, bnd_keys=None):
    """Remove stereo information (atom and bond parities) from this graph

    For TS graphs, product and TS stereochemistry will be removed as well

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
    # ^^^ If one of these is set, assume we want to wipe out *only* those atom
    # or bond keys (not, e.g., a selecton of bonds and all of the atoms)

    atm_keys = atom_keys(gra) if atm_keys is None else atm_keys
    bnd_keys = bond_keys(gra) if bnd_keys is None else bnd_keys
    atm_ste_par_dct = dict_.by_key({}, atm_keys, fill_val=None)
    bnd_ste_par_dct = dict_.by_key({}, bnd_keys, fill_val=None)
    gra = set_atom_stereo_parities(gra, atm_ste_par_dct)
    gra = set_bond_stereo_parities(gra, bnd_ste_par_dct)
    return gra


def explicit(gra, atm_keys=None):
    """Make implicit hydrogens in this graph explicit

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_keys: Optionally, restrict this operation to a subset of atoms
    :type atm_keys: list[int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    if not atoms(gra):
        return gra

    atm_keys = backbone_keys(gra) if atm_keys is None else atm_keys
    atm_keys = sorted(atm_keys)
    atm_imp_hyd_dct = dict_.by_key(atom_implicit_hydrogens(gra), atm_keys)

    atm_exp_hyd_keys_dct = {}
    next_atm_key = max(atom_keys(gra)) + 1
    for atm_key in atm_keys:
        imp_hyd = atm_imp_hyd_dct[atm_key]
        atm_exp_hyd_keys_dct[atm_key] = set(range(next_atm_key, next_atm_key + imp_hyd))
        next_atm_key += imp_hyd

    gra = set_atom_implicit_hydrogens(gra, dict_.by_key({}, atm_keys, fill_val=0))
    gra = add_atom_explicit_hydrogens(gra, atm_exp_hyd_keys_dct)
    return gra


def implicit(gra, atm_keys=None):
    """Make explicit hydrogens in this graph implicit

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atm_keys: Optionally, restrict this operation to a subset of atoms
    :type atm_keys: list[int]
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    atm_keys = backbone_keys(gra) if atm_keys is None else atm_keys

    atm_exp_hyd_keys_dct = dict_.by_key(atom_nonbackbone_hydrogen_keys(gra), atm_keys)

    inc_imp_hyd_keys_dct = dict_.transform_values(atm_exp_hyd_keys_dct, len)
    gra = change_implicit_hydrogens(gra, inc_imp_hyd_keys_dct)

    exp_hyd_keys = set(itertools.chain(*atm_exp_hyd_keys_dct.values()))
    gra = remove_atoms(gra, exp_hyd_keys, stereo=True)
    return gra


# # unions
def union(gra1, gra2, check=True):
    """Get the union of two molecular graphs

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
    """Get the union of a sequence of graphs

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
    """Get the subgraph of a subset of the atoms in this graph

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
    """Get the subgraph of a subset of bonds in this graph

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


def atom_neighborhood(
    gra, atm_key, bnd_keys=None, stereo=False, ts_=True, second_degree: bool = False
):
    """Get the neighborhood subgraph of a specific atom

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
    :param second_degree: Include second-degree neighbors?
    :type second_degree: bool, optional
    :returns: A molecular graph
    :rtype: automol graph data structure
    """
    gra = gra if ts_ else ts_reagents_graph_without_stereo(gra)
    bnd_keys = bond_keys(gra, ts_=ts_) if bnd_keys is None else bnd_keys
    nbh_bnd_keys = set(k for k in bnd_keys if atm_key in k)

    # If we are allowing second-degree neighbors, extend this to include them
    if second_degree:
        nbh_atm_keys = set(itertools.chain(*nbh_bnd_keys))
        nbh_bnd_keys = set(k for k in bnd_keys if nbh_atm_keys & k)

    if nbh_bnd_keys:
        nbh = bond_induced_subgraph(gra, nbh_bnd_keys, stereo=stereo)
    else:
        nbh = subgraph(gra, {atm_key}, stereo=stereo)

    return nbh


def atom_neighborhoods(gra, bnd_keys=None, stereo=False, ts_=True):
    """Get the neighborhood subgraphs of each atom, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param bnd_keys: optionally, restrict this to a subset of the bond keys
    :type bnd_keys: tuple[frozenset[int]]
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: neighborhood subgraphs, by atom key
    :rtype: dict[int: automol graph data structure]
    """
    bnd_keys = bond_keys(gra, ts_=ts_) if bnd_keys is None else bnd_keys

    def _neighborhood(atm_key):
        return atom_neighborhood(
            gra, atm_key, bnd_keys=bnd_keys, stereo=stereo, ts_=ts_
        )

    atm_keys = list(atom_keys(gra))
    atm_nbh_dct = dict(zip(atm_keys, map(_neighborhood, atm_keys)))
    return atm_nbh_dct


def bond_neighborhood(gra, bnd_key, bnd_keys=None, stereo=False, ts_=True):
    """Get the neighborhood subgraph of a specific bond

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
    bnd_keys = bond_keys(gra, ts_=ts_) if bnd_keys is None else bnd_keys
    nbh_bnd_keys = set(filter(lambda x: bnd_key & x, bnd_keys))
    nbh = bond_induced_subgraph(gra, nbh_bnd_keys, stereo=stereo)
    return nbh


def bond_neighborhoods(gra, bnd_keys=None, stereo=False, ts_=True):
    """Get the neighborhood subgraphs of each bond, as a dictionary

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
        return bond_neighborhood(
            gra, bnd_key, bnd_keys=bnd_keys, stereo=stereo, ts_=ts_
        )

    bnd_nbh_dct = dict(zip(bnd_keys, map(_neighborhood, bnd_keys)))
    return bnd_nbh_dct


def atom_neighbor_atom_keys(
    gra,
    atm_key,
    bnd_keys=None,
    symb=None,
    excl_symbs=(),
    ts_=True,
    second_degree: bool = False,
    include_self: bool = False,
):
    """Get keys for an atom's neighbors

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
    :param second_degree: Include second-degree neighbors? defaults to False
    :type second_degree: bool, optional
    :param include_self: Include the atom itself in the returned keys? defaults to False
    :type include_self: bool, optional
    :returns: The keys of neighboring atoms
    :rtype: frozenset[int]
    """
    atm_nbh = atom_neighborhood(
        gra, atm_key, bnd_keys=bnd_keys, ts_=ts_, second_degree=second_degree
    )
    atm_nkeys = atom_keys(atm_nbh, symb=symb, excl_symbs=excl_symbs)

    if not include_self:
        atm_nkeys -= {atm_key}

    return frozenset(atm_nkeys)


def local_stereo_priorities(gra) -> Dict[int, int]:
    """Generate a local priority dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The local priorities, by atom key
    :rtype: Dict[int, int]
    """
    loc_pri_dct = {}
    loc_pri_dct.update({k: k for k in backbone_keys(gra, hyd=False)})
    loc_pri_dct.update({k: -abs(k) for k in backbone_hydrogen_keys(gra)})
    loc_pri_dct.update({k: -999 for k in nonbackbone_hydrogen_keys(gra)})
    return loc_pri_dct


def atom_stereo_sorted_neighbor_keys(gra, key, self_apex=False, pri_dct=None):
    """Get keys for the neighbors of an atom that are relevant for atom
    stereochemistry, sorted by priority (if requested)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param key: the atom key
    :type key: int
    :param self_apex: If there are only 3 neighbors, put this atom as the apex?
    :type self_apex: bool, optional
    :param pri_dct: Priorities to sort by (optional)
    :type pri_dct: dict[int: int]
    :returns: The keys of neighboring atoms
    :rtype: tuple[int]
    """
    gra = without_dummy_atoms(gra)
    pri_dct = local_stereo_priorities(gra) if pri_dct is None else pri_dct

    nkeys = atom_neighbor_atom_keys(gra, key, ts_=True)
    nkeys_no_brk = atom_neighbor_atom_keys(gra, key, ts_=False)
    nlps = atom_lone_pairs(gra)[key]
    # For Sn2 reactions, don't include the breaking bond key
    valence = len(nkeys) + nlps
    if valence > 4:
        assert valence == 5 and len(nkeys_no_brk) == len(nkeys) - 1, (
            f"Unanticipated valence {valence} at key {key} is not resoved by "
            f"dropping breaking bonds:\n{gra}"
        )
        nkeys = nkeys_no_brk

    # Sort them by priority
    nkeys = sorted(nkeys, key=pri_dct.__getitem__)

    # If there are only three groups, use the stereo atom itself as
    # the top apex of the tetrahedron.
    if self_apex and len(nkeys) < 4:
        assert len(nkeys) == 3
        nkeys = [key] + list(nkeys)

    return tuple(nkeys)


def bond_stereo_sorted_neighbor_keys(gra, key1, key2, pri_dct=None):
    """Get keys for the neighbors of a bond that are relevant for bond
    stereochemistry, sorted by priority (if requested)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param key1: the first atom in the bond
    :type key1: int
    :param key2: the second atom in the bond
    :type key2: int
    :param pri_dct: Priorities to sort by (optional)
    :type pri_dct: dict[int: int]
    :returns: The keys of neighboring atoms for the first and second atoms
    :rtype: tuple[int], tuple[int]
    """
    gra = without_dummy_atoms(gra)
    pri_dct = local_stereo_priorities(gra) if pri_dct is None else pri_dct

    keys = {key1, key2}
    gra_no_rxbs = without_bonds_by_orders(gra, [0.1, 0.9])

    def _neighbor_keys(key):
        nkeys = atom_neighbor_atom_keys(gra, key) - keys
        nkeys_no_rbs = atom_neighbor_atom_keys(gra_no_rxbs, key) - keys
        nlps = atom_lone_pairs(gra)[key]
        # Deal with cases of the form VC(W)=C(X)Y + Z <=> V[C](W)C(X)(Y)Z
        valence = len(nkeys) + nlps
        if valence > 2:
            assert valence == 3 and len(nkeys_no_rbs) == len(nkeys) - 1, (
                f"Unanticipated valence {valence} at key {key} is not resoved "
                f"by dropping breaking bonds:\n{gra}"
            )
            nkeys = nkeys_no_rbs
        # Sort them by priority
        nkeys = sorted(nkeys, key=pri_dct.__getitem__)
        return tuple(nkeys)

    return (_neighbor_keys(key1), _neighbor_keys(key2))


def atoms_neighbor_atom_keys(gra, ts_=True):
    """Get the keys for each atom's neighboring atoms, as a dictionary

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
        atom_neighborhoods(gra, ts_=ts_), _neighbor_keys
    )
    return atm_ngb_keys_dct


def atoms_sorted_neighbor_atom_keys(
    gra,
    symbs_first=("C",),
    symbs_last=("H",),
    ords_last=(0.1,),
    prioritize_keys=(),
    ts_=True,
):
    """Get keys for each atom's neighbors, sorted in a particular order

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
        ords = [-1 if o not in ords_last else ords_last.index(o) for o in ords]
        symbs = list(map(atm_symb_dct.__getitem__, keys))
        pris = [0 if k in prioritize_keys else 1 for k in keys]
        srt_vals = list(zip(ords, pris, symbs))
        srt = automol.formula.argsort_symbols(srt_vals, symbs_first, symbs_last, idx=2)
        keys = tuple(map(keys.__getitem__, srt))
        return keys

    atm_ngb_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(gra, ts_=ts_), _neighbor_keys
    )
    return atm_ngb_keys_dct


def atom_sorted_neighbor_atom_keys(
    gra,
    atm_key,
    excl_atm_keys=(),
    incl_atm_keys=None,
    symbs_first=("C",),
    symbs_last=("H",),
    ts_=True,
):
    """Get keys for this atom's neighbors, sorted in a particular order

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


def atom_neighbor_atom_key(
    gra,
    atm_key,
    excl_atm_keys=(),
    incl_atm_keys=None,
    symbs_first=("C",),
    symbs_last=("H",),
):
    """Get a key for one of an atom's neighbors

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
        gra,
        atm_key,
        excl_atm_keys=excl_atm_keys,
        incl_atm_keys=incl_atm_keys,
        symbs_first=symbs_first,
        symbs_last=symbs_last,
    )
    return atm_keys[0] if atm_keys else None


def atom_bond_keys(gra, atm_key, ts_=True):
    """Get the bond keys of a specific atom

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
    """Get the bond keys of each atom, as a dictionary

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: The bond keys of each atom, as a dictionary
    :rtype: dict[int: frozenset]
    """
    atm_nbhs = atom_neighborhoods(gra, ts_=ts_)
    return dict_.transform_values(atm_nbhs, bond_keys)


def dummy_parent_dict(gra, ts_=True):
    """Get the atoms that are connected to dummy atoms, by dummy atom key
    (Requires that each dummy atom only be connected to one neighbor)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param ts_: If this is a TS graph, treat it as such?
    :type ts_: bool
    :returns: The atoms that are connected to dummy atoms, by dummy atom key
    :rtype: dict[int: int]
    """
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra, ts_=ts_)
    dummy_atm_keys = atom_keys(gra, symb="X")

    dummy_parent_dct = {}
    for key in dummy_atm_keys:
        ngb_keys = atm_ngb_keys_dct[key]
        assert len(ngb_keys) == 1, "Dummy atoms should only be connected to one atom!"
        (ngb_key,) = ngb_keys
        dummy_parent_dct[key] = ngb_key

    return dummy_parent_dct


def bonds_neighbor_atom_keys(gra, ts_=True):
    """Get the keys of each bond's neighboring atoms, as a dictionary

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
        bond_neighborhoods(gra, ts_=ts_), _neighbor_keys
    )
    return bnd_ngb_keys_dct


def bonds_neighbor_bond_keys(gra, ts_=True):
    """Get the keys of each bond's neighboring bonds, as a dictionary

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
        bond_neighborhoods(gra, ts_=ts_), _neighbor_keys
    )
    return bnd_ngb_keys_dct
