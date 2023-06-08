""" transition state graph data structure

Data structure:
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

Differences from the ordinary molecular graph data structure:

1. Extra atom and bond properties for marking product and TS stereo parities.
   (The first stereo parity property contains reactant stereo parities.)
2. Forming bonds are encoded as 0.1-order bonds and breaking bonds are encoded
   as 0.9-order bonds (for single bonds).
3. For resonance graphs, breaking double bonds will be encoded as 1.9 and
   breaking triple bonds as 2.9. Note, however, that "breaking" still refers to
   a total separation of the two atoms, not simply a bond order decrement.
   (Bond order decrements are not and should not be considered.)

However, ordinary graph functions should still work for TS graphs as well.

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
# # Non-TS functions from _core
# from automol.graph.base._core import from_data
# from automol.graph.base._core import atom_symbols
# from automol.graph.base._core import bond_keys
# from automol.graph.base._core import atom_implicit_hydrogens
# from automol.graph.base._core import atom_stereo_parities
# from automol.graph.base._core import bond_orders
# from automol.graph.base._core import bond_stereo_parities
# from automol.graph.base._core import add_bonds
# from automol.graph.base._core import is_ts_graph
from automol.util import dict_
from automol.graph.base._core import bond_orders
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import set_bond_orders
from automol.graph.base._core import set_atom_stereo_parities
from automol.graph.base._core import set_bond_stereo_parities
from automol.graph.base._core import has_pi_bonds
from automol.graph.base._core import string
# TS-specific functions in _core
from automol.graph.base._core import ts_graph
from automol.graph.base._core import ts_forming_bond_keys
from automol.graph.base._core import ts_breaking_bond_keys
from automol.graph.base._core import ts_reacting_bonds
from automol.graph.base._core import ts_reacting_atoms
from automol.graph.base._core import ts_reactants_graph
from automol.graph.base._core import ts_products_graph

# Rename TS functions from core
graph = ts_graph
forming_bond_keys = ts_forming_bond_keys
breaking_bond_keys = ts_breaking_bond_keys
reacting_bonds = ts_reacting_bonds
reacting_atoms = ts_reacting_atoms
reactants_graph = ts_reactants_graph
products_graph = ts_products_graph


def reverse(tsg):
    """ Reverse the reaction direction of a TS graph

    This is done by swapping (a.) breaking and forming bonds, and (b.) product
    and reactant stereo parities with each other, and by (c.) transforming any
    TS stereo parities as needed.

    :param tsg: TS graph
    :type tsg: automol TS graph data structure
    :returns: A TS graph for the reverse reaction
    :rtype: automol TS graph data structure
    """
    # Step 1: Swap bond orders for breaking and forming bonds
    assert not has_pi_bonds(tsg), (
        f"Cannot reverse TS graph with pi bonds:\n{string(tsg)}")
    bnd_ord_dct = {k: (0.9 if round(o, 1) == 0.1
                       else 0.1 if round(o, 1) == 0.9
                       else o) for k, o in bond_orders(tsg).items()}
    tsg = set_bond_orders(tsg, bnd_ord_dct)

    # Step 2: Swap reactant and product stereo parities
    r_atm_ste_par_dct = atom_stereo_parities(tsg, ts_select='R')
    p_atm_ste_par_dct = atom_stereo_parities(tsg, ts_select='P')
    r_bnd_ste_par_dct = bond_stereo_parities(tsg, ts_select='R')
    p_bnd_ste_par_dct = bond_stereo_parities(tsg, ts_select='P')
    tsg = set_atom_stereo_parities(tsg, p_atm_ste_par_dct, ts_select='R')
    tsg = set_atom_stereo_parities(tsg, r_atm_ste_par_dct, ts_select='P')
    tsg = set_bond_stereo_parities(tsg, p_bnd_ste_par_dct, ts_select='R')
    tsg = set_bond_stereo_parities(tsg, r_bnd_ste_par_dct, ts_select='P')

    # Step 3: Invert TS stereo where necessary
    atm_ts_ste_par_dct = dict_.filter_by_value(
        atom_stereo_parities(tsg, ts_select='T'), lambda x: x is not None)
    bnd_ts_ste_par_dct = dict_.filter_by_value(
        bond_stereo_parities(tsg, ts_select='T'), lambda x: x is not None)
    if atm_ts_ste_par_dct or bnd_ts_ste_par_dct:
        raise NotImplementedError("Not yet implemented for TS parities.")

    return tsg
