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
# TS-specific functions in _core
from automol.graph.base._core import ts_forming_bond_keys
from automol.graph.base._core import ts_breaking_bond_keys
from automol.graph.base._core import ts_reacting_bonds
from automol.graph.base._core import ts_reacting_atoms
from automol.graph.base._core import ts_atom_product_stereo_parities
from automol.graph.base._core import ts_atom_fleeting_stereo_parities
from automol.graph.base._core import ts_bond_product_stereo_parities
from automol.graph.base._core import ts_bond_fleeting_stereo_parities
from automol.graph.base._core import ts_reactants_graph
from automol.graph.base._core import ts_set_atom_product_stereo_parities
from automol.graph.base._core import ts_set_atom_fleeting_stereo_parities
from automol.graph.base._core import ts_set_bond_product_stereo_parities
from automol.graph.base._core import ts_set_bond_fleeting_stereo_parities
# Non-TS functions from _core
from automol.graph.base._core import from_data
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import bond_keys
from automol.graph.base._core import atom_implicit_hydrogens
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_orders
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import add_bonds
from automol.graph.base._core import is_ts_graph

# Rename TS functions from core
forming_bond_keys = ts_forming_bond_keys
breaking_bond_keys = ts_breaking_bond_keys
reacting_bonds = ts_reacting_bonds
reacting_atoms = ts_reacting_atoms
atom_product_stereo_parities = ts_atom_product_stereo_parities
atom_fleeting_stereo_parities = ts_atom_fleeting_stereo_parities
bond_product_stereo_parities = ts_bond_product_stereo_parities
bond_fleeting_stereo_parities = ts_bond_fleeting_stereo_parities
reactants_graph = ts_reactants_graph
set_atom_product_stereo_parities = ts_set_atom_product_stereo_parities
set_atom_fleeting_stereo_parities = ts_set_atom_fleeting_stereo_parities
set_bond_product_stereo_parities = ts_set_bond_product_stereo_parities
set_bond_fleeting_stereo_parities = ts_set_bond_fleeting_stereo_parities


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
        bnd_ord_dct=bnd_ord_dct, bnd_ste_par_dct=bond_stereo_parities(gra),
        bnd_prd_ste_par_dct=bnd_prd_ste_par_dct,
        bnd_ts_ste_par_dct=bnd_ts_ste_par_dct,
        ts_=True
    )

    # Encode forming and breaking bonds
    frm_bnd_keys = frozenset(map(frozenset, frm_bnd_keys))
    brk_bnd_keys = frozenset(map(frozenset, brk_bnd_keys))

    frm_ord_dct = {k: 0.1 for k in frm_bnd_keys}
    brk_ord_dct = {k: bnd_ord_dct[k] - 0.1 if k in bnd_ord_dct else 0.9
                   for k in brk_bnd_keys}

    tsg = add_bonds(tsg, frm_bnd_keys, ord_dct=frm_ord_dct, check=False)
    tsg = add_bonds(tsg, brk_bnd_keys, ord_dct=brk_ord_dct, check=False)
    return tsg
