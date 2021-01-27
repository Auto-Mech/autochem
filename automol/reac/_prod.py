"""
  Build graphs for the products of reactions
"""

import itertools
import automol.graph
from automol.graph._graph_base import atom_symbol_idxs
from automol.graph._graph import atom_keys
from automol.graph._graph import union
from automol.graph._graph import full_isomorphism
from automol.graph._graph import add_bonds
from automol.graph._graph import remove_atoms
from automol.graph._graph import remove_bonds
from automol.graph._graph import unsaturated_atom_keys
from automol.graph._graph import add_atom_explicit_hydrogen_keys
from automol.graph._graph import atoms_neighbor_atom_keys
from automol.graph._res import resonance_dominant_radical_atom_keys
from automol.graph._func_group import chem_unique_atoms_of_type
from automol.graph._func_group import bonds_of_order


def prod_hydrogen_abstraction(x_gra, y_gra):
    """ products of hydrogen loss (generalize to group loss?)
    """

    prod_gras = tuple()

    # Build the H atom graph
    num_y_gra = len(automol.graph.atoms(y_gra)) + 1
    h_gra = ({num_y_gra: ('H', 0, None)}, {})

    # Form graphs where rct1 loses H and rct2 gains H
    x_gras_hloss = _prod_group_loss(x_gra, grp='H')
    y_gras_hadd = prod_addition(y_gra, h_gra)

    for gra1 in x_gras_hloss:
        for gra2 in y_gras_hadd:
            prod_gras += ((gra1, gra2),)

    return prod_gras


def prod_addition(x_gra, y_gra):
    """ products of addition
    """

    prod_gras = tuple()

    shift = len(automol.graph.atoms(x_gra))
    y_gra = automol.graph.transform_keys(y_gra, lambda x: x+shift)

    x_keys = unsaturated_atom_keys(x_gra)
    y_keys = unsaturated_atom_keys(y_gra)

    for x_key, y_key in itertools.product(x_keys, y_keys):
        xy_gra = add_bonds(union(x_gra, y_gra), [{x_key, y_key}])
        prod_gras += ((xy_gra,),)

    return _unique_gras(prod_gras)


def prod_hydrogen_migration(gra):
    """ products of hydrogen migration
    """

    prod_gras = tuple()

    keys = atom_keys(gra)

    num_keys = len(keys)
    if num_keys > 2:
        rad_idxs = resonance_dominant_radical_atom_keys(gra)
        uni_h_idxs = chem_unique_atoms_of_type(gra, 'H')

        h_atm_key = max(keys) + 1

        for h_idx in uni_h_idxs:
            for rad_idx in rad_idxs:
                gra2 = remove_atoms(gra, [h_idx])
                gra2_h = add_atom_explicit_hydrogen_keys(
                    gra2, {rad_idx: [h_atm_key]})
                if not full_isomorphism(gra, gra2_h):
                    prod_gras += ((gra2_h,),)

    return _unique_gras(prod_gras)


def prod_beta_scission(gra):
    """ products of beta scission
    """

    prod_gras = tuple()

    rad_idxs = resonance_dominant_radical_atom_keys(gra)
    single_bonds = bonds_of_order(gra, mbond=1)

    for rad_idx in rad_idxs:
        rad_neighs = atoms_neighbor_atom_keys(gra)[rad_idx]
        for single_bond in single_bonds:
            bond = frozenset(single_bond)
            if rad_neighs & bond and rad_idx not in bond:
                gra2 = remove_bonds(gra, [bond])
                disconn_gras = automol.graph.connected_components(gra2)
                prod_gras += (disconn_gras,)

    return _unique_gras(prod_gras)


def prod_homolytic_scission(gra):
    """ products of homolytic single bond scission
    """

    prod_gras = tuple()

    single_bonds = bonds_of_order(gra, mbond=1)
    for bond in single_bonds:
        gra2 = remove_bonds(gra, [frozenset(bond)])
        disconn_gras = automol.graph.connected_components(gra2)
        prod_gras += (disconn_gras,)

    return _unique_gras(prod_gras)


def _prod_group_loss(gra, grp='H'):
    """ products of hydrogen loss. Need to generalize to group loss
    """

    prod_gras = tuple()

    symb_idx_dct = atom_symbol_idxs(gra)
    h_idxs = symb_idx_dct[grp]

    for idx in h_idxs:
        prod_gras += ((remove_atoms(gra, [idx]),),)

    return _unique_gras(prod_gras)


def _unique_gras(gra_lst):
    """ Determine all of the unique gras deals with gras with multiple components
    """

    uni_gras = tuple()
    if gra_lst:

        # Initialize list with first element
        uni_gras += (gra_lst[0],)

        # Test if the del_gra is isomorphic to any of the uni_del_gras
        for gra in gra_lst[1:]:
            new_uni = True
            for uni_gra in uni_gras:
                if len(gra) == 1:
                    isodct = full_isomorphism(gra[0], uni_gra[0])
                else:
                    isodct = full_isomorphism(union(*gra), union(*uni_gra))
                if isodct:
                    new_uni = False
                    break

            # Add graph and idx to lst if del gra is unique
            if new_uni:
                uni_gras += (gra,)

    return uni_gras
