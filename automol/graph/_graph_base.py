""" base molecular graph library
"""
import yaml
from automol.util import dict_
import automol.create.graph as _create
from automol.graph._graph_dep import atom_keys
from automol.graph._graph_dep import atom_symbols
from automol.graph._graph_dep import set_atom_stereo_parities
from automol.graph._graph_dep import set_bond_stereo_parities
from automol.graph._graph_dep import relabel
# from automol.graph._graph_dep import ATM_SYM_POS
# from automol.graph._graph_dep import ATM_IMP_HYD_VLC_POS
# from automol.graph._graph_dep import ATM_STE_PAR_POS
# from automol.graph._graph_dep import BND_SYM_POS
# from automol.graph._graph_dep import BND_STE_PAR_POS

ATM_PROP_NAMES = ('symbol', 'implicit_hydrogen_valence', 'stereo_parity')
BND_PROP_NAMES = ('order', 'stereo_parity')


# setters
def atom_symbol_idxs(gra):
    """ determine the indices for each atom symbol
        :return: dict[symb] = [idxs]
    """

    idx_symb_dct = atom_symbols(gra)

    symb_idx_dct = {}
    for idx, symb in idx_symb_dct.items():
        if symb not in symb_idx_dct:
            symb_idx_dct[symb] = [idx]
        else:
            symb_idx_dct[symb].append(idx)

    return symb_idx_dct


def remove_atom_stereo_parities(gra, atm_keys):
    """ Remove stereo parities for certain atoms

    :param gra: the graph
    :param atm_keys: the keys for which to remove stereo parities
    :param gra: the resulting graph
    """
    return set_atom_stereo_parities(gra, {k: None for k in atm_keys})


def remove_bond_stereo_parities(gra, bnd_keys):
    """ Remove stereo parities for certain bonds

    :param gra: the graph
    :param bnd_keys: the keys for which to remove stereo parities
    :param gra: the resulting graph
    """
    return set_bond_stereo_parities(gra, {k: None for k in bnd_keys})


# I/O
def from_string(gra_str, one_indexed=True):
    """ read the graph from a string
    """
    yaml_gra_dct = yaml.load(gra_str, Loader=yaml.FullLoader)
    gra = from_yaml_dictionary(yaml_gra_dct, one_indexed=one_indexed)
    return gra


def from_yaml_dictionary(yaml_gra_dct, one_indexed=True):
    """ read the graph from a yaml dictionary
    """
    atm_dct = yaml_gra_dct['atoms']
    bnd_dct = yaml_gra_dct['bonds']

    atm_dct = dict_.transform_values(
        atm_dct, lambda x: tuple(map(x.__getitem__, ATM_PROP_NAMES)))

    bnd_dct = dict_.transform_keys(
        bnd_dct, lambda x: frozenset(map(int, x.split('-'))))

    bnd_dct = dict_.transform_values(
        bnd_dct, lambda x: tuple(map(x.__getitem__, BND_PROP_NAMES)))

    gra = _create.from_atoms_and_bonds(atm_dct, bnd_dct)

    if one_indexed:
        # revert one-indexing if the input is one-indexed
        atm_key_dct = {atm_key: atm_key-1 for atm_key in atom_keys(gra)}
        gra = relabel(gra, atm_key_dct)

    return gra
