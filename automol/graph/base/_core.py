""" core graph functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools
import functools
import operator
import numpy
import yaml
import future.moves.itertools as fmit
from phydat import ptab
from phydat import phycon
import automol.formula
from automol import util
from automol.util import dict_
import automol.util.dict_.multi as mdict

ATM_SYM_POS = 0
ATM_IMP_HYD_VLC_POS = 1
ATM_STE_PAR_POS = 2

BND_ORD_POS = 0
BND_STE_PAR_POS = 1

VALENCE_DCT = {
    None: 0,
    1: 1,   # H
    2: 2,   # Be
    13: 3,  # B
    14: 4,  # C
    15: 3,  # N
    16: 2,  # O
    17: 1,  # F
    18: 0,  # He
}

LONE_PAIR_COUNTS_DCT = {
    None: 0,
    1: 0,   # H
    2: 0,   # Be
    13: 0,  # B
    14: 0,  # C
    15: 1,  # N
    16: 2,  # O
    17: 3,  # F
    18: 4,  # He
}

ATM_PROP_NAMES = ('symbol', 'implicit_hydrogen_valence', 'stereo_parity')
BND_PROP_NAMES = ('order', 'stereo_parity')


# # constructors
def from_data(atm_symb_dct, bnd_keys, atm_imp_hyd_vlc_dct=None,
              atm_ste_par_dct=None, bnd_ord_dct=None, bnd_ste_par_dct=None):
    """ Construct a molecular graph from data

    format:
        gra = (atm_dct, bnd_dct)
        atm_dct := {atm_key: (atm_sym, atm_imp_hyd_vlc, atm_ste_par), ...}
        bnd_dct := {bnd_key: (bnd_ord, bnd_ste_par), ...}
        [where bnd_key := frozenset({atm1_key, atm2_key})]

    :param atm_symb_dct: atomic symbols, by atom key
    :type atm_symb_dct: dict
    :param bnd_keys: bond keys
    :type bnd_keys: set
    :param atm_imp_hyd_vlc_dct: the number of implicit hydrogens associated
        with each atom, by atom key
    :type atm_imp_hyd_vlc_dct: dict
    :param atm_ste_par_dct: atom stereo parities, by atom key
    :type atm_ste_par_dct: dict
    :param bnd_ord_dct: bond orders, by bond key
    :type bnd_ord_dct: dict
    :param bnd_ste_par_dct: bond stereo parities, by bond key
    :type bnd_ste_par_dct: dict
    """

    atm_dct = atoms_from_data(
        atm_symb_dct=atm_symb_dct,
        atm_imp_hyd_vlc_dct=atm_imp_hyd_vlc_dct,
        atm_ste_par_dct=atm_ste_par_dct)

    bnd_dct = bonds_from_data(
        bnd_keys=bnd_keys,
        bnd_ord_dct=bnd_ord_dct,
        bnd_ste_par_dct=bnd_ste_par_dct)

    gra = _from_atoms_and_bonds(atm_dct=atm_dct, bnd_dct=bnd_dct)

    return gra


def atoms_from_data(atm_symb_dct, atm_imp_hyd_vlc_dct=None,
                    atm_ste_par_dct=None):
    """ Construct an atom dictionary from constituent data.

        format:
            atm_dct := {atm_key: (atm_sym, atm_imp_hyd_vlc, atm_ste_par), ...}

        :param atm_symb_dct: atomic symbols, by atom key
        :type atm_symb_dct: dict
        :param atm_imp_hyd_vlc_dct: the no. of implicit hydrogens associated
            with each atom, by atom key
        :type atm_imp_hyd_vlc_dct: dict
        :param atm_ste_par_dct: atom stereo parities, by atom key
        :type atm_ste_par_dct: dict
    """

    keys = sorted(atm_symb_dct.keys())
    symbs = dict_.values_by_key(atm_symb_dct, keys)
    vlcs = dict_.values_by_key(
        dict_.empty_if_none(atm_imp_hyd_vlc_dct), keys, fill_val=0)
    pars = dict_.values_by_key(
        dict_.empty_if_none(atm_ste_par_dct), keys, fill_val=None)

    natms = len(symbs)
    vlcs = [0] * natms if vlcs is None else list(vlcs)
    pars = [None] * natms if pars is None else list(pars)

    assert len(vlcs) == natms
    assert len(pars) == natms

    symbs = list(map(ptab.to_symbol, symbs))
    vlcs = list(map(int, vlcs))

    assert all(par in (None, False, True) for par in pars)
    pars = [bool(par) if par is not None else par for par in pars]

    atm_dct = dict(zip(keys, zip(symbs, vlcs, pars)))

    return atm_dct


def bonds_from_data(bnd_keys, bnd_ord_dct=None, bnd_ste_par_dct=None):
    """ construct bond dictionary graph from data

        format:
            bnd_dct := {bnd_key: (bnd_ord, bnd_ste_par), ...}
            [where bnd_key := frozenset({atm1_key, atm2_key})]

        :param bnd_keys: bond keys
        :type bnd_keys: set
        :param bnd_ord_dct: bond orders, by bond key
        :type bnd_ord_dct: dict
        :param bnd_ste_par_dct: bond stereo parities, by bond key
        :type bnd_ste_par_dct: dict
        :rtype: dict[frozenset({int}): tuple(str, str))]
    """

    keys = sorted(bnd_keys)
    assert all(len(key) == 2 for key in keys)
    ords = dict_.values_by_key(
        dict_.empty_if_none(bnd_ord_dct), keys, fill_val=1)
    pars = dict_.values_by_key(
        dict_.empty_if_none(bnd_ste_par_dct), keys, fill_val=None)

    nbnds = len(keys)

    ords = [1] * nbnds if ords is None else list(ords)
    pars = [None] * nbnds if pars is None else list(pars)

    assert len(ords) == nbnds
    assert len(pars) == nbnds

    keys = list(map(frozenset, keys))
    ords = [int(o) if round(o) == o else float(round(o, 1)) for o in ords]

    assert all(par in (None, False, True) for par in pars)
    pars = [bool(par) if par is not None else par for par in pars]

    bnd_dct = dict(zip(keys, zip(ords, pars)))

    return bnd_dct


def from_atoms_and_bonds(atm_dct, bnd_dct):
    """ Construct a molecular graph from atom and bond dictionaries.

        format:
            gra = (atm_dct, bnd_dct)

        :param atm_dct: atom dictionary
        :type atm_dct: dict
        :param bnd_dct: bond dictionary
        :type bnd_dct: dict
        :rtype: (dict, dict)
    """

    atm_dct = dict(atm_dct)
    bnd_dct = dict(bnd_dct)

    atm_symb_dct = dict_.transform_values(atm_dct, lambda x: x[0])
    atm_imp_hyd_vlc_dct = dict_.transform_values(atm_dct, lambda x: x[1])
    atm_ste_par_dct = dict_.transform_values(atm_dct, lambda x: x[2])

    bnd_ord_dct = dict_.transform_values(bnd_dct, lambda x: x[0])
    bnd_ste_par_dct = dict_.transform_values(bnd_dct, lambda x: x[1])

    return from_data(
        atm_symb_dct, bnd_dct.keys(),
        atm_imp_hyd_vlc_dct=atm_imp_hyd_vlc_dct,
        atm_ste_par_dct=atm_ste_par_dct, bnd_ord_dct=bnd_ord_dct,
        bnd_ste_par_dct=bnd_ste_par_dct)


# # getters
def atoms(gra):
    """ atoms, as a dictionary
    """
    atm_dct, _ = gra
    return atm_dct


def bonds(gra):
    """ bonds, as a dictionary
    """
    # print(gra)
    _, bnd_dct = gra
    return bnd_dct


def atom_keys(gra, sym=None, excl_syms=()):
    """ atom keys
    """
    atm_keys = frozenset(atoms(gra).keys())
    if sym is not None:
        atm_sym_dct = atom_symbols(gra)
        atm_keys = frozenset(k for k in atm_keys
                             if atm_sym_dct[k] == sym and
                             atm_sym_dct[k] not in excl_syms)
    return atm_keys


def bond_keys(gra):
    """ bond keys
    """
    return frozenset(bonds(gra).keys())


def atom_symbols(gra):
    """ atom symbols, as a dictionary
    """
    return mdict.by_key_by_position(atoms(gra), atom_keys(gra), ATM_SYM_POS)


def bond_orders(gra):
    """ bond orders, as a dictionary
    """
    return mdict.by_key_by_position(bonds(gra), bond_keys(gra), BND_ORD_POS)


def atom_implicit_hydrogen_valences(gra):
    """ atom implicit hydrogen valences, as a dictionary
    """
    return mdict.by_key_by_position(atoms(gra), atom_keys(gra),
                                    ATM_IMP_HYD_VLC_POS)


def atom_stereo_parities(gra):
    """ atom parities, as a dictionary
    """
    return mdict.by_key_by_position(atoms(gra), atom_keys(gra),
                                    ATM_STE_PAR_POS)


def bond_stereo_parities(gra):
    """ bond parities, as a dictionary
    """
    return mdict.by_key_by_position(bonds(gra), bond_keys(gra),
                                    BND_STE_PAR_POS)


# # setters
def set_atom_symbols(sgr, atm_symb_dct):
    """ set atom parities
    """
    atm_dct = mdict.set_by_key_by_position(atoms(sgr), atm_symb_dct,
                                           ATM_SYM_POS)
    return from_atoms_and_bonds(atm_dct, bonds(sgr))


def set_bond_orders(rgr, bnd_ord_dct):
    """ set bond orders
    """
    bnd_dct = mdict.set_by_key_by_position(bonds(rgr), bnd_ord_dct,
                                           BND_ORD_POS)
    return from_atoms_and_bonds(atoms(rgr), bnd_dct)


def set_atom_implicit_hydrogen_valences(gra, atm_imp_hyd_vlc_dct):
    """ set atom implicit hydrogen valences
    """
    atm_dct = mdict.set_by_key_by_position(atoms(gra), atm_imp_hyd_vlc_dct,
                                           ATM_IMP_HYD_VLC_POS)
    bnd_dct = bonds(gra)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def set_atom_stereo_parities(sgr, atm_par_dct):
    """ set atom parities
    """
    atm_dct = mdict.set_by_key_by_position(atoms(sgr), atm_par_dct,
                                           ATM_STE_PAR_POS)
    return from_atoms_and_bonds(atm_dct, bonds(sgr))


def set_bond_stereo_parities(sgr, bnd_par_dct):
    """ set bond parities
    """
    bnd_dct = mdict.set_by_key_by_position(bonds(sgr), bnd_par_dct,
                                           BND_STE_PAR_POS)
    return from_atoms_and_bonds(atoms(sgr), bnd_dct)


# # I/O
def string(gra, one_indexed=True):
    """ write the graph to a string
    """
    yaml_gra_dct = yaml_dictionary(gra, one_indexed=one_indexed)
    gra_str = yaml.dump(yaml_gra_dct, default_flow_style=None, sort_keys=False)
    return gra_str


def yaml_dictionary(gra, one_indexed=True):
    """ generate a YAML dictionary representing a given graph
    """
    if one_indexed:
        # shift to one-indexing when we print
        atm_key_dct = {atm_key: atm_key+1 for atm_key in atom_keys(gra)}
        gra = relabel(gra, atm_key_dct)

    yaml_atm_dct = atoms(gra)
    yaml_bnd_dct = bonds(gra)

    # prepare the atom dictionary
    yaml_atm_dct = dict(sorted(yaml_atm_dct.items()))
    yaml_atm_dct = dict_.transform_values(
        yaml_atm_dct, lambda x: dict(zip(ATM_PROP_NAMES, x)))

    # perpare the bond dictionary
    yaml_bnd_dct = dict_.transform_keys(
        yaml_bnd_dct, lambda x: tuple(sorted(x)))
    yaml_bnd_dct = dict(sorted(yaml_bnd_dct.items()))
    yaml_bnd_dct = dict_.transform_keys(
        yaml_bnd_dct, lambda x: '-'.join(map(str, x)))
    yaml_bnd_dct = dict_.transform_values(
        yaml_bnd_dct, lambda x: dict(zip(BND_PROP_NAMES, x)))

    yaml_gra_dct = {'atoms': yaml_atm_dct, 'bonds': yaml_bnd_dct}
    return yaml_gra_dct


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

    gra = from_atoms_and_bonds(atm_dct, bnd_dct)

    if one_indexed:
        # revert one-indexing if the input is one-indexed
        atm_key_dct = {atm_key: atm_key-1 for atm_key in atom_keys(gra)}
        gra = relabel(gra, atm_key_dct)

    return gra


# # conversions
def frozen(gra):
    """ hashable, sortable, immutable container of graph data
    """
    atm_keys = sorted(atom_keys(gra))
    bnd_keys = sorted(bond_keys(gra), key=sorted)

    # make it sortable by replacing Nones with -infinity
    atm_vals = numpy.array(dict_.values_by_key(atoms(gra), atm_keys),
                           dtype=numpy.object)
    bnd_vals = numpy.array(dict_.values_by_key(bonds(gra), bnd_keys),
                           dtype=numpy.object)
    atm_vals[numpy.equal(atm_vals, None)] = -numpy.inf
    bnd_vals[numpy.equal(bnd_vals, None)] = -numpy.inf

    frz_atms = tuple(zip(atm_keys, map(tuple, atm_vals)))
    frz_bnds = tuple(zip(bnd_keys, map(tuple, bnd_vals)))
    return (frz_atms, frz_bnds)


def formula(gra):
    """ Generate a stoichiometric formula dictionary from a molecular graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :type: dict[str: int]
    """

    gra = explicit(gra)
    syms = atom_symbols(gra).values()
    fml = util.formula_from_symbols(syms)

    return fml


# # properties
def atom_count(gra, dummy=False, with_implicit=True):
    """ count the number of atoms in this molecule

    by default, this includes implicit hydrogens and excludes dummy atoms
    """
    if not dummy:
        gra = without_dummy_atoms(gra)
    natms = len(atoms(gra))
    if with_implicit:
        atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(gra)
        natms += sum(atm_imp_hyd_vlc_dct.values())
    return natms


def atom_count_by_type(gra, sym, keys=None):
    """ count the number of atoms with a given type (symbol)

    :param gra: the graph
    :param sym: the symbol
    :param keys: optionally, restrict the count to a subset of keys
    """
    keys = atom_keys(gra) if keys is None else keys
    symb_dct = atom_symbols(gra)
    symbs = list(map(symb_dct.__getitem__, keys))
    return symbs.count(sym)


def heavy_atom_count(gra, dummy=False):
    """ the number of heavy atoms
    """
    if not dummy:
        gra = without_dummy_atoms(gra)
    atm_symb_dct = atom_symbols(gra)
    nhvy_atms = sum(ptab.to_number(sym) != 1 for sym in atm_symb_dct.values())
    return nhvy_atms


def electron_count(gra, charge=0):
    """ the number of electrons in the molecule
    """
    atm_symb_dct = atom_symbols(explicit(gra))
    nelec = sum(map(ptab.to_number, atm_symb_dct.values())) - charge
    return nelec


def atom_stereo_keys(sgr):
    """ keys to atom stereo-centers
    """
    atm_ste_keys = dict_.keys_by_value(atom_stereo_parities(sgr),
                                       lambda x: x in [True, False])
    return atm_ste_keys


def bond_stereo_keys(sgr):
    """ keys to bond stereo-centers
    """
    bnd_ste_keys = dict_.keys_by_value(bond_stereo_parities(sgr),
                                       lambda x: x in [True, False])
    return bnd_ste_keys


def has_stereo(gra):
    """ does this graph have stereo of any kind?
    """
    return bool(atom_stereo_keys(gra) or bond_stereo_keys(gra))


def atom_element_valences(gra):
    """ element valences (# possible single bonds), by atom
    """
    atm_symb_dct = atom_symbols(gra)
    atm_group_idx_dct = dict_.transform_values(atm_symb_dct, ptab.to_group)
    atm_elem_vlc_dct = dict_.transform_values(atm_group_idx_dct,
                                              VALENCE_DCT.__getitem__)
    return atm_elem_vlc_dct


def atom_lone_pair_counts(gra):
    """ lone pair counts, by atom
    """
    atm_symb_dct = atom_symbols(gra)
    atm_group_idx_dct = dict_.transform_values(atm_symb_dct, ptab.to_group)
    atm_lpc_dct = dict_.transform_values(atm_group_idx_dct,
                                         LONE_PAIR_COUNTS_DCT.__getitem__)
    atm_lpc_dct = dict_.transform_values(atm_lpc_dct, int)
    return atm_lpc_dct


def atom_van_der_waals_radius(gra, key):
    """ van der waals radius for an atom in the graph (in angstroms)
    """

    symb_dct = atom_symbols(gra)
    symb = symb_dct[key]
    rad = ptab.van_der_waals_radius(symb) * phycon.BOHR2ANG

    return rad


def atom_bond_valences(gra, bond_order=True):
    """ bond count (bond valence), by atom
    """
    atm_keys = list(atom_keys(gra))
    gra = explicit(gra)
    if not bond_order:
        gra = without_bond_orders(gra)

    atm_nbhs = dict_.values_by_key(atom_neighborhoods(gra), atm_keys)
    atm_bnd_vlcs = [sum(bond_orders(nbh).values()) for nbh in atm_nbhs]
    atm_bnd_vlc_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_bnd_vlcs)), int)
    return atm_bnd_vlc_dct


def atom_unsaturated_valences(gra, bond_order=True):
    """ unsaturated valences, by atom

    element valences minus bonding valences = pi sites and radical electrons
    """
    atm_keys = list(atom_keys(gra))
    if not bond_order:
        gra = without_bond_orders(gra)

    atm_bnd_vlcs = dict_.values_by_key(atom_bond_valences(gra), atm_keys)
    atm_tot_vlcs = dict_.values_by_key(atom_element_valences(gra), atm_keys)
    atm_rad_vlcs = numpy.subtract(atm_tot_vlcs, atm_bnd_vlcs)
    atm_unsat_vlc_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_rad_vlcs)), int)
    return atm_unsat_vlc_dct


def atom_explicit_hydrogen_valences(gra):
    """ explicit hydrogen valences, by atom
    """
    return dict_.transform_values(atom_explicit_hydrogen_keys(gra), len)


def atom_hybridizations(rgr):
    """ atom hybridizations, by atom
    """
    rgr = without_fractional_bonds(rgr)
    atm_keys = list(atom_keys(rgr))
    atm_unsat_vlc_dct = atom_unsaturated_valences(rgr, bond_order=True)
    atm_bnd_vlc_dct = atom_bond_valences(rgr, bond_order=False)     # note!!
    atm_unsat_vlcs = numpy.array(
        dict_.values_by_key(atm_unsat_vlc_dct, atm_keys))
    atm_bnd_vlcs = numpy.array(dict_.values_by_key(atm_bnd_vlc_dct, atm_keys))
    atm_lpcs = numpy.array(
        dict_.values_by_key(atom_lone_pair_counts(rgr), atm_keys))
    atm_hybs = atm_unsat_vlcs + atm_bnd_vlcs + atm_lpcs - 1
    atm_hyb_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_hybs)), int)
    return atm_hyb_dct


def maximum_spin_multiplicity(gra, bond_order=True):
    """ the highest possible spin multiplicity for this molecular graph
    """
    atm_rad_vlc_dct = atom_unsaturated_valences(gra, bond_order=bond_order)
    return sum(atm_rad_vlc_dct.values()) + 1


def possible_spin_multiplicities(gra, bond_order=True):
    """ possible spin multiplicities for this molecular graph
    """
    mult_max = maximum_spin_multiplicity(gra, bond_order=bond_order)
    mult_min = 2 if mult_max % 2 == 0 else 1
    mults = tuple(range(mult_min, mult_max+1, 2))
    return mults


def atom_symbol_keys(gra):
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


def backbone_keys(gra):
    """ backbone atom keys
    """
    bbn_keys = atom_keys(gra) - explicit_hydrogen_keys(gra)
    return bbn_keys


def atom_explicit_hydrogen_keys(gra):
    """ explicit hydrogen valences, by atom
    """
    exp_hyd_keys = explicit_hydrogen_keys(gra)
    atm_exp_hyd_keys_dct = dict_.transform_values(
        atoms_neighbor_atom_keys(gra), lambda x: x & exp_hyd_keys)
    return atm_exp_hyd_keys_dct


def explicit_hydrogen_keys(gra):
    """ explicit hydrogen keys (H types: explicit, implicit, backbone)
    """
    hyd_keys = dict_.keys_by_value(atom_symbols(gra), lambda x: x == 'H')
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    def _is_backbone(hyd_key):
        is_h2 = all(ngb_key in hyd_keys and hyd_key < ngb_key
                    for ngb_key in atm_ngb_keys_dct[hyd_key])
        is_multivalent = len(atm_ngb_keys_dct[hyd_key]) > 1
        return is_h2 or is_multivalent

    exp_hyd_keys = frozenset(fmit.filterfalse(_is_backbone, hyd_keys))
    return exp_hyd_keys


def terminal_heavy_atom_keys(gra):
    """ terminal heavy atoms, sorted by atom type and hydrogen count
    """
    gra = implicit(gra)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(gra)
    atm_keys = [key for key, ngb_keys in atoms_neighbor_atom_keys(gra).items()
                if len(ngb_keys) == 1]
    atm_keys = sorted(atm_keys, key=atm_imp_hyd_vlc_dct.__getitem__,
                      reverse=True)
    atm_symbs = dict_.values_by_key(atom_symbols(gra), atm_keys)
    srt = automol.formula.argsort_symbols(atm_symbs, symbs_first=('C',))
    atm_keys = tuple(map(atm_keys.__getitem__, srt))
    return atm_keys


def unsaturated_atom_keys(gra):
    """ keys of unsaturated (radical or pi-bonded) atoms
    """
    atm_unsat_vlc_dct = atom_unsaturated_valences(gra, bond_order=False)
    unsat_atm_keys = frozenset(dict_.keys_by_value(atm_unsat_vlc_dct, bool))
    return unsat_atm_keys


def angle_keys(gra):
    """ triples of keys for pairs of adjacent bonds, with the central atom in
    the middle
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

    return tuple(ang_keys)


# # relabeling and changing keys
def relabel(gra, atm_key_dct):
    """ relabel the graph with new atom keys
    """
    orig_atm_keys = atom_keys(gra)
    assert set(atm_key_dct.keys()) <= orig_atm_keys, (
        '{}\n{}'.format(set(atm_key_dct.keys()), orig_atm_keys)
    )

    new_atm_key_dct = dict(zip(orig_atm_keys, orig_atm_keys))
    new_atm_key_dct.update(atm_key_dct)

    _relabel_atom_key = new_atm_key_dct.__getitem__

    def _relabel_bond_key(bnd_key):
        return frozenset(map(_relabel_atom_key, bnd_key))

    atm_dct = dict_.transform_keys(atoms(gra), _relabel_atom_key)
    bnd_dct = dict_.transform_keys(bonds(gra), _relabel_bond_key)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def transform_keys(gra, atm_key_func):
    """ transform atom keys with a function
    """
    atm_keys = atom_keys(gra)
    atm_key_dct = dict(zip(atm_keys, map(atm_key_func, atm_keys)))
    return relabel(gra, atm_key_dct)


def standard_keys(gra):
    """ replace the current atom keys with standard indices, counting from zero
    """
    atm_key_dct = dict(map(reversed, enumerate(sorted(atom_keys(gra)))))
    return relabel(gra, atm_key_dct)


def standard_keys_for_sequence(gras):
    """ assigns non-overlapping keys to a sequence of graphs

    (returns a series of key maps for each)
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
    """ relabel a geometry graph to line up with a z-matrix

    Inserts dummy atoms and sorts/relabels in z-matrix order.

    Graph keys should correspond to the geometry used for conversion.

    :param gra: the graph
    :param zma_keys: graph keys in the order they appear in the z-matrix
    :param dummy_key_dct: dummy keys introduced on z-matrix conversion, by atom
        they are attached to
    """
    gra = add_dummy_atoms(gra, dummy_key_dct)
    key_dct = dict(map(reversed, enumerate(zma_keys)))
    gra = relabel(gra, key_dct)
    return gra


def relabel_for_geometry(gra):
    """ relabel a z-matrix graph to line up with a geometry

    The result will line up with a geometry converted from the z-matrix, with
    dummy atoms removed.

    Removes dummy atoms and relabels in geometry order.

    Graph keys should correspond to the z-matrix used for conversion.

    :param gra: the graph
    """
    dummy_keys = sorted(atom_keys(gra, sym='X'))

    for dummy_key in reversed(dummy_keys):
        gra = _shift_remove_dummy_atom(gra, dummy_key)
    return gra


def _shift_remove_dummy_atom(gra, dummy_key):
    keys = sorted(atom_keys(gra))
    idx = keys.index(dummy_key)
    key_dct = {}
    key_dct.update({k: k for k in keys[:idx]})
    key_dct.update({k: k-1 for k in keys[(idx+1):]})
    gra = remove_atoms(gra, [dummy_key])
    gra = relabel(gra, key_dct)
    return gra


# # add/remove/insert/without
def add_atoms(gra, symb_dct, imp_hyd_vlc_dct=None, ste_par_dct=None):
    """ add atoms to this molecular graph, setting their keys
    """
    atm_keys = atom_keys(gra)
    atm_symb_dct = atom_symbols(gra)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(gra)
    atm_ste_par_dct = atom_stereo_parities(gra)

    keys = set(symb_dct.keys())
    imp_hyd_vlc_dct = {} if imp_hyd_vlc_dct is None else imp_hyd_vlc_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct

    assert not keys & atm_keys
    assert set(imp_hyd_vlc_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys

    atm_symb_dct.update(symb_dct)
    atm_imp_hyd_vlc_dct.update(imp_hyd_vlc_dct)
    atm_ste_par_dct.update(ste_par_dct)

    atm_dct = atoms_from_data(
        atm_symb_dct=atm_symb_dct, atm_imp_hyd_vlc_dct=atm_imp_hyd_vlc_dct,
        atm_ste_par_dct=atm_ste_par_dct)
    bnd_dct = bonds(gra)
    gra = from_atoms_and_bonds(atm_dct=atm_dct, bnd_dct=bnd_dct)
    return gra


def add_bonds(gra, keys, ord_dct=None, ste_par_dct=None, check=True):
    """ add bonds to this molecular graph
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
        assert not keys & bnd_keys, (
            '{} and {} have a non-empty intersection'.format(keys, bnd_keys)
        )

    assert set(ord_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys

    bnd_keys.update(keys)
    bnd_ord_dct.update(ord_dct)
    bnd_ste_par_dct.update(ste_par_dct)

    atm_dct = atoms(gra)
    bnd_dct = bonds_from_data(
        bnd_keys=bnd_keys, bnd_ord_dct=bnd_ord_dct,
        bnd_ste_par_dct=bnd_ste_par_dct)

    gra = from_atoms_and_bonds(atm_dct=atm_dct, bnd_dct=bnd_dct)
    return gra


def remove_atoms(gra, atm_keys, check=True, stereo=False):
    """ remove atoms from the molecular graph
    """
    all_atm_keys = atom_keys(gra)
    atm_keys = set(atm_keys)

    if check:
        assert atm_keys <= all_atm_keys

    atm_keys_left = all_atm_keys - atm_keys
    return subgraph(gra, atm_keys_left, stereo=stereo)


def remove_bonds(gra, bnd_keys, check=True):
    """ remove bonds from the molecular graph
    """
    all_bnd_keys = bond_keys(gra)
    bnd_keys = set(map(frozenset, bnd_keys))

    if check:
        assert bnd_keys <= all_bnd_keys

    bnd_keys = all_bnd_keys - bnd_keys
    atm_dct = atoms(gra)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


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


def add_atom_implicit_hydrogen_valences(gra, inc_atm_imp_hyd_vlc_dct):
    """ add atom imlicit hydrogen valences

    (increments can be positive or negative)
    """
    atm_keys = list(inc_atm_imp_hyd_vlc_dct.keys())
    atm_imp_hyd_vlcs = numpy.add(
        dict_.values_by_key(atom_implicit_hydrogen_valences(gra), atm_keys),
        dict_.values_by_key(inc_atm_imp_hyd_vlc_dct, atm_keys))
    assert all(atm_imp_hyd_vlc >= 0 for atm_imp_hyd_vlc in atm_imp_hyd_vlcs)
    atm_imp_hyd_vlc_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_imp_hyd_vlcs)), int)
    return set_atom_implicit_hydrogen_valences(gra, atm_imp_hyd_vlc_dct)


def add_atom_explicit_hydrogen_keys(gra, atm_exp_hyd_keys_dct):
    """ add explicit hydrogens by atom
    """
    assert set(atm_exp_hyd_keys_dct.keys()) <= atom_keys(gra), (
        '{} !<= {}'.format(
            set(atm_exp_hyd_keys_dct.keys()), atom_keys(gra))
    )
    for atm_key, atm_exp_hyd_keys in atm_exp_hyd_keys_dct.items():
        assert not set(atm_exp_hyd_keys) & atom_keys(gra)
        atm_exp_hyd_bnd_keys = {frozenset({atm_key, atm_exp_hyd_key})
                                for atm_exp_hyd_key in atm_exp_hyd_keys}
        atm_exp_hyd_symb_dct = dict_.by_key({}, atm_exp_hyd_keys, fill_val='H')
        gra = add_atoms(gra, atm_exp_hyd_symb_dct)
        gra = add_bonds(gra, atm_exp_hyd_bnd_keys)
    return gra


def add_bonded_atom(gra, sym, atm_key, bnd_atm_key=None, imp_hyd_vlc=None,
                    atm_ste_par=None, bnd_ord=None, bnd_ste_par=None):
    """ add a single atom with a bond to an atom already in the graph
    """
    atm_keys = atom_keys(gra)

    bnd_atm_key = max(atm_keys) + 1 if bnd_atm_key is None else bnd_atm_key

    symb_dct = {bnd_atm_key: sym}
    imp_hyd_vlc_dct = ({bnd_atm_key: imp_hyd_vlc}
                       if imp_hyd_vlc is not None else None)
    atm_ste_par_dct = ({bnd_atm_key: atm_ste_par}
                       if atm_ste_par is not None else None)

    gra = add_atoms(gra, symb_dct, imp_hyd_vlc_dct=imp_hyd_vlc_dct,
                    ste_par_dct=atm_ste_par_dct)

    bnd_key = frozenset({bnd_atm_key, atm_key})
    bnd_ord_dct = {bnd_key: bnd_ord} if bnd_ord is not None else None
    bnd_ste_par_dct = ({bnd_key: bnd_ste_par}
                       if bnd_ste_par is not None else None)

    gra = add_bonds(gra, [bnd_key], ord_dct=bnd_ord_dct,
                    ste_par_dct=bnd_ste_par_dct)

    return gra


def insert_bonded_atom(gra, sym, atm_key, bnd_atm_key=None, imp_hyd_vlc=None,
                       atm_ste_par=None, bnd_ord=None, bnd_ste_par=None):
    """ insert a single atom with a bond to an atom already in the graph

    Keys will be standardized upon insertion
    """
    keys = sorted(atom_keys(gra))
    bnd_atm_key_ = max(keys) + 1

    gra = add_bonded_atom(gra, sym, atm_key, bnd_atm_key=bnd_atm_key_,
                          imp_hyd_vlc=imp_hyd_vlc, atm_ste_par=atm_ste_par,
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
    """ add dummy atoms to the graph, with dummy bonds to particular atoms

    :param dummy_key_dct: keys are atoms in the graph on which to place a dummy
        atom; values are the desired keys of the dummy atoms themselves, which
        must not overlap with already existing atoms
    """
    atm_keys = atom_keys(gra)
    assert set(dummy_key_dct.keys()) <= atm_keys, (
        "Keys must be existing atoms in the graph.")
    assert not set(dummy_key_dct.values()) & atm_keys, (
        "Dummy atom keys cannot overlap with existing atoms.")

    for key, dummy_key in sorted(dummy_key_dct.items()):
        gra = add_bonded_atom(gra, 'X', key, bnd_atm_key=dummy_key, bnd_ord=0)

    return gra


def insert_dummy_atoms(gra, dummy_key_dct):
    """ add dummy atoms to the graph, with dummy bonds to particular atoms

    :param dummy_key_dct: keys are atoms in the graph on which to place a dummy
        atom; values are the desired keys of the dummy atoms themselves, which
        must not overlap with already existing atoms
    """
    for key, dummy_key in reversed(sorted(dummy_key_dct.items())):
        # If the dummy key comes first, the anchor key will be offset by 1
        if dummy_key < key:
            key -= 1

        gra = insert_bonded_atom(
            gra, 'X', key, bnd_atm_key=dummy_key, bnd_ord=0)

    return gra


def standard_keys_without_dummy_atoms(gra):
    """ remove dummy atoms and standardize keys, returning the dummy key
    dictionary for converting back

    Requires that graph follows z-matrix ordering (this is checked)
    """
    dummy_ngb_key_dct = dummy_atoms_neighbor_atom_key(gra)

    dummy_keys_dct = {}
    last_dummy_key = -1
    decr = 0
    for dummy_key, key in sorted(dummy_ngb_key_dct.items()):
        assert last_dummy_key <= key, (
            "{:d} must follow previous dummy {:d}"
            .format(key, last_dummy_key))

        dummy_keys_dct[key-decr] = dummy_key-decr
        gra = remove_atoms(gra, [dummy_key])

        decr += 1
        last_dummy_key = dummy_key-decr

    gra = standard_keys(gra)
    return gra, dummy_keys_dct


def without_bond_orders(gra):
    """ resonance graph with maximum spin (i.e. no pi bonds)
    """
    bnd_keys = list(bond_keys(gra))
    # don't set dummy bonds to one!
    bnd_ord_dct = bond_orders(gra)
    bnd_vals = [1 if v != 0 else 0
                for v in map(bnd_ord_dct.__getitem__, bnd_keys)]
    bnd_ord_dct = dict(zip(bnd_keys, bnd_vals))
    return set_bond_orders(gra, bnd_ord_dct)


def without_dummy_atoms(gra):
    """ remove dummy atoms from the graph
    """
    atm_symb_dct = atom_symbols(gra)
    atm_keys = [key for key, sym in atm_symb_dct.items()
                if ptab.to_number(sym)]
    return subgraph(gra, atm_keys, stereo=True)


def without_fractional_bonds(gra):
    """ rounds fractional bonds in the graph
    """
    ord_dct = dict_.transform_values(bond_orders(gra), func=round)
    gra = set_bond_orders(gra, ord_dct)
    return gra


def without_dummy_bonds(gra):
    """ remove 0-order bonds from the graph
    """
    ord_dct = dict_.filter_by_value(bond_orders(gra), func=lambda x: x == 0)
    gra = remove_bonds(gra, ord_dct.keys())
    return gra


def without_stereo_parities(gra):
    """ graph with stereo assignments wiped out
    """
    atm_ste_par_dct = dict_.by_key({}, atom_keys(gra), fill_val=None)
    bnd_ste_par_dct = dict_.by_key({}, bond_keys(gra), fill_val=None)
    gra = set_atom_stereo_parities(gra, atm_ste_par_dct)
    gra = set_bond_stereo_parities(gra, bnd_ste_par_dct)
    return gra


def explicit(gra, atm_keys=None):
    """ make the hydrogens at these atoms explicit
    """
    atm_keys = backbone_keys(gra) if atm_keys is None else atm_keys
    atm_keys = sorted(atm_keys)
    atm_imp_hyd_vlc_dct = dict_.by_key(
        atom_implicit_hydrogen_valences(gra), atm_keys)

    atm_exp_hyd_keys_dct = {}
    next_atm_key = max(atom_keys(gra)) + 1
    for atm_key in atm_keys:
        imp_hyd_vlc = atm_imp_hyd_vlc_dct[atm_key]
        atm_exp_hyd_keys_dct[atm_key] = set(
            range(next_atm_key, next_atm_key+imp_hyd_vlc))
        next_atm_key += imp_hyd_vlc

    gra = set_atom_implicit_hydrogen_valences(
        gra, dict_.by_key({}, atm_keys, fill_val=0))
    gra = add_atom_explicit_hydrogen_keys(gra, atm_exp_hyd_keys_dct)
    return gra


def implicit(gra, atm_keys=None):
    """ make the hydrogens at these atoms implicit
    """
    atm_keys = backbone_keys(gra) if atm_keys is None else atm_keys

    atm_exp_hyd_keys_dct = dict_.by_key(
        atom_explicit_hydrogen_keys(gra), atm_keys)

    inc_imp_hyd_keys_dct = dict_.transform_values(atm_exp_hyd_keys_dct, len)
    gra = add_atom_implicit_hydrogen_valences(gra, inc_imp_hyd_keys_dct)

    exp_hyd_keys = set(itertools.chain(*atm_exp_hyd_keys_dct.values()))
    gra = remove_atoms(gra, exp_hyd_keys, stereo=True)
    return gra


# # unions
def union(gra1, gra2, check=True):
    """ a union of two graphs
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


def union_from_sequence(gras, check=True):
    """ a union of all parts of a sequence of graphs
    """
    def _union(gra1, gra2):
        return union(gra1, gra2, check=check)

    return tuple(functools.reduce(_union, gras))


# # subgraphs and neighborhoods
def subgraph(gra, atm_keys, stereo=False):
    """ the subgraph induced by a subset of the atoms

    :param gra: the graph
    :param atm_keys: the atom keys to be included in the subgraph
    :param stereo: whether or not to include stereo in the subgraph
    :returns: the subgraph
    """
    atm_keys = set(atm_keys)
    assert atm_keys <= atom_keys(gra)
    bnd_keys = set(filter(lambda x: x <= atm_keys, bond_keys(gra)))
    atm_dct = dict_.by_key(atoms(gra), atm_keys)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    sub = from_atoms_and_bonds(atm_dct, bnd_dct)
    if not stereo:
        sub = without_stereo_parities(sub)
    return sub


def bond_induced_subgraph(gra, bnd_keys, stereo=False):
    """ the subgraph induced by a subset of the bonds
    """
    atm_keys = set(itertools.chain(*bnd_keys))
    bnd_keys = set(bnd_keys)
    assert atm_keys <= atom_keys(gra)
    atm_dct = dict_.by_key(atoms(gra), atm_keys)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    sub = from_atoms_and_bonds(atm_dct, bnd_dct)
    if not stereo:
        sub = without_stereo_parities(sub)
    return sub


def atom_neighborhood(gra, atm_key, bnd_keys=None, stereo=False):
    """ neighborhood subgraph for a specific atom
    """
    bnd_keys = bond_keys(gra) if bnd_keys is None else bnd_keys
    nbh_bnd_keys = set(k for k in bnd_keys if atm_key in k)
    nbh = bond_induced_subgraph(gra, nbh_bnd_keys, stereo=stereo)
    return nbh


def atom_neighborhoods(gra, stereo=False):
    """ neighborhood subgraphs, by atom
    """
    bnd_keys = bond_keys(gra)

    def _neighborhood(atm_key):
        return atom_neighborhood(gra, atm_key, bnd_keys=bnd_keys,
                                 stereo=stereo)

    atm_keys = list(atom_keys(gra))
    atm_nbh_dct = dict(zip(atm_keys, map(_neighborhood, atm_keys)))
    return atm_nbh_dct


def bond_neighborhood(gra, bnd_key, bnd_keys=None, stereo=False):
    """ neighborhood subgraph for a specific bond

    :param gra: the graph
    :param bnd_key: the bond key
    :type bnd_key: frozenset[int]
    :param bnd_keys: optionally, restrict this to a subset of the bond keys
    :tupe bnd_keys: tuple[frozenset[int]]
    :returns: the neighborhood subgraph
    """
    bnd_keys = bond_keys(gra) if bnd_keys is None else bnd_keys
    nbh_bnd_keys = set(filter(lambda x: bnd_key & x, bnd_keys))
    nbh = bond_induced_subgraph(gra, nbh_bnd_keys, stereo=stereo)
    return nbh


def bond_neighborhoods(gra, stereo=False):
    """ neighborhood subgraphs, by bond
    """
    bnd_keys = list(bond_keys(gra))

    def _neighborhood(bnd_key):
        return bond_neighborhood(gra, bnd_key, bnd_keys=bnd_keys,
                                 stereo=stereo)

    bnd_nbh_dct = dict(zip(bnd_keys, map(_neighborhood, bnd_keys)))
    return bnd_nbh_dct


def atom_neighbor_atom_key(gra, atm_key, excl_atm_keys=(), incl_atm_keys=None,
                           symbs_first=('C',), symbs_last=('H',)):
    """ get the next in a sorted list of neighbor keys, excluding some
    """
    atm_keys = atom_sorted_neighbor_atom_keys(
        gra, atm_key, excl_atm_keys=excl_atm_keys, incl_atm_keys=incl_atm_keys,
        symbs_first=symbs_first, symbs_last=symbs_last)
    return atm_keys[0] if atm_keys else None


def atoms_neighbor_atom_keys(gra):
    """ keys of neighboring atoms, by atom
    """
    def _neighbor_keys(atm_key, atm_nbh):
        return frozenset(atom_keys(atm_nbh) - {atm_key})

    atm_ngb_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(gra), _neighbor_keys)
    return atm_ngb_keys_dct


def atom_sorted_neighbor_atom_keys(gra, atm_key, excl_atm_keys=(),
                                   incl_atm_keys=None, symbs_first=('C',),
                                   symbs_last=('H',)):
    """ get the next in a sorted list of neighbor keys, excluding some
    """
    atm_symb_dct = atom_symbols(gra)
    incl_atm_keys = atom_keys(gra) if incl_atm_keys is None else incl_atm_keys

    atm_nbh = atom_neighborhood(gra, atm_key)
    atm_keys = sorted(atom_keys(atm_nbh) - {atm_key} - set(excl_atm_keys))
    atm_keys = [k for k in atm_keys if k in incl_atm_keys]

    symbs = list(map(atm_symb_dct.__getitem__, atm_keys))
    srt = automol.formula.argsort_symbols(symbs, symbs_first, symbs_last)
    atm_keys = tuple(map(atm_keys.__getitem__, srt))
    return atm_keys


def atoms_sorted_neighbor_atom_keys(gra, symbs_first=('C',), symbs_last=('H',),
                                    ords_last=(0.1,), prioritize_keys=()):
    """ keys of neighboring atoms, by atom

    :param gra: the graph
    :param symbs_first: atomic symbols to put put first in the sort order
    :param symbs_last: atomic symbols to put last in the sort order
    :param ords_last: neighors connected with a bond of this order will be put
        last in the sort order
    :param prioritize_keys: keys to put first no matter what
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
        atom_neighborhoods(gra), _neighbor_keys)
    return atm_ngb_keys_dct


def atoms_second_degree_neighbor_atom_keys(gra):
    """ keys of second-degree neighboring atoms, by atom

    That is, atoms that are connected through an intermediate atom
    """
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    atm_ngb2_keys_dct = {}
    for atm_key, atm_ngb_keys in atm_ngb_keys_dct.items():
        # Union of neighbors of neighbors
        atm_ngb2_keys = functools.reduce(
            operator.or_, map(atm_ngb_keys_dct.__getitem__, atm_ngb_keys))
        # Subtract of the atom itself and its neighbors (in case of 3-rings)
        atm_ngb2_keys -= {atm_key} | atm_ngb_keys

        atm_ngb2_keys_dct[atm_key] = frozenset(atm_ngb2_keys)
    return atm_ngb2_keys_dct


def atoms_bond_keys(gra):
    """ bond keys, by atom
    """
    return dict_.transform_values(atom_neighborhoods(gra), bond_keys)


def dummy_atoms_neighbor_atom_key(gra):
    """ Atoms that are connected to dummy atoms, by dummy atom key

    (Requires that each dummy atom only be connected to one neighbor)
    """
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    dummy_atm_keys = atom_keys(gra, sym='X')

    dummy_ngb_key_dct = {}
    for key in dummy_atm_keys:
        ngb_keys = atm_ngb_keys_dct[key]
        assert len(ngb_keys) == 1, (
            "Dummy atoms should only be connected to one atom!")
        ngb_key, = ngb_keys
        dummy_ngb_key_dct[key] = ngb_key

    return dummy_ngb_key_dct


def bonds_neighbor_atom_keys(gra):
    """ keys of neighboring atoms, by bond
    """
    def _neighbor_keys(bnd_key, bnd_nbh):
        return frozenset(atom_keys(bnd_nbh) - bnd_key)

    bnd_ngb_keys_dct = dict_.transform_items_to_values(
        bond_neighborhoods(gra), _neighbor_keys)
    return bnd_ngb_keys_dct


def bonds_neighbor_bond_keys(gra):
    """ keys of neighboring bonds, by bond
    """
    def _neighbor_keys(bnd_key, bnd_nbh):
        bnd_keys = bond_keys(bnd_nbh)
        bnd_keys -= {bnd_key}
        bnd_keys = frozenset(key for key in bnd_keys if key & bnd_key)
        return bnd_keys

    bnd_ngb_keys_dct = dict_.transform_items_to_values(
        bond_neighborhoods(gra), _neighbor_keys)
    return bnd_ngb_keys_dct


# # helpers
def _from_atoms_and_bonds(atm_dct, bnd_dct):
    """ Construct a molecular graph from atom and bond dictionaries.

        format:
            gra = (atm_dct, bnd_dct)

        :param atm_dct: atom dictionary
        :type atm_dct: dict
        :param bnd_dct: bond dictionary
        :type bnd_dct: dict
        :rtype: (dict, dict)
    """

    atm_dct = dict(atm_dct)
    bnd_dct = dict(bnd_dct)
    atm_keys = set(atm_dct.keys())
    bnd_keys = set(bnd_dct.keys())
    assert all(bnd_key <= atm_keys for bnd_key in bnd_keys)

    return (atm_dct, bnd_dct)
