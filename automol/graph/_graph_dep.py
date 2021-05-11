""" objects that everything else depends on
"""
import itertools
import functools
import numpy
import yaml
import future.moves.itertools as fmit
from phydat import ptab
import automol.util.dict_.multi as mdict
from automol.util import dict_
import automol.create.graph as _create
from automol.formula import argsort_symbols

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


# Basic graph objects
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


# Transformations
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


def remove_atoms(gra, atm_keys, check=True, stereo=False):
    """ remove atoms from the molecular graph
    """
    all_atm_keys = atom_keys(gra)
    atm_keys = set(atm_keys)

    if check:
        assert atm_keys <= all_atm_keys

    atm_keys_left = all_atm_keys - atm_keys
    return subgraph(gra, atm_keys_left, stereo=stereo)


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

    atm_dct = _create.atoms_from_data(
        atom_symbols=atm_symb_dct,
        atom_implicit_hydrogen_valences=atm_imp_hyd_vlc_dct,
        atom_stereo_parities=atm_ste_par_dct)
    bnd_dct = bonds(gra)
    gra = _create.from_atoms_and_bonds(atoms=atm_dct, bonds=bnd_dct)
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
    bnd_dct = _create.bonds_from_data(
        bond_keys=bnd_keys, bond_orders=bnd_ord_dct,
        bond_stereo_parities=bnd_ste_par_dct)

    gra = _create.from_atoms_and_bonds(atoms=atm_dct, bonds=bnd_dct)
    return gra


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
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


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


def set_bond_orders(rgr, bnd_ord_dct):
    """ set bond orders
    """
    bnd_dct = mdict.set_by_key_by_position(bonds(rgr), bnd_ord_dct,
                                           BND_ORD_POS)
    return _create.from_atoms_and_bonds(atoms(rgr), bnd_dct)


def atom_explicit_hydrogen_keys(gra):
    """ explicit hydrogen valences, by atom
    """
    exp_hyd_keys = explicit_hydrogen_keys(gra)
    atm_exp_hyd_keys_dct = dict_.transform_values(
        atoms_neighbor_atom_keys(gra), lambda x: x & exp_hyd_keys)
    return atm_exp_hyd_keys_dct


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


# # other properties
def backbone_keys(gra):
    """ backbone atom keys
    """
    bbn_keys = atom_keys(gra) - explicit_hydrogen_keys(gra)
    return bbn_keys


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


def maximum_spin_multiplicity(gra, bond_order=True):
    """ the highest possible spin multiplicity for this molecular graph
    """
    atm_rad_vlc_dct = atom_unsaturated_valences(gra, bond_order=bond_order)
    return sum(atm_rad_vlc_dct.values()) + 1


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


def without_stereo_parities(gra):
    """ graph with stereo assignments wiped out
    """
    atm_ste_par_dct = dict_.by_key({}, atom_keys(gra), fill_val=None)
    bnd_ste_par_dct = dict_.by_key({}, bond_keys(gra), fill_val=None)
    gra = set_atom_stereo_parities(gra, atm_ste_par_dct)
    gra = set_bond_stereo_parities(gra, bnd_ste_par_dct)
    return gra


def set_atom_stereo_parities(sgr, atm_par_dct):
    """ set atom parities
    """
    atm_dct = mdict.set_by_key_by_position(atoms(sgr), atm_par_dct,
                                           ATM_STE_PAR_POS)
    return _create.from_atoms_and_bonds(atm_dct, bonds(sgr))


def set_bond_stereo_parities(sgr, bnd_par_dct):
    """ set bond parities
    """
    bnd_dct = mdict.set_by_key_by_position(bonds(sgr), bnd_par_dct,
                                           BND_STE_PAR_POS)
    return _create.from_atoms_and_bonds(atoms(sgr), bnd_dct)


# neighboorhoods and subcomponents
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
    srt = argsort_symbols(symbs, symbs_first, symbs_last)
    atm_keys = tuple(map(atm_keys.__getitem__, srt))
    return atm_keys


def atom_element_valences(gra):
    """ element valences (# possible single bonds), by atom
    """
    atm_symb_dct = atom_symbols(gra)
    atm_group_idx_dct = dict_.transform_values(atm_symb_dct, ptab.to_group)
    atm_elem_vlc_dct = dict_.transform_values(atm_group_idx_dct,
                                              VALENCE_DCT.__getitem__)
    return atm_elem_vlc_dct


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


def bond_induced_subgraph(gra, bnd_keys, stereo=False):
    """ the subgraph induced by a subset of the bonds
    """
    atm_keys = set(itertools.chain(*bnd_keys))
    bnd_keys = set(bnd_keys)
    assert atm_keys <= atom_keys(gra)
    atm_dct = dict_.by_key(atoms(gra), atm_keys)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    sub = _create.from_atoms_and_bonds(atm_dct, bnd_dct)
    if not stereo:
        sub = without_stereo_parities(sub)
    return sub


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
    sub = _create.from_atoms_and_bonds(atm_dct, bnd_dct)
    if not stereo:
        sub = without_stereo_parities(sub)
    return sub


def set_atom_implicit_hydrogen_valences(gra, atm_imp_hyd_vlc_dct):
    """ set atom implicit hydrogen valences
    """
    atm_dct = mdict.set_by_key_by_position(atoms(gra), atm_imp_hyd_vlc_dct,
                                           ATM_IMP_HYD_VLC_POS)
    bnd_dct = bonds(gra)
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


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


# stereo
def has_stereo(gra):
    """ does this graph have stereo of any kind?
    """
    return bool(atom_stereo_keys(gra) or bond_stereo_keys(gra))


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


# resonances
def dominant_resonance(rgr):
    """ *a* dominant (minimum spin/maximum pi) resonance graph
    """
    return next(iter(dominant_resonances(rgr)))


def dominant_resonances(rgr):
    """ all dominant (minimum spin/maximum pi) resonance graphs
    """
    rgr = without_fractional_bonds(rgr)
    rgrs = resonances(rgr)
    mult_min = min(map(maximum_spin_multiplicity, rgrs))
    dom_rgrs = tuple(
        rgr for rgr in rgrs if maximum_spin_multiplicity(rgr) == mult_min)
    return dom_rgrs


def resonances(rgr):
    """ all resonance graphs with this connectivity
    """
    return subresonances(without_bond_orders(rgr))


def subresonances(rgr):
    """ this connected graph and its lower-spin (more pi-bonded) resonances
    """
    rgr = without_fractional_bonds(rgr)
    # get the bond capacities (room for increasing bond order), filtering out
    # the negative ones to avoid complications with hypervalent atoms in TSs
    bnd_cap_dct = dict_.by_value(_bond_capacities(rgr), lambda x: x > 0)

    ret_rgrs = []
    if bnd_cap_dct:
        bnd_keys, bnd_caps = zip(*bnd_cap_dct.items())
        atm_keys = list(functools.reduce(frozenset.union, bnd_keys))

        # Loop over all possible combinations of bond order increments (amounts
        # by which to increase the bond order), filtering out combinations that
        # exceed the valences of the atoms involved.
        # (Note that we are only testing the bonds with available pi electrons,
        # so this is compatible with having hypervalent atoms elsewhere in the
        # molecule)
        bnd_ord_inc_ranges = [range(bnd_cap+1) for bnd_cap in bnd_caps]
        for bnd_ord_incs in itertools.product(*bnd_ord_inc_ranges):
            bnd_ord_inc_dct = dict(zip(bnd_keys, bnd_ord_incs))
            ret_rgr = _add_pi_bonds(rgr, bnd_ord_inc_dct)

            max_bnd_ord = max(bond_orders(ret_rgr).values())

            atm_unsat_vlcs = dict_.values_by_key(
                atom_unsaturated_valences(ret_rgr), atm_keys)

            if not any(atm_unsat_vlc < 0 for atm_unsat_vlc in atm_unsat_vlcs):
                if max_bnd_ord < 4:
                    ret_rgrs.append(ret_rgr)

    if not ret_rgrs:
        ret_rgrs = (rgr,)
    else:
        ret_rgrs = tuple(ret_rgrs)

    return ret_rgrs


def _bond_capacities(rgr):
    """ the number of electron pairs available for further pi-bonding, by bond
    """
    rgr = without_dummy_bonds(rgr)
    atm_unsat_vlc_dct = atom_unsaturated_valences(rgr)

    def _pi_capacities(bnd_key):
        return min(map(atm_unsat_vlc_dct.__getitem__, bnd_key))

    bnd_keys = list(bond_keys(rgr))
    bnd_caps = tuple(map(_pi_capacities, bnd_keys))
    bnd_cap_dct = dict(zip(bnd_keys, bnd_caps))
    return bnd_cap_dct


def _add_pi_bonds(rgr, bnd_ord_inc_dct):
    """ add pi bonds to this graph
    """
    bnd_keys = bond_keys(rgr)
    assert set(bnd_ord_inc_dct.keys()) <= bnd_keys

    bnd_keys = list(bnd_keys)
    bnd_ords = dict_.values_by_key(bond_orders(rgr), bnd_keys)
    bnd_ord_incs = dict_.values_by_key(bnd_ord_inc_dct, bnd_keys, fill_val=0)
    new_bnd_ords = numpy.add(bnd_ords, bnd_ord_incs)
    bnd_ord_dct = dict(zip(bnd_keys, new_bnd_ords))
    rgr = set_bond_orders(rgr, bnd_ord_dct)
    return rgr


def bond_neighborhoods(gra, stereo=False):
    """ neighborhood subgraphs, by bond
    """
    bnd_keys = list(bond_keys(gra))

    def _neighborhood(bnd_key):
        nbh_bnd_keys = set(filter(lambda x: bnd_key & x, bnd_keys))
        return bond_induced_subgraph(gra, nbh_bnd_keys, stereo=stereo)

    bnd_nbh_dct = dict(zip(bnd_keys, map(_neighborhood, bnd_keys)))
    return bnd_nbh_dct


def sp2_bond_keys(gra):
    """ determine the sp2 bonds in this graph
    """
    gra = without_bond_orders(gra)
    bnd_keys = dict_.keys_by_value(
        resonance_dominant_bond_orders(gra), lambda x: 2 in x)

    # make sure both ends are sp^2 (excludes cumulenes)
    atm_hyb_dct = resonance_dominant_atom_hybridizations(gra)
    sp2_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 2)
    bnd_keys = frozenset({bnd_key for bnd_key in bnd_keys
                          if bnd_key <= sp2_atm_keys})
    return bnd_keys


def resonance_dominant_bond_orders(rgr):
    """ resonance-dominant bond orders, by bond
    """
    rgr = without_fractional_bonds(rgr)
    bnd_keys = list(bond_keys(rgr))
    bnd_ords_by_res = [
        dict_.values_by_key(bond_orders(dom_rgr), bnd_keys)
        for dom_rgr in dominant_resonances(rgr)]
    bnd_ords_lst = list(map(frozenset, zip(*bnd_ords_by_res)))
    bnd_dom_res_ords_dct = dict(zip(bnd_keys, bnd_ords_lst))
    return bnd_dom_res_ords_dct


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


def resonance_dominant_atom_hybridizations(rgr):
    """ resonance-dominant atom hybridizations, by atom
    """
    rgr = without_fractional_bonds(rgr)
    atm_keys = list(atom_keys(rgr))
    atm_hybs_by_res = [
        dict_.values_by_key(atom_hybridizations(dom_rgr), atm_keys)
        for dom_rgr in dominant_resonances(rgr)]
    atm_hybs = [min(hybs) for hybs in zip(*atm_hybs_by_res)]
    atm_hyb_dct = dict(zip(atm_keys, atm_hybs))
    return atm_hyb_dct


def atom_lone_pair_counts(gra):
    """ lone pair counts, by atom
    """
    atm_symb_dct = atom_symbols(gra)
    atm_group_idx_dct = dict_.transform_values(atm_symb_dct, ptab.to_group)
    atm_lpc_dct = dict_.transform_values(atm_group_idx_dct,
                                         LONE_PAIR_COUNTS_DCT.__getitem__)
    atm_lpc_dct = dict_.transform_values(atm_lpc_dct, int)
    return atm_lpc_dct


def atom_stereo_sorted_neighbor_atom_keys(gra, atm_key, atm_ngb_keys):
    """ get the neighbor keys of an atom sorted by stereo priority
    """
    atm_ngb_keys = list(atm_ngb_keys)

    # explicitly create an object array because otherwise the argsort
    # interprets [()] as []
    atm_pri_vecs = numpy.empty(len(atm_ngb_keys), dtype=numpy.object_)
    atm_pri_vecs[:] = [stereo_priority_vector(gra, atm_key, atm_ngb_key)
                       for atm_ngb_key in atm_ngb_keys]

    sort_idxs = numpy.argsort(atm_pri_vecs)
    sorted_atm_ngb_keys = tuple(map(atm_ngb_keys.__getitem__, sort_idxs))
    return sorted_atm_ngb_keys


def stereo_priority_vector(gra, atm_key, atm_ngb_key):
    """ generates a sortable one-to-one representation of the branch extending
    from `atm_key` through its bonded neighbor `atm_ngb_key`
    """
    bbn_keys = backbone_keys(gra)
    exp_hyd_keys = explicit_hydrogen_keys(gra)

    if atm_ngb_key not in bbn_keys:
        assert atm_ngb_key in exp_hyd_keys
        assert frozenset({atm_key, atm_ngb_key}) in bonds(gra)
        pri_vec = ()
    else:
        gra = implicit(gra)
        atm_dct = atoms(gra)
        bnd_dct = bonds(gra)
        assert atm_key in bbn_keys
        assert frozenset({atm_key, atm_ngb_key}) in bnd_dct

        # here, switch to an implicit graph
        atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

        def _priority_vector(atm1_key, atm2_key, seen_keys):
            # we keep a list of seen keys to cut off cycles, avoiding infinite
            # loops

            bnd_val = bnd_dct[frozenset({atm1_key, atm2_key})]
            atm_val = atm_dct[atm2_key]

            bnd_val = _replace_nones_with_negative_infinity(bnd_val)
            atm_val = _replace_nones_with_negative_infinity(atm_val)

            if atm2_key in seen_keys:
                ret = (bnd_val,)
            else:
                seen_keys.update({atm1_key, atm2_key})
                atm3_keys = atm_ngb_keys_dct[atm2_key] - {atm1_key}
                if atm3_keys:
                    next_vals, seen_keys = zip(*[
                        _priority_vector(atm2_key, atm3_key, seen_keys)
                        for atm3_key in atm3_keys])
                    ret = (bnd_val, atm_val) + next_vals
                else:
                    ret = (bnd_val, atm_val)

            return ret, seen_keys

        pri_vec, _ = _priority_vector(atm_key, atm_ngb_key, set())

    return pri_vec


def _replace_nones_with_negative_infinity(seq):
    return [-numpy.inf if val is None else val for val in seq]


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


ATM_PROP_NAMES = ('symbol', 'implicit_hydrogen_valence', 'stereo_parity')
BND_PROP_NAMES = ('order', 'stereo_parity')


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
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


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
        srt = argsort_symbols(
            srt_vals, symbs_first, symbs_last, idx=2)
        keys = tuple(map(keys.__getitem__, srt))
        return keys

    atm_ngb_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(gra), _neighbor_keys)
    return atm_ngb_keys_dct


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
