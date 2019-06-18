""" molecular graph
"""
import itertools
import functools
import numpy
import future.moves.itertools as fmit
from qcelemental import periodictable as pt
from automol import dict_
from automol.graph import _networkx
import automol.dict_.multi as mdict
import automol.create.graph as _create

ATM_SYM_POS = 0
ATM_IMP_HYD_VLC_POS = 1
ATM_STE_PAR_POS = 2

BND_ORD_POS = 0
BND_STE_PAR_POS = 1


# getters
def atoms(xgr):
    """ atoms, as a dictionary
    """
    atm_dct, _ = xgr
    return atm_dct


def bonds(xgr):
    """ bonds, as a dictionary
    """
    _, bnd_dct = xgr
    return bnd_dct


def atom_keys(xgr):
    """ sorted atom keys
    """
    return frozenset(atoms(xgr).keys())


def bond_keys(xgr):
    """ sorted bond keys
    """
    return frozenset(bonds(xgr).keys())


def atom_symbols(xgr):
    """ atom symbols, as a dictionary
    """
    return mdict.by_key_by_position(atoms(xgr), atom_keys(xgr), ATM_SYM_POS)


def atom_implicit_hydrogen_valences(xgr):
    """ atom implicit hydrogen valences, as a dictionary
    """
    return mdict.by_key_by_position(atoms(xgr), atom_keys(xgr),
                                    ATM_IMP_HYD_VLC_POS)


def atom_stereo_parities(sgr):
    """ atom parities, as a dictionary
    """
    return mdict.by_key_by_position(atoms(sgr), atom_keys(sgr),
                                    ATM_STE_PAR_POS)


def bond_orders(rgr):
    """ bond orders, as a dictionary
    """
    return mdict.by_key_by_position(bonds(rgr), bond_keys(rgr), BND_ORD_POS)


def bond_stereo_parities(sgr):
    """ bond parities, as a dictionary
    """
    return mdict.by_key_by_position(bonds(sgr), bond_keys(sgr),
                                    BND_STE_PAR_POS)


# setters
def relabel(xgr, atm_key_dct):
    """ relabel the graph with new atom keys
    """
    orig_atm_keys = atom_keys(xgr)
    assert set(atm_key_dct.keys()) <= orig_atm_keys

    new_atm_key_dct = dict(zip(orig_atm_keys, orig_atm_keys))
    new_atm_key_dct.update(atm_key_dct)

    _relabel_atom_key = new_atm_key_dct.__getitem__

    def _relabel_bond_key(bnd_key):
        return frozenset(map(_relabel_atom_key, bnd_key))

    atm_dct = dict_.transform_keys(atoms(xgr), _relabel_atom_key)
    bnd_dct = dict_.transform_keys(bonds(xgr), _relabel_bond_key)
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


def standard_keys(xgr):
    """ replace the current atom keys with standard indices, counting from zero
    """
    atm_key_dct = dict(enumerate(sorted(atom_keys(xgr))))
    return relabel(xgr, atm_key_dct)


def transform_keys(xgr, atm_key_func):
    """ transform atom keys with a function
    """
    atm_keys = atom_keys(xgr)
    atm_key_dct = dict(zip(atm_keys, map(atm_key_func, atm_keys)))
    return relabel(xgr, atm_key_dct)


def set_atom_implicit_hydrogen_valences(xgr, atm_imp_hyd_vlc_dct):
    """ set atom implicit hydrogen valences
    """
    atm_dct = mdict.set_by_key_by_position(atoms(xgr), atm_imp_hyd_vlc_dct,
                                           ATM_IMP_HYD_VLC_POS)
    bnd_dct = bonds(xgr)
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


def set_atom_stereo_parities(sgr, atm_par_dct):
    """ set atom parities
    """
    atm_dct = mdict.set_by_key_by_position(atoms(sgr), atm_par_dct,
                                           ATM_STE_PAR_POS)
    return _create.from_atoms_and_bonds(atm_dct, bonds(sgr))


def set_bond_orders(rgr, bnd_ord_dct):
    """ set bond orders
    """
    bnd_dct = mdict.set_by_key_by_position(bonds(rgr), bnd_ord_dct,
                                           BND_ORD_POS)
    return _create.from_atoms_and_bonds(atoms(rgr), bnd_dct)


def set_bond_stereo_parities(sgr, bnd_par_dct):
    """ set bond parities
    """
    bnd_dct = mdict.set_by_key_by_position(bonds(sgr), bnd_par_dct,
                                           BND_STE_PAR_POS)
    return _create.from_atoms_and_bonds(atoms(sgr), bnd_dct)


def add_atom_implicit_hydrogen_valences(xgr, inc_atm_imp_hyd_vlc_dct):
    """ add atom imlicit hydrogen valences

    (increments can be positive or negative)
    """
    atm_keys = list(inc_atm_imp_hyd_vlc_dct.keys())
    atm_imp_hyd_vlcs = numpy.add(
        dict_.values_by_key(atom_implicit_hydrogen_valences(xgr), atm_keys),
        dict_.values_by_key(inc_atm_imp_hyd_vlc_dct, atm_keys))
    assert all(atm_imp_hyd_vlc >= 0 for atm_imp_hyd_vlc in atm_imp_hyd_vlcs)
    atm_imp_hyd_vlc_dct = dict(zip(atm_keys, atm_imp_hyd_vlcs))
    return set_atom_implicit_hydrogen_valences(xgr, atm_imp_hyd_vlc_dct)


def without_bond_orders(xgr):
    """ resonance graph with maximum spin (i.e. no pi bonds)
    """
    bnd_ord_dct = dict_.by_key({}, bond_keys(xgr), fill_val=1)
    return set_bond_orders(xgr, bnd_ord_dct)


def without_stereo_parities(xgr):
    """ graph with stereo assignments wiped out
    """
    atm_ste_par_dct = dict_.by_key({}, atom_keys(xgr), fill_val=None)
    bnd_ste_par_dct = dict_.by_key({}, bond_keys(xgr), fill_val=None)
    xgr = set_atom_stereo_parities(xgr, atm_ste_par_dct)
    xgr = set_bond_stereo_parities(xgr, bnd_ste_par_dct)
    return xgr


def add_atoms(xgr, sym_dct, imp_hyd_vlc_dct=None, ste_par_dct=None):
    """ add atoms to this molecular graph
    """
    atm_keys = atom_keys(xgr)
    atm_sym_dct = atom_symbols(xgr)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(xgr)
    atm_ste_par_dct = atom_stereo_parities(xgr)

    keys = set(sym_dct.keys())
    imp_hyd_vlc_dct = {} if imp_hyd_vlc_dct is None else imp_hyd_vlc_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct

    assert not keys & atm_keys
    assert set(imp_hyd_vlc_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys

    atm_sym_dct.update(sym_dct)
    atm_imp_hyd_vlc_dct.update(imp_hyd_vlc_dct)
    atm_ste_par_dct.update(ste_par_dct)

    atm_dct = _create.atoms_from_data(
        atom_symbols=atm_sym_dct,
        atom_implicit_hydrogen_valences=atm_imp_hyd_vlc_dct,
        atom_stereo_parities=atm_ste_par_dct)
    bnd_dct = bonds(xgr)
    xgr = _create.from_atoms_and_bonds(atoms=atm_dct, bonds=bnd_dct)
    return xgr


def add_bonds(xgr, keys, ord_dct=None, ste_par_dct=None):
    """ add bonds to this molecular graph
    """
    bnd_keys = set(bond_keys(xgr))
    bnd_ord_dct = bond_orders(xgr)
    bnd_ste_par_dct = bond_stereo_parities(xgr)

    keys = set(map(frozenset, keys))
    ord_dct = {} if ord_dct is None else ord_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct

    assert not keys & bnd_keys
    assert set(ord_dct.keys()) <= keys
    assert set(ste_par_dct.keys()) <= keys

    bnd_keys.update(keys)
    bnd_ord_dct.update(ord_dct)
    bnd_ste_par_dct.update(ste_par_dct)

    atm_dct = atoms(xgr)
    bnd_dct = _create.bonds_from_data(
        bond_keys=bnd_keys, bond_orders=bnd_ord_dct,
        bond_stereo_parities=bnd_ste_par_dct)

    xgr = _create.from_atoms_and_bonds(atoms=atm_dct, bonds=bnd_dct)
    return xgr


def frozen(xgr):
    """ hashable, sortable, immutable container of graph data
    """
    atm_keys = sorted(atom_keys(xgr))
    bnd_keys = sorted(bond_keys(xgr), key=sorted)

    # make it sortable by replacing Nones with -infinity
    atm_vals = numpy.array(dict_.values_by_key(atoms(xgr), atm_keys))
    bnd_vals = numpy.array(dict_.values_by_key(bonds(xgr), bnd_keys))
    atm_vals[numpy.equal(atm_vals, None)] = -numpy.inf
    bnd_vals[numpy.equal(bnd_vals, None)] = -numpy.inf

    frz_atms = tuple(zip(atm_keys, map(tuple, atm_vals)))
    frz_bnds = tuple(zip(bnd_keys, map(tuple, bnd_vals)))
    return (frz_atms, frz_bnds)


# graph theory library
# # atom properties
def atom_neighbor_keys(xgr):
    """ keys of neighboring atoms, by atom
    """
    def _neighbor_keys(atm_key, atm_nbh):
        return frozenset(atom_keys(atm_nbh) - {atm_key})

    atm_ngb_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(xgr), _neighbor_keys)
    return atm_ngb_keys_dct


def atom_bond_keys(xgr):
    """ bond keys, by atom
    """
    return dict_.transform_values(atom_neighborhoods(xgr), bond_keys)


def atom_neighborhoods(xgr):
    """ neighborhood subgraphs, by atom
    """
    bnd_keys = bond_keys(xgr)

    def _neighborhood(atm_key):
        nbh_bnd_keys = set(filter(lambda x: atm_key in x, bnd_keys))
        return bond_induced_subgraph(xgr, nbh_bnd_keys)

    atm_keys = list(atom_keys(xgr))
    atm_nbh_dct = dict(zip(atm_keys, map(_neighborhood, atm_keys)))
    return atm_nbh_dct


# # bond properties
def bond_neighbor_keys(xgr):
    """ keys of neighboring bonds, by bond
    """
    def _neighbor_keys(bnd_key, bnd_nbh):
        return frozenset(bond_keys(bnd_nbh) - {bnd_key})

    bnd_ngb_keys_dct = dict_.transform_items_to_values(
        bond_neighborhoods(xgr), _neighbor_keys)
    return bnd_ngb_keys_dct


def bond_neighborhoods(xgr):
    """ neighborhood subgraphs, by bond
    """
    bnd_keys = list(bond_keys(xgr))

    def _neighborhood(bnd_key):
        nbh_bnd_keys = set(filter(lambda x: bnd_key & x, bnd_keys))
        return bond_induced_subgraph(xgr, nbh_bnd_keys)

    bnd_nbh_dct = dict(zip(bnd_keys, map(_neighborhood, bnd_keys)))
    return bnd_nbh_dct


# # other properties
def branch(xgr, atm_key, bnd_key):
    """ branch extending along `bnd_key` away from `atm_key`
    """
    return bond_induced_subgraph(xgr, branch_bond_keys(xgr, atm_key, bnd_key))


def branch_bond_keys(xgr, atm_key, bnd_key):
    """ keys for branch extending along `bnd_key` away from `atm_key`
    """
    bnd_key = frozenset(bnd_key)
    assert atm_key in bnd_key
    assert bnd_key in bond_keys(xgr)

    atm_bnd_keys_dct = atom_bond_keys(xgr)

    bnch_bnd_keys = {bnd_key}
    seen_bnd_keys = set()
    excl_bnd_keys = atm_bnd_keys_dct[atm_key] - {bnd_key}

    new_bnd_keys = {bnd_key}

    bnd_ngb_keys_dct = bond_neighbor_keys(xgr)
    while new_bnd_keys:
        new_bnd_ngb_keys = set(
            itertools.chain(
                *dict_.values_by_key(bnd_ngb_keys_dct, new_bnd_keys)))
        bnch_bnd_keys.update(new_bnd_ngb_keys - excl_bnd_keys)
        seen_bnd_keys.update(new_bnd_keys)
        new_bnd_keys = bnch_bnd_keys - seen_bnd_keys

    return frozenset(bnch_bnd_keys)


def rings(xgr):
    """ rings in the graph (minimal basis)
    """
    xgrs = [bond_induced_subgraph(xgr, bnd_keys)
            for bnd_keys in rings_bond_keys(xgr)]
    return tuple(sorted(xgrs, key=frozen))


def rings_sorted_atom_keys(xgr):
    """ atom keys for each ring in the graph sorted by connectivity (minimal basis)
    """
    def _sorted_ring_atom_keys(rng_bnd_keys):
        rng_bnd_keys = list(rng_bnd_keys)
        bnd_key = min(rng_bnd_keys, key=sorted)
        first_atm_key, atm_key = sorted(bnd_key)
        rng_bnd_keys.remove(bnd_key)
        rng_atm_keys = [first_atm_key, atm_key]
        while rng_bnd_keys:
            bnd_key = next(filter(lambda x: atm_key in x, rng_bnd_keys))
            rng_bnd_keys.remove(bnd_key)
            bnd_key = set(bnd_key)
            bnd_key.remove(atm_key)
            atm_key = next(iter(bnd_key))
            rng_atm_keys.append(atm_key)
        rng_atm_keys.pop(-1)
        rng_atm_keys = tuple(rng_atm_keys)
        return rng_atm_keys

    rng_atm_keys_lst = frozenset(
        map(_sorted_ring_atom_keys, rings_bond_keys(xgr)))
    return rng_atm_keys_lst


def rings_bond_keys(xgr):
    """ bond keys for each ring in the graph (minimal basis)
    """
    bnd_keys = bond_keys(xgr)

    def _ring_bond_keys(rng_atm_keys):
        return frozenset(filter(lambda x: x <= rng_atm_keys, bnd_keys))

    nxg = _networkx.from_graph(xgr)
    rng_atm_keys_lst = _networkx.minimum_cycle_basis(nxg)
    rng_bnd_keys_lst = frozenset(map(_ring_bond_keys, rng_atm_keys_lst))
    return rng_bnd_keys_lst


def connected_components(xgr):
    """ connected components in the graph
    """
    cmp_xgr_atm_keys_lst = connected_components_atom_keys(xgr)
    cmp_xgrs = tuple(subgraph(xgr, cmp_xgr_atm_keys)
                     for cmp_xgr_atm_keys in cmp_xgr_atm_keys_lst)
    return cmp_xgrs


def connected_components_atom_keys(xgr):
    """ atom keys for each connected component in the graph
    """
    nxg = _networkx.from_graph(xgr)
    cmp_xgr_atm_keys_lst = _networkx.connected_component_atom_keys(nxg)
    return cmp_xgr_atm_keys_lst


def union(xgr1, xgr2):
    """ a union of two graphs
    """
    assert not atom_keys(xgr1) & atom_keys(xgr2)
    atm_dct = {}
    atm_dct.update(atoms(xgr1))
    atm_dct.update(atoms(xgr2))

    bnd_dct = {}
    bnd_dct.update(bonds(xgr1))
    bnd_dct.update(bonds(xgr2))
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


def subgraph(xgr, atm_keys):
    """ the subgraph induced by a subset of the atoms
    """
    atm_keys = set(atm_keys)
    assert atm_keys <= atom_keys(xgr)
    bnd_keys = set(filter(lambda x: x <= atm_keys, bond_keys(xgr)))
    atm_dct = dict_.by_key(atoms(xgr), atm_keys)
    bnd_dct = dict_.by_key(bonds(xgr), bnd_keys)
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


def bond_induced_subgraph(xgr, bnd_keys):
    """ the subgraph induced by a subset of the bonds
    """
    atm_keys = set(itertools.chain(*bnd_keys))
    bnd_keys = set(bnd_keys)
    assert atm_keys <= atom_keys(xgr)
    assert bnd_keys <= bond_keys(xgr)
    atm_dct = dict_.by_key(atoms(xgr), atm_keys)
    bnd_dct = dict_.by_key(bonds(xgr), bnd_keys)
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


# # transformations
def remove_atoms(xgr, atm_keys):
    """ remove atoms from the molecular graph
    """
    all_atm_keys = atom_keys(xgr)
    atm_keys = set(atm_keys)
    assert atm_keys <= all_atm_keys
    atm_keys_left = all_atm_keys - atm_keys
    return subgraph(xgr, atm_keys_left)


def remove_bonds(xgr, bnd_keys):
    """ remove bonds from the molecular graph
    """
    all_bnd_keys = bond_keys(xgr)
    bnd_keys = set(bnd_keys)
    assert bnd_keys <= all_bnd_keys
    bnd_keys = all_bnd_keys - bnd_keys
    atm_dct = atoms(xgr)
    bnd_dct = dict_.by_key(bonds(xgr), bnd_keys)
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


def without_ghost_atoms(xgr):
    """ remove ghost atoms from the graph
    """
    atm_sym_dct = atom_symbols(xgr)
    atm_keys = [key for key, sym in atm_sym_dct.items() if pt.to_Z(sym)]
    return subgraph(xgr, atm_keys)


# implicit/explicit hydrogen functions
# # atom properties
def atom_explicit_hydrogen_valences(xgr):
    """ explicit hydrogen valences, by atom
    """
    return dict_.transform_values(atom_explicit_hydrogen_keys(xgr), len)


def atom_explicit_hydrogen_keys(xgr):
    """ explicit hydrogen valences, by atom
    """
    def _explicit_hydrogen_keys(atm_key, atm_nbh):
        return frozenset(explicit_hydrogen_keys(atm_nbh) - {atm_key})

    atm_exp_hyd_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(xgr), _explicit_hydrogen_keys)
    return atm_exp_hyd_keys_dct


# # other properties
def backbone_keys(xgr):
    """ backbone atom keys
    """
    bbn_keys = atom_keys(xgr) - explicit_hydrogen_keys(xgr)
    return bbn_keys


def explicit_hydrogen_keys(xgr):
    """ explicit hydrogen keys (H types: explicit, implicit, backbone)
    """
    hyd_keys = dict_.keys_by_value(atom_symbols(xgr), lambda x: x == 'H')
    atm_ngb_keys_dct = atom_neighbor_keys(xgr)

    def _is_backbone(hyd_key):
        return all(ngb_key in hyd_keys and hyd_key < ngb_key
                   for ngb_key in atm_ngb_keys_dct[hyd_key])

    exp_hyd_keys = frozenset(fmit.filterfalse(_is_backbone, hyd_keys))
    return exp_hyd_keys


# # transformations
def add_atom_explicit_hydrogen_keys(xgr, atm_exp_hyd_keys_dct):
    """ add explicit hydrogens by atom
    """
    assert set(atm_exp_hyd_keys_dct.keys()) <= atom_keys(xgr)
    for atm_key, atm_exp_hyd_keys in atm_exp_hyd_keys_dct.items():
        assert not set(atm_exp_hyd_keys) & atom_keys(xgr)
        atm_exp_hyd_bnd_keys = {frozenset({atm_key, atm_exp_hyd_key})
                                for atm_exp_hyd_key in atm_exp_hyd_keys}
        atm_exp_hyd_sym_dct = dict_.by_key({}, atm_exp_hyd_keys, fill_val='H')
        xgr = add_atoms(xgr, atm_exp_hyd_sym_dct)
        xgr = add_bonds(xgr, atm_exp_hyd_bnd_keys)
    return xgr


def implicit(xgr, atm_keys=None):
    """ make the hydrogens at these atoms implicit
    """
    atm_keys = backbone_keys(xgr) if atm_keys is None else atm_keys

    atm_exp_hyd_keys_dct = dict_.by_key(
        atom_explicit_hydrogen_keys(xgr), atm_keys)

    inc_imp_hyd_keys_dct = dict_.transform_values(atm_exp_hyd_keys_dct, len)
    xgr = add_atom_implicit_hydrogen_valences(xgr, inc_imp_hyd_keys_dct)

    exp_hyd_keys = set(itertools.chain(*atm_exp_hyd_keys_dct.values()))
    xgr = remove_atoms(xgr, exp_hyd_keys)
    return xgr


def explicit(xgr, atm_keys=None):
    """ make the hydrogens at these atoms explicit
    """
    atm_keys = backbone_keys(xgr) if atm_keys is None else atm_keys
    atm_keys = sorted(atm_keys)
    atm_imp_hyd_vlc_dct = dict_.by_key(
        atom_implicit_hydrogen_valences(xgr), atm_keys)

    atm_exp_hyd_keys_dct = {}
    next_atm_key = max(atom_keys(xgr)) + 1
    for atm_key in atm_keys:
        imp_hyd_vlc = atm_imp_hyd_vlc_dct[atm_key]
        atm_exp_hyd_keys_dct[atm_key] = set(
            range(next_atm_key, next_atm_key+imp_hyd_vlc))
        next_atm_key += imp_hyd_vlc

    xgr = set_atom_implicit_hydrogen_valences(
        xgr, dict_.by_key({}, atm_keys, fill_val=0))
    xgr = add_atom_explicit_hydrogen_keys(xgr, atm_exp_hyd_keys_dct)
    return xgr


# # comparisons
def full_isomorphism(xgr1, xgr2):
    """ full graph isomorphism
    """
    assert xgr1 == explicit(xgr1) and xgr2 == explicit(xgr2)
    nxg1 = _networkx.from_graph(xgr1)
    nxg2 = _networkx.from_graph(xgr2)
    iso_dct = _networkx.isomorphism(nxg1, nxg2)
    return iso_dct


def backbone_isomorphic(xgr1, xgr2):
    """ are these molecular graphs backbone isomorphic?
    """
    return backbone_isomorphism(xgr1, xgr2) is not None


def backbone_isomorphism(xgr1, xgr2):
    """ graph backbone isomorphism

    for implicit graphs, this is the relabeling of `xgr1` to produce `xgr2`
    for other graphs, it gives the correspondences between backbone atoms
    """
    xgr1 = implicit(xgr1)
    xgr2 = implicit(xgr2)
    nxg1 = _networkx.from_graph(xgr1)
    nxg2 = _networkx.from_graph(xgr2)
    iso_dct = _networkx.isomorphism(nxg1, nxg2)
    return iso_dct


def backbone_unique(xgrs):
    """ unique non-isomorphic graphs from a series
    """
    xgrs = _unique(xgrs, equiv=backbone_isomorphic)
    return xgrs


def _unique(itms, equiv):
    """ unique items from a list, according to binary comparison `equiv`
    """
    uniq_itms = []
    for itm in itms:
        if not any(map(functools.partial(equiv, itm), uniq_itms)):
            uniq_itms.append(itm)

    return tuple(uniq_itms)


# chemistry library
VALENCE_DCT = {
    None: None,
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
    None: None,
    1: 0,   # H
    2: 0,   # Be
    13: 0,  # B
    14: 0,  # C
    15: 1,  # N
    16: 2,  # O
    17: 3,  # F
    18: 4,  # He
}


# # atom properties
def atom_element_valences(xgr):
    """ element valences (# possible single bonds), by atom
    """
    atm_sym_dct = atom_symbols(xgr)
    atm_group_idx_dct = dict_.transform_values(atm_sym_dct, pt.to_group)
    atm_elem_vlc_dct = dict_.transform_values(atm_group_idx_dct,
                                              VALENCE_DCT.__getitem__)
    return atm_elem_vlc_dct


def atom_lone_pair_counts(xgr):
    """ lone pair counts, by atom
    """
    atm_sym_dct = atom_symbols(xgr)
    atm_group_idx_dct = dict_.transform_values(atm_sym_dct, pt.to_group)
    atm_lpc_dct = dict_.transform_values(atm_group_idx_dct,
                                         LONE_PAIR_COUNTS_DCT.__getitem__)
    return atm_lpc_dct


def atom_bond_valences(xgr, bond_order=True):
    """ bond count (bond valence), by atom
    """
    atm_keys = list(atom_keys(xgr))
    xgr = explicit(xgr)
    if not bond_order:
        xgr = without_bond_orders(xgr)

    atm_nbhs = dict_.values_by_key(atom_neighborhoods(xgr), atm_keys)
    atm_bnd_vlcs = [sum(bond_orders(nbh).values()) for nbh in atm_nbhs]
    atm_bnd_vlc_dct = dict(zip(atm_keys, atm_bnd_vlcs))
    return atm_bnd_vlc_dct


def atom_unsaturated_valences(xgr, bond_order=True):
    """ unsaturated valences, by atom

    (element valences minus bonding valences -- pi sites and radical electrons)
    """
    atm_keys = list(atom_keys(xgr))
    if not bond_order:
        xgr = without_bond_orders(xgr)

    atm_bnd_vlcs = dict_.values_by_key(atom_bond_valences(xgr), atm_keys)
    atm_tot_vlcs = dict_.values_by_key(atom_element_valences(xgr), atm_keys)
    atm_rad_vlcs = numpy.subtract(atm_tot_vlcs, atm_bnd_vlcs)
    return dict(zip(atm_keys, atm_rad_vlcs))


def unsaturated_atom_keys(xgr):
    """ keys of unsaturated (radical or pi-bonded) atoms
    """
    atm_unsat_vlc_dct = atom_unsaturated_valences(xgr, bond_order=False)
    unsat_atm_keys = frozenset(dict_.keys_by_value(atm_unsat_vlc_dct, bool))
    return unsat_atm_keys


# # other properties
def maximum_spin_multiplicity(xgr, bond_order=True):
    """ the highest possible spin multiplicity for this molecular graph
    """
    atm_rad_vlc_dct = atom_unsaturated_valences(xgr, bond_order=bond_order)
    return sum(atm_rad_vlc_dct.values()) + 1


def possible_spin_multiplicities(xgr, bond_order=True):
    """ possible spin multiplicities for this molecular graph
    """
    mult_max = maximum_spin_multiplicity(xgr, bond_order=bond_order)
    mult_min = 2 if mult_max % 2 == 0 else 1
    mults = tuple(range(mult_min, mult_max+1, 2))
    return mults


# miscellaneous
# # bond properties
def bond_symmetry_numbers(xgr):
    """ symmetry numbers, by bond

    the (approximate) symmetry number of the torsional potential for this bond,
    based on the hydrogen counts for each atom
    """
    imp_xgr = implicit(xgr)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(imp_xgr)

    bnd_keys = bond_keys(imp_xgr)
    bnd_max_hyd_vlcs = [max(map(atm_imp_hyd_vlc_dct.__getitem__, bnd_key))
                        for bnd_key in bnd_keys]
    bnd_sym_nums = [3 if vlc == 3 else 1 for vlc in bnd_max_hyd_vlcs]
    bnd_sym_num_dct = dict(zip(bnd_keys, bnd_sym_nums))

    # fill in the rest of the bonds for completeness
    bnd_sym_num_dct = dict_.by_key(bnd_sym_num_dct, bond_keys(xgr), fill_val=1)
    return bnd_sym_num_dct
