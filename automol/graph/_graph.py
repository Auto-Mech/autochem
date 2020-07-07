""" molecular graph
"""
import itertools
import functools
import numpy
import future.moves.itertools as fmit
from qcelemental import periodictable as pt
from automol import dict_
from automol.graph._graph_base import atoms
from automol.graph._graph_base import bonds
from automol.graph._graph_base import atom_keys
from automol.graph._graph_base import bond_keys
from automol.graph._graph_base import atom_symbols
from automol.graph._graph_base import atom_implicit_hydrogen_valences
from automol.graph._graph_base import atom_stereo_parities
from automol.graph._graph_base import bond_orders
from automol.graph._graph_base import bond_stereo_parities
from automol.graph._graph_base import set_atom_implicit_hydrogen_valences
from automol.graph._graph_base import set_atom_stereo_parities
from automol.graph._graph_base import set_bond_orders
from automol.graph._graph_base import set_bond_stereo_parities
from automol.graph._graph_base import relabel
from automol.graph import _networkx
import automol.create.graph as _create


# setters
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
        natms = atom_count(gra, with_dummy=True, with_implicit=False)

        atm_key_dct = {atm_key: idx+shift
                       for idx, atm_key in enumerate(sorted(atom_keys(gra)))}
        atm_key_dcts.append(atm_key_dct)

        shift += natms

    gras = [relabel(gra, atm_key_dct)
            for gra, atm_key_dct in zip(gras, atm_key_dcts)]

    return gras, atm_key_dcts


def transform_keys(gra, atm_key_func):
    """ transform atom keys with a function
    """
    atm_keys = atom_keys(gra)
    atm_key_dct = dict(zip(atm_keys, map(atm_key_func, atm_keys)))
    return relabel(gra, atm_key_dct)


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


def without_stereo_parities(gra):
    """ graph with stereo assignments wiped out
    """
    atm_ste_par_dct = dict_.by_key({}, atom_keys(gra), fill_val=None)
    bnd_ste_par_dct = dict_.by_key({}, bond_keys(gra), fill_val=None)
    gra = set_atom_stereo_parities(gra, atm_ste_par_dct)
    gra = set_bond_stereo_parities(gra, bnd_ste_par_dct)
    return gra


def frozen(gra):
    """ hashable, sortable, immutable container of graph data
    """
    atm_keys = sorted(atom_keys(gra))
    bnd_keys = sorted(bond_keys(gra), key=sorted)

    # make it sortable by replacing Nones with -infinity
    atm_vals = numpy.array(dict_.values_by_key(atoms(gra), atm_keys))
    bnd_vals = numpy.array(dict_.values_by_key(bonds(gra), bnd_keys))
    atm_vals[numpy.equal(atm_vals, None)] = -numpy.inf
    bnd_vals[numpy.equal(bnd_vals, None)] = -numpy.inf

    frz_atms = tuple(zip(atm_keys, map(tuple, atm_vals)))
    frz_bnds = tuple(zip(bnd_keys, map(tuple, bnd_vals)))
    return (frz_atms, frz_bnds)


# graph theory library
# # atom properties
def electron_count(gra, charge=0):
    """ the number of electrons in the molecule
    """
    atm_sym_dct = atom_symbols(explicit(gra))
    nelec = sum(map(pt.to_Z, atm_sym_dct.values())) - charge
    return nelec


def atom_count(gra, with_dummy=False, with_implicit=True):
    """ count the number of atoms in this molecule

    by default, this includes implicit hydrogens and excludes dummy atoms
    """
    if not with_dummy:
        gra = without_dummy_atoms(gra)
    natms = len(atoms(gra))
    if with_implicit:
        atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(gra)
        natms += sum(atm_imp_hyd_vlc_dct.values())
    return natms


def heavy_atom_count(gra, with_dummy=False):
    """ the number of heavy atoms
    """
    if not with_dummy:
        gra = without_dummy_atoms(gra)
    atm_sym_dct = atom_symbols(gra)
    nhvy_atms = sum(pt.to_Z(sym) != 1 for sym in atm_sym_dct.values())
    return nhvy_atms


def atom_neighbor_keys(gra):
    """ keys of neighboring atoms, by atom
    """
    def _neighbor_keys(atm_key, atm_nbh):
        return frozenset(atom_keys(atm_nbh) - {atm_key})

    atm_ngb_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(gra), _neighbor_keys)
    return atm_ngb_keys_dct


def atom_bond_keys(gra):
    """ bond keys, by atom
    """
    return dict_.transform_values(atom_neighborhoods(gra), bond_keys)


def atom_neighborhoods(gra):
    """ neighborhood subgraphs, by atom
    """
    bnd_keys = bond_keys(gra)

    def _neighborhood(atm_key):
        nbh_bnd_keys = set(filter(lambda x: atm_key in x, bnd_keys))
        return bond_induced_subgraph(gra, nbh_bnd_keys)

    atm_keys = list(atom_keys(gra))
    atm_nbh_dct = dict(zip(atm_keys, map(_neighborhood, atm_keys)))
    return atm_nbh_dct


# # bond properties
def bond_neighbor_keys(gra):
    """ keys of neighboring bonds, by bond
    """
    def _neighbor_keys(bnd_key, bnd_nbh):
        return frozenset(bond_keys(bnd_nbh) - {bnd_key})
    bnd_ngb_keys_dct = dict_.transform_items_to_values(
        bond_neighborhoods(gra), _neighbor_keys)
    return bnd_ngb_keys_dct


def bond_neighbor_bonds(bnd_key, gra):
    """ keys of neighboring bonds, by bond
    """
    atmi, atmj = list(bnd_key)
    ngb_atm_dct = atom_neighbor_keys(gra)
    bnds = []
    for atm in [atmi, atmj]:
        alpha_atms = ngb_atm_dct[atm]
        for alpha_atm in alpha_atms:
            if alpha_atm not in [atmi, atmj]:
                bnds.append(frozenset({atm, alpha_atm}))
    return bnds


def bond_neighborhoods(gra):
    """ neighborhood subgraphs, by bond
    """
    bnd_keys = list(bond_keys(gra))

    def _neighborhood(bnd_key):
        nbh_bnd_keys = set(filter(lambda x: bnd_key & x, bnd_keys))
        return bond_induced_subgraph(gra, nbh_bnd_keys)

    bnd_nbh_dct = dict(zip(bnd_keys, map(_neighborhood, bnd_keys)))
    return bnd_nbh_dct


# # other properties
def branch(gra, atm_key, bnd_key, saddle=False, ts_bnd=None):
    """ branch extending along `bnd_key` away from `atm_key`
    """
    return bond_induced_subgraph(
        gra,
        branch_bond_keys(gra, atm_key, bnd_key, saddle=saddle, ts_bnd=ts_bnd),
        saddle=saddle)


def branch_atom_keys(gra, atm_key, bnd_key, saddle=False, ts_bnd=None):
    """ atom keys for branch extending along `bnd_key` away from `atm_key`
    """
    bnch_atm_keys = atom_keys(
        branch(gra, atm_key, bnd_key, saddle=saddle, ts_bnd=ts_bnd))
    return bnch_atm_keys - {atm_key}


def branch_bond_keys(gra, atm_key, bnd_key, saddle=False, ts_bnd=None):
    """ bond keys for branch extending along `bnd_key` away from `atm_key`
    """

    # bnd_key is the set of atom indices for the bond of interest
    # atm_bnd_keys_dct is a dictionary of atoms that are connected to each atom
    bnd_key = frozenset(bnd_key)
    assert atm_key in bnd_key
    if not saddle:
        assert bnd_key in bond_keys(gra)

    atm_bnd_keys_dct = atom_bond_keys(gra)

    bnch_bnd_keys = {bnd_key}
    seen_bnd_keys = set()
    excl_bnd_keys = atm_bnd_keys_dct[atm_key] - {bnd_key}

    new_bnd_keys = {bnd_key}

    bnd_ngb_keys_dct = bond_neighbor_keys(gra)
    if ts_bnd:
        fts_bnd = frozenset(ts_bnd)
        bnd_ngb_keys_dct[fts_bnd] = bond_neighbor_bonds(ts_bnd, gra)
    while new_bnd_keys:
        new_bnd_ngb_keys = set(
            itertools.chain(
                *dict_.values_by_key(bnd_ngb_keys_dct, new_bnd_keys)))
        bnch_bnd_keys.update(new_bnd_ngb_keys - excl_bnd_keys)
        seen_bnd_keys.update(new_bnd_keys)
        new_bnd_keys = bnch_bnd_keys - seen_bnd_keys

    return frozenset(bnch_bnd_keys)


def rings(gra):
    """ rings in the graph (minimal basis)
    """
    gras = [bond_induced_subgraph(gra, bnd_keys)
            for bnd_keys in rings_bond_keys(gra)]
    return tuple(sorted(gras, key=frozen))


def rings_sorted_atom_keys(gra):
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
        map(_sorted_ring_atom_keys, rings_bond_keys(gra)))
    return rng_atm_keys_lst


def rings_bond_keys(gra):
    """ bond keys for each ring in the graph (minimal basis)
    """
    bnd_keys = bond_keys(gra)

    def _ring_bond_keys(rng_atm_keys):
        return frozenset(filter(lambda x: x <= rng_atm_keys, bnd_keys))

    nxg = _networkx.from_graph(gra)
    rng_atm_keys_lst = _networkx.minimum_cycle_basis(nxg)
    rng_bnd_keys_lst = frozenset(map(_ring_bond_keys, rng_atm_keys_lst))
    return rng_bnd_keys_lst


def connected_components(gra):
    """ connected components in the graph
    """
    cmp_gra_atm_keys_lst = connected_components_atom_keys(gra)
    cmp_gras = tuple(subgraph(gra, cmp_gra_atm_keys)
                     for cmp_gra_atm_keys in cmp_gra_atm_keys_lst)
    return cmp_gras


def connected_components_atom_keys(gra):
    """ atom keys for each connected component in the graph
    """
    nxg = _networkx.from_graph(gra)
    cmp_gra_atm_keys_lst = _networkx.connected_component_atom_keys(nxg)
    return cmp_gra_atm_keys_lst


def longest_chain(gra):
    """ longest chain in the graph
    """
    atm_keys = atom_keys(gra)

    max_chain = max((_longest_chain(gra, atm_key) for atm_key in atm_keys),
                    key=len)
    return max_chain


def atom_longest_chains(gra):
    """ longest chains, by atom
    """
    atm_keys = atom_keys(gra)

    long_chain_dct = {atm_key: _longest_chain(gra, atm_key)
                      for atm_key in atm_keys}
    return long_chain_dct


def _longest_chain(gra, atm_key):
    atm_ngb_keys_dct = atom_neighbor_keys(gra)
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    chains_lst = []
    if atm_ngb_keys:
        next_chains_lst = [
            [atm_key, atm_ngb_key] for atm_ngb_key in atm_ngb_keys]

        while True:
            chains_lst = next_chains_lst
            next_chains_lst = []
            for chain in chains_lst:
                atm_ngb_keys = atm_ngb_keys_dct[chain[-1]]
                next_atm_keys = sorted(atm_ngb_keys - set(chain))
                for next_atm_key in next_atm_keys:
                    next_chains_lst.append(chain + [next_atm_key])

            if not next_chains_lst:
                break

        max_chain = tuple(chains_lst[0])
    else:
        max_chain = tuple((atm_key,))
    return max_chain


def union(gra1, gra2):
    """ a union of two graphs
    """
    assert not atom_keys(gra1) & atom_keys(gra2)
    atm_dct = {}
    atm_dct.update(atoms(gra1))
    atm_dct.update(atoms(gra2))

    bnd_dct = {}
    bnd_dct.update(bonds(gra1))
    bnd_dct.update(bonds(gra2))
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


def union_from_sequence(gras):
    """ a union of all parts of a sequence of graphs
    """
    return tuple(functools.reduce(union, gras))


def subgraph(gra, atm_keys):
    """ the subgraph induced by a subset of the atoms
    """
    atm_keys = set(atm_keys)
    assert atm_keys <= atom_keys(gra)
    bnd_keys = set(filter(lambda x: x <= atm_keys, bond_keys(gra)))
    atm_dct = dict_.by_key(atoms(gra), atm_keys)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


def bond_induced_subgraph(gra, bnd_keys, saddle=False):
    """ the subgraph induced by a subset of the bonds
    """
    atm_keys = set(itertools.chain(*bnd_keys))
    bnd_keys = set(bnd_keys)
    assert atm_keys <= atom_keys(gra)
    if not saddle:
        assert bnd_keys <= bond_keys(gra)
    atm_dct = dict_.by_key(atoms(gra), atm_keys)
    bnd_dct = dict_.by_key(bonds(gra), bnd_keys)
    return _create.from_atoms_and_bonds(atm_dct, bnd_dct)


# # transformations
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


def add_atoms(gra, sym_dct, imp_hyd_vlc_dct=None, ste_par_dct=None):
    """ add atoms to this molecular graph, setting their keys
    """
    atm_keys = atom_keys(gra)
    atm_sym_dct = atom_symbols(gra)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(gra)
    atm_ste_par_dct = atom_stereo_parities(gra)

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
    bnd_dct = bonds(gra)
    gra = _create.from_atoms_and_bonds(atoms=atm_dct, bonds=bnd_dct)
    return gra


def add_bonded_atom(gra, sym, bnd_atm_key, imp_hyd_vlc=None, atm_ste_par=None,
                    bnd_ord=None, bnd_ste_par=None):
    """ add a single atom with a bond to an atom already in the graph
    """
    atm_keys = atom_keys(gra)

    atm_key = max(atm_keys) + 1

    sym_dct = {atm_key: sym}
    imp_hyd_vlc_dct = ({atm_key: imp_hyd_vlc}
                       if imp_hyd_vlc is not None else None)
    atm_ste_par_dct = ({atm_key: atm_ste_par}
                       if atm_ste_par is not None else None)

    gra = add_atoms(gra, sym_dct, imp_hyd_vlc_dct=imp_hyd_vlc_dct,
                    ste_par_dct=atm_ste_par_dct)

    bnd_key = frozenset({atm_key, bnd_atm_key})
    bnd_ord_dct = {bnd_key: bnd_ord} if bnd_ord is not None else None
    bnd_ste_par_dct = ({bnd_key: bnd_ste_par}
                       if bnd_ste_par is not None else None)

    gra = add_bonds(gra, [bnd_key], ord_dct=bnd_ord_dct,
                    ste_par_dct=bnd_ste_par_dct)

    return gra, atm_key


def remove_atoms(gra, atm_keys, check=True):
    """ remove atoms from the molecular graph
    """
    all_atm_keys = atom_keys(gra)
    atm_keys = set(atm_keys)

    if check:
        assert atm_keys <= all_atm_keys

    atm_keys_left = all_atm_keys - atm_keys
    return subgraph(gra, atm_keys_left)


def add_bonds(gra, keys, ord_dct=None, ste_par_dct=None, check=True):
    """ add bonds to this molecular graph
    """
    bnd_keys = set(bond_keys(gra))
    bnd_ord_dct = bond_orders(gra)
    bnd_ste_par_dct = bond_stereo_parities(gra)

    keys = set(map(frozenset, keys))
    ord_dct = {} if ord_dct is None else ord_dct
    ste_par_dct = {} if ste_par_dct is None else ste_par_dct

    if check:
        assert not keys & bnd_keys

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


def without_dummy_atoms(gra):
    """ remove dummy atoms from the graph
    """
    atm_sym_dct = atom_symbols(gra)
    atm_keys = [key for key, sym in atm_sym_dct.items() if pt.to_Z(sym)]
    return subgraph(gra, atm_keys)


# implicit/explicit hydrogen functions
# # atom properties
def atom_explicit_hydrogen_valences(gra):
    """ explicit hydrogen valences, by atom
    """
    return dict_.transform_values(atom_explicit_hydrogen_keys(gra), len)


def atom_explicit_hydrogen_keys(gra):
    """ explicit hydrogen valences, by atom
    """
    def _explicit_hydrogen_keys(atm_key, atm_nbh):
        return frozenset(explicit_hydrogen_keys(atm_nbh) - {atm_key})

    atm_exp_hyd_keys_dct = dict_.transform_items_to_values(
        atom_neighborhoods(gra), _explicit_hydrogen_keys)
    return atm_exp_hyd_keys_dct


# # other properties
def backbone_keys(gra):
    """ backbone atom keys
    """
    bbn_keys = atom_keys(gra) - explicit_hydrogen_keys(gra)
    return bbn_keys


def explicit_hydrogen_keys(gra):
    """ explicit hydrogen keys (H types: explicit, implicit, backbone)
    """
    hyd_keys = dict_.keys_by_value(atom_symbols(gra), lambda x: x == 'H')
    atm_ngb_keys_dct = atom_neighbor_keys(gra)

    def _is_backbone(hyd_key):
        return all(ngb_key in hyd_keys and hyd_key < ngb_key
                   for ngb_key in atm_ngb_keys_dct[hyd_key])

    exp_hyd_keys = frozenset(fmit.filterfalse(_is_backbone, hyd_keys))
    return exp_hyd_keys


# # transformations
def add_atom_explicit_hydrogen_keys(gra, atm_exp_hyd_keys_dct):
    """ add explicit hydrogens by atom
    """
    assert set(atm_exp_hyd_keys_dct.keys()) <= atom_keys(gra)
    for atm_key, atm_exp_hyd_keys in atm_exp_hyd_keys_dct.items():
        assert not set(atm_exp_hyd_keys) & atom_keys(gra)
        atm_exp_hyd_bnd_keys = {frozenset({atm_key, atm_exp_hyd_key})
                                for atm_exp_hyd_key in atm_exp_hyd_keys}
        atm_exp_hyd_sym_dct = dict_.by_key({}, atm_exp_hyd_keys, fill_val='H')
        gra = add_atoms(gra, atm_exp_hyd_sym_dct)
        gra = add_bonds(gra, atm_exp_hyd_bnd_keys)
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
    gra = remove_atoms(gra, exp_hyd_keys)
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


# # comparisons
def full_isomorphism(gra1, gra2):
    """ full graph isomorphism
    """
    assert gra1 == explicit(gra1) and gra2 == explicit(gra2)
    nxg1 = _networkx.from_graph(gra1)
    nxg2 = _networkx.from_graph(gra2)
    iso_dct = _networkx.isomorphism(nxg1, nxg2)
    return iso_dct


def full_subgraph_isomorphism(gra1, gra2):
    """ gra2 is fully isomorphic to a subgraph of gra1
    """
    assert gra1 == explicit(gra1) and gra2 == explicit(gra2)
    nxg1 = _networkx.from_graph(gra1)
    nxg2 = _networkx.from_graph(gra2)
    iso_dct = _networkx.subgraph_isomorphism(nxg1, nxg2)
    return iso_dct


def backbone_isomorphic(gra1, gra2):
    """ are these molecular graphs backbone isomorphic?
    """
    return backbone_isomorphism(gra1, gra2) is not None


def backbone_isomorphism(gra1, gra2):
    """ graph backbone isomorphism

    for implicit graphs, this is the relabeling of `gra1` to produce `gra2`
    for other graphs, it gives the correspondences between backbone atoms
    """
    gra1 = implicit(gra1)
    gra2 = implicit(gra2)
    nxg1 = _networkx.from_graph(gra1)
    nxg2 = _networkx.from_graph(gra2)
    iso_dct = _networkx.isomorphism(nxg1, nxg2)
    return iso_dct


def backbone_unique(gras):
    """ unique non-isomorphic graphs from a series
    """
    gras = _unique(gras, equiv=backbone_isomorphic)
    return gras


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


# # atom properties
def atom_element_valences(gra):
    """ element valences (# possible single bonds), by atom
    """
    atm_sym_dct = atom_symbols(gra)
    atm_group_idx_dct = dict_.transform_values(atm_sym_dct, pt.to_group)
    atm_elem_vlc_dct = dict_.transform_values(atm_group_idx_dct,
                                              VALENCE_DCT.__getitem__)
    return atm_elem_vlc_dct


def atom_lone_pair_counts(gra):
    """ lone pair counts, by atom
    """
    atm_sym_dct = atom_symbols(gra)
    atm_group_idx_dct = dict_.transform_values(atm_sym_dct, pt.to_group)
    atm_lpc_dct = dict_.transform_values(atm_group_idx_dct,
                                         LONE_PAIR_COUNTS_DCT.__getitem__)
    atm_lpc_dct = dict_.transform_values(atm_lpc_dct, int)
    return atm_lpc_dct


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


def unsaturated_atom_keys(gra):
    """ keys of unsaturated (radical or pi-bonded) atoms
    """
    atm_unsat_vlc_dct = atom_unsaturated_valences(gra, bond_order=False)
    unsat_atm_keys = frozenset(dict_.keys_by_value(atm_unsat_vlc_dct, bool))
    return unsat_atm_keys


# # other properties
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


# miscellaneous
# # bond properties
def bond_symmetry_numbers(gra, frm_bnd_key=None, brk_bnd_key=None):
    """ symmetry numbers, by bond

    the (approximate) symmetry number of the torsional potential for this bond,
    based on the hydrogen counts for each atom
    It is reduced to 1 if one of the H atoms in the torsional bond is a
    neighbor to the special bonding atom (the atom that is being transferred)
    """
    imp_gra = implicit(gra)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(imp_gra)

    bnd_keys = bond_keys(imp_gra)

    tfr_atm = None
    if frm_bnd_key and brk_bnd_key:
        for atm_f in list(frm_bnd_key):
            for atm_b in list(brk_bnd_key):
                if atm_f == atm_b:
                    tfr_atm = atm_f

        if tfr_atm:
            neighbor_dct = atom_neighbor_keys(gra)
            nei_tfr = neighbor_dct[tfr_atm]

            gra = gra[0]
            all_hyds = []
            for atm in gra:
                if gra[atm][0] == 'H':
                    all_hyds.append(atm)
        else:
            nei_tfr = {}

    bnd_sym_num_dct = {}
    bnd_sym_nums = []
    for bnd_key in bnd_keys:
        bnd_sym = 1
        vlc = max(map(atm_imp_hyd_vlc_dct.__getitem__, bnd_key))
        if vlc == 3:
            bnd_sym = 3
            if tfr_atm:
                for atm in nei_tfr:
                    nei_s = neighbor_dct[atm]
                    h_nei = 0
                    for nei in nei_s:
                        if nei in all_hyds:
                            h_nei += 1
                    if h_nei == 3:
                        bnd_sym = 1
        bnd_sym_nums.append(bnd_sym)

    bnd_sym_num_dct = dict(zip(bnd_keys, bnd_sym_nums))

    # fill in the rest of the bonds for completeness
    bnd_sym_num_dct = dict_.by_key(bnd_sym_num_dct, bond_keys(gra), fill_val=1)
    return bnd_sym_num_dct


if __name__ == '__main__':
    print(standard_keys(({5: ('H', 0, None)}, {})))
