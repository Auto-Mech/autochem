""" graph functions associated with kekule (resonance) structures

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
from typing import List

import numpy

from automol.graph.base._0core import (
    atom_bond_counts,
    atom_keys,
    atom_lone_pairs,
    atom_unpaired_electrons,
    atoms,
    atoms_bond_keys,
    atoms_neighbor_atom_keys,
    bond_keys,
    bond_orders,
    bond_stereo_keys,
    bond_unpaired_electrons,
    dummy_source_dict,
    has_atom_stereo,
    implicit,
    is_ts_graph,
    set_bond_orders,
    subgraph,
    tetrahedral_atom_keys,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reacting_atom_keys,
    ts_reagents_graph_without_stereo,
    ts_transferring_atoms,
    without_dummy_atoms,
    without_pi_bonds,
)
from automol.graph.base._2algo import (
    branches,
    connected_components,
    connected_components_atom_keys,
    rings_atom_keys,
)
from automol.util import dict_


# # core functions
def kekule(gra, max_stereo_overlap=True):
    """One low-spin kekule graph, ignoring current bond orders

    Low-spin kekule graphs have double and triple bonds assigned to
    minimize the number of unpaired electrons.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param max_stereo_overlap: optionally, request as many stereo bonds as
        possible to have a bond order of 2
    :type max_stereo_overlap: bool
    :returns: a kekule graph
    """
    ste_bnd_keys = bond_stereo_keys(gra)

    def _count_stereo_double_bonds(gra):
        bnd_ord_dct = bond_orders(gra)
        count = sum(bnd_ord_dct[k] == 2 for k in ste_bnd_keys)
        return count

    gras = kekules(gra)
    if not max_stereo_overlap:
        gra = next(iter(gras))
    else:
        gra = max(gras, key=_count_stereo_double_bonds)

    return gra


def kekules(gra):
    """All possible low-spin kekule graphs, ignoring current bond orders

    Low-spin kekule graphs have double and triple bonds assigned to
    minimize the number of unpaired electrons.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: all possible low-spin kekule graphs
    """
    orig_gra = gra
    bnd_ord_dcts = kekules_bond_orders(gra)
    gras = []
    for bnd_ord_dct in bnd_ord_dcts:
        gra = set_bond_orders(orig_gra, bnd_ord_dct)
        gras.append(gra)
    return tuple(gras)


def kekule_bond_orders(gra, max_stereo_overlap=True):
    """Bond orders for one low-spin kekule graph

    Low-spin kekule graphs have double and triple bonds assigned to
    minimize the number of unpaired electrons.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param max_stereo_overlap: optionally, request as many stereo bonds as
        possible to have a bond order of 2
    :type max_stereo_overlap: bool
    :returns: a bond order dictionary
    """
    return bond_orders(kekule(gra, max_stereo_overlap=max_stereo_overlap))


def kekules_bond_orders(gra):
    """Bond orders for all possible low-spin kekule graphs

    Low-spin kekule graphs have double and triple bonds assigned to
    minimize the number of unpaired electrons.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: bond orders for all possible low-spin kekule graphs
    :rtype: tuple[dict]
    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    orig_bnd_ord_dct = bond_orders(gra)
    gra = implicit(gra)
    gra = without_pi_bonds(gra)

    # identify all of the independent pi systems and assign kekules to each
    pi_keys_lst = pi_system_atom_keys(gra)
    pi_bord_dcts_lst = [
        pi_system_kekules_bond_orders_brute_force(gra, pi_keys)
        for pi_keys in pi_keys_lst
    ]

    bnd_ord_dcts = []
    # combine the kekules from each pi system together in all possible ways
    for bord_dcts in itertools.product(*pi_bord_dcts_lst):
        bord_dct = orig_bnd_ord_dct.copy()
        for dct in bord_dcts:
            bord_dct.update(dct)
        bnd_ord_dcts.append(bord_dct)

    if bnd_ord_dcts:
        bnd_ord_dcts = tuple(bnd_ord_dcts)
    else:
        bnd_ord_dcts = (orig_bnd_ord_dct,)

    return bnd_ord_dcts


def kekules_bond_orders_collated(gra):
    """Bond orders for all possible low-spin kekule graphs, collated into a
    single dictionary

    For TS graphs, collates possible bond orders from both reactants and
    products.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: bond orders for all possible low-spin kekule graphs
    :rtype: tuple[dict]
    """
    if is_ts_graph(gra):
        gras = [
            ts_reagents_graph_without_stereo(gra, prod=False),
            ts_reagents_graph_without_stereo(gra, prod=True),
        ]
    else:
        gras = [gra]

    bnd_keys = list(bond_keys(gra))
    bnd_ords_lst = list(
        dict_.values_by_key(d, bnd_keys, fill_val=0)
        for g in gras
        for d in kekules_bond_orders(g)
    )
    bnd_ords_dct = dict(zip(bnd_keys, zip(*bnd_ords_lst)))
    return bnd_ords_dct


def kekules_bond_orders_averaged(gra):
    """Bond orders for all possible low-spin kekule graphs, averaged

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: bond orders for all possible low-spin kekule graphs
    :rtype: tuple[dict]
    """
    bnd_ords_dct = kekules_bond_orders_collated(gra)
    avg_bnd_ord_dct = {k: sum(v) / len(v) for k, v in bnd_ords_dct.items()}
    return avg_bnd_ord_dct


# # derived properties
def linear_atom_keys(gra, dummy=True):
    """Atoms forming linear bonds, based on their hybridization

    For TS graphs, includes atoms that are linear for *either* the reactants
    *or* the products. This both simplifies the way Reaction objects can be
    handled and anticipates cases where the TS structure if close to either
    reactants or products.

    :param gra: the graph
    :param dummy: whether or not to consider atoms connected to dummy atoms as
        linear, if different from what would be predicted based on their
        hybridization
    :returns: the linear atom keys
    :rtype: tuple[int]
    """
    ts_ = is_ts_graph(gra)
    if ts_:
        gras = [
            ts_reagents_graph_without_stereo(gra, prod=False),
            ts_reagents_graph_without_stereo(gra, prod=True),
        ]
    else:
        gras = [gra]

    lin_atm_keys = set()
    for gra_ in gras:
        gra_ = ts_reagents_graph_without_stereo(gra_)
        nkeys_dct = atoms_neighbor_atom_keys(gra_)
        # To be linear, an atom must be both (a.) sp1 hybridized and (b.) have more than
        # 1 neighbor (exactly 2)
        atm_hyb_dct = atom_hybridizations_from_kekule(implicit(kekule(gra_)))
        sp1_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 1)
        lin_atm_keys |= set(k for k in sp1_atm_keys if len(nkeys_dct[k]) > 1)

    # If requested, include all keys associated with dummy atoms
    if dummy:
        dum_ngb_key_dct = dummy_source_dict(gra, dir_=False)
        lin_atm_keys |= set(dum_ngb_key_dct.values())

    if ts_:
        lin_atm_keys |= ts_linear_reacting_atom_keys(gra, ring=False)

    return frozenset(lin_atm_keys)


def linear_segments_atom_keys(gra, lin_keys=None):
    """atom keys for linear segments in the graph"""
    ngb_keys_dct = atoms_neighbor_atom_keys(without_dummy_atoms(gra))

    lin_keys = linear_atom_keys(gra, dummy=True) if lin_keys is None else lin_keys
    lin_segs = connected_components(subgraph(gra, lin_keys))

    lin_keys_lst = []
    for lin_seg in lin_segs:
        lin_seg_keys = atom_keys(lin_seg)
        if len(lin_seg_keys) == 1:
            (key,) = lin_seg_keys
            lin_keys_lst.append([key])
        else:
            end_key1, end_key2 = sorted(
                [
                    key
                    for key, ngb_keys in atoms_neighbor_atom_keys(lin_seg).items()
                    if len(ngb_keys) == 1
                ]
            )
            ngb_keys_dct = atoms_neighbor_atom_keys(lin_seg)

            key = None
            keys = [end_key1]
            while key != end_key2:
                (key,) = ngb_keys_dct[keys[-1]] - set(keys)
                keys.append(key)
            lin_keys_lst.append(keys)

    lin_keys_lst = tuple(map(tuple, lin_keys_lst))
    return lin_keys_lst


def unneeded_dummy_atom_keys(gra) -> frozenset:
    """Get the keys of dummy atoms which are not connected to linear atoms

    :param gra: A graph
    :type gra: automol graph data structure
    :returns: Keys for the unneeded dummy atoms
    :rtype: frozenset[int]
    """
    dum_lin_key_dct = dummy_source_dict(gra, dir_=False)
    lin_keys = linear_atom_keys(gra, dummy=False)
    return frozenset(dk for dk, lk in dum_lin_key_dct.items() if lk not in lin_keys)


def atom_hybridizations(gra):
    """resonance-dominant atom hybridizations, by atom"""
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_keys = list(atom_keys(gra))
    atm_hybs_by_res = [
        dict_.values_by_key(atom_hybridizations_from_kekule(g), atm_keys)
        for g in kekules(gra)
    ]
    atm_hybs = [min(hybs) for hybs in zip(*atm_hybs_by_res)]
    atm_hyb_dct = dict(zip(atm_keys, atm_hybs))
    return atm_hyb_dct


def atom_hybridizations_from_kekule(gra):
    """atom hybridizations, by atom"""
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_keys = list(atom_keys(gra))
    atm_unsat_dct = atom_unpaired_electrons(gra, bond_order=True)
    atm_bnd_vlc_dct = atom_bond_counts(gra, bond_order=False)  # note!!
    atm_unsats = numpy.array(dict_.values_by_key(atm_unsat_dct, atm_keys))
    atm_bnd_vlcs = numpy.array(dict_.values_by_key(atm_bnd_vlc_dct, atm_keys))
    atm_lpcs = numpy.array(dict_.values_by_key(atom_lone_pairs(gra), atm_keys))
    atm_hybs = atm_unsats + atm_bnd_vlcs + atm_lpcs - 1
    atm_hyb_dct = dict_.transform_values(dict(zip(atm_keys, atm_hybs)), int)
    return atm_hyb_dct


def radical_atom_keys(gra, sing_res=False, min_valence=1.0):
    """Radical atom keys for this molecular graph

    Radical atoms are based on the lowest-spin resonance structures for this
    graph. If the `sing_res` flag is set, a single low-spin resonance
    structure will be chosen when there are multiple such structures.

    If you wish to identify radical atom keys based on the bond orders already
    in `gra`, this can be done by using the `radical_atom_keys_from_kekule`
    function.

    :param gra: the molecular graph
    :param sing_res: only include radical keys for a single (arbitrary)
        resonance structure, or include all atoms that are radicals in any of
        the low-spin resonance structures?
    :type sing_res: bool
    :param min_valence: optionally, specify that only sites with at least a
        certain number of radical electrons be included
    :type min_valence: int
    :returns: the radical atom keys
    :rtype: frozenset[int]

    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_keys = list(atom_keys(gra))

    if sing_res:
        atm_rad_vlcs = dict_.values_by_key(
            atom_unpaired_electrons(kekule(gra)), atm_keys
        )
    else:
        atm_rad_vlcs_by_res = [
            dict_.values_by_key(atom_unpaired_electrons(dom_gra), atm_keys)
            for dom_gra in kekules(gra)
        ]
        atm_rad_vlcs = [max(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]

    atm_rad_keys = frozenset(
        atm_key
        for atm_key, atm_rad_vlc in zip(atm_keys, atm_rad_vlcs)
        if atm_rad_vlc >= min_valence
    )
    return atm_rad_keys


def radical_atom_keys_from_kekule(gra, min_valence=1.0):
    """Radical atom keys for a particular kekule graph

    Assumes the graph already has assigned bond orders

    :param gra: a resonance-structure molecular graph
    :param min_valence: optionally, specify that only sites with at least a
        certain number of radical electrons be included
    :type min_valence: int
    :returns: the radical atom keys
    :rtype: frozenset[int]

    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_keys = list(atom_keys(gra))

    atm_rad_vlcs = dict_.values_by_key(atom_unpaired_electrons(gra), atm_keys)

    atm_rad_keys = frozenset(
        atm_key
        for atm_key, atm_rad_vlc in zip(atm_keys, atm_rad_vlcs)
        if atm_rad_vlc >= min_valence
    )
    return atm_rad_keys


def nonresonant_radical_atom_keys(gra):
    """keys for radical atoms that are not in resonance"""
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_keys = list(atom_keys(gra))
    atm_rad_vlcs_by_res = [
        dict_.values_by_key(atom_unpaired_electrons(g), atm_keys) for g in kekules(gra)
    ]
    atm_rad_vlcs = [min(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]
    atm_rad_keys = frozenset(
        atm_key for atm_key, atm_rad_vlc in zip(atm_keys, atm_rad_vlcs) if atm_rad_vlc
    )
    return atm_rad_keys


def vinyl_radical_atom_bond_keys(gra):
    """Vinyl radical atom keys for this molecular graph

    :param gra: the molecular graph
    :returns: the vinyl radical atom keys
    :rtype: frozenset[int]
    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_rad_keys = nonresonant_radical_atom_keys(gra)
    bnd_ords_dct = kekules_bond_orders_collated(gra)
    atm_bnd_keys_dct = atoms_bond_keys(gra)
    vin_dct = {}
    for atm_key in atm_rad_keys:
        for bnd_key in atm_bnd_keys_dct[atm_key]:
            if 2 in bnd_ords_dct[bnd_key]:
                vin_dct[atm_key] = bnd_key
                break
    return vin_dct


def sigma_radical_atom_bond_keys(gra):
    """keys for sigma radical atoms

    :param gra: the molecular graph
    :returns: the sigma radical atom keys
    :rtype: frozenset[int]
    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_rad_keys = nonresonant_radical_atom_keys(gra)
    bnd_ords_dct = kekules_bond_orders_collated(gra)
    atm_bnd_keys_dct = atoms_bond_keys(gra)
    sig_dct = {}
    for atm_key in atm_rad_keys:
        for bnd_key in atm_bnd_keys_dct[atm_key]:
            if 3 in bnd_ords_dct[bnd_key]:
                sig_dct[atm_key] = bnd_key
                break
    return sig_dct


def vinyl_radical_atom_keys(gra):
    """Vinyl radical atom keys for this molecular graph

    :param gra: the molecular graph
    :returns: the vinyl radical atom keys
    :rtype: frozenset[int]
    """
    return frozenset(vinyl_radical_atom_bond_keys(gra))


def sigma_radical_atom_keys(gra):
    """keys for sigma radical atoms

    :param gra: the molecular graph
    :returns: the sigma radical atom keys
    :rtype: frozenset[int]
    """
    return frozenset(sigma_radical_atom_bond_keys(gra))


def ts_reacting_electron_direction(tsg, key: int):
    """Determine the reacting electron direction at one end of a forming bond

    Does *not* account for stereochemistry

    The direction is determined as follows:
        1. One bond, defining the 'x' axis direction
        2. Another bond, defining the 'y' axis direction
        3. An angle, describing how far to rotate the 'x' axis bond about a right-handed
        'z'-axis in order to arrive at the appropriate direction

    The 'y'-axis bond is `None` if the direction is parallel or antiparallel
    to the 'x'-axis bond, or if the orientation doesn't matter.

    Both bonds are `None` if the direction is perpendicular to the 'x-y' plane.

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param key: The key of a bond-forming atom
    :type key: int
    :returns: Two directed bond keys (x and y, respectively) and an angle
    :rtype: (Tuple[int, int], Tuple[int, int], float)
    """
    assert key in ts_reacting_atom_keys(tsg), f"Atom {key} is not a reacting atom:{tsg}"
    rcts_gra = ts_reagents_graph_without_stereo(tsg)
    tra_dct = ts_transferring_atoms(tsg)
    nkeys_dct = atoms_neighbor_atom_keys(rcts_gra)
    vin_dct = vinyl_radical_atom_bond_keys(rcts_gra)
    sig_dct = sigma_radical_atom_bond_keys(rcts_gra)

    if key in tra_dct:
        # key1 = transferring atom key
        # key2 = donor atom
        dkey, _ = tra_dct[key]
        xbnd_key = (key, dkey)
        ybnd_key = None
        phi = numpy.pi
    elif key in vin_dct:
        # key1 = this key
        # key2 = opposite end of the vinyl bond
        (opp_key,) = vin_dct[key] - {key}
        nkey = next(iter(nkeys_dct[key] - {key, opp_key}), None)
        xbnd_key = (key, opp_key)
        ybnd_key = None if nkey is None else (key, nkey)
        phi = 4.0 * numpy.pi / 3.0
    elif key in sig_dct:
        # key1 = attacking atom key
        # key2 = neighbor
        (nkey,) = sig_dct[key] - {key}
        xbnd_key = (key, nkey)
        ybnd_key = None
        phi = numpy.pi
    else:
        xbnd_key = None
        ybnd_key = None
        phi = None

    return xbnd_key, ybnd_key, phi


def ts_linear_reacting_atom_keys(
    tsg, breaking: bool = True, ring: bool = False
) -> List[int]:
    """Identify linear reacting atoms in a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param breaking: Include breaking bonds?; default `True`
    :type breaking: bool, optional
    :param ring: Include atoms in rings?; default `False`
    :type ring: bool, optional
    :returns: `True` if it is, `False` if it isn't
    :rtype: bool
    """
    keys = set(itertools.chain(*ts_forming_bond_keys(tsg)))
    if breaking:
        keys |= set(itertools.chain(*ts_breaking_bond_keys(tsg)))

    if not ring:
        keys -= set(itertools.chain(*rings_atom_keys(tsg)))

    lin_keys = set()
    for key in keys:
        _, _, phi = ts_reacting_electron_direction(tsg, key)
        if phi is not None and numpy.allclose(phi, numpy.pi):
            lin_keys.add(key)

    return frozenset(lin_keys)


def has_separated_radical_sites(gra):
    """does this radical have two or more separated radical sites?

    The identification is performed based on one of its lowest-spin resonance
    structures. It shouldn't matter which of the low-spin resonance structures
    is used -- if one of them has separated radical sites, they all should.

    This identifies polyradical molecules, but excludes things like carbenes.

    :param gra: the graph
    :returns: True if it has, False if not
    :rtype: bool
    """
    rad_atm_keys = radical_atom_keys(gra, sing_res=True)
    return len(rad_atm_keys) > 1


def resonance_bond_stereo_keys(gra):
    """does this graph have stereo at a resonance bond?

    :param gra: the molecular graph
    :rtype: bool
    """
    ste_bnd_keys = bond_stereo_keys(gra)
    res_bnd_ords_dct = kekules_bond_orders_collated(gra)

    res_bnd_keys = []
    for bnd_key in ste_bnd_keys:
        if bnd_key in res_bnd_ords_dct and 1 in res_bnd_ords_dct[bnd_key]:
            res_bnd_keys.append(bnd_key)

    res_bnd_keys = tuple(res_bnd_keys)
    return res_bnd_keys


def vinyl_bond_stereo_keys(gra):
    """does this graph have stereo at a resonance bond?

    :param gra: the molecular graph
    :rtype: bool
    """
    ste_bnd_keys = bond_stereo_keys(gra)
    vin_atm_keys = vinyl_radical_atom_keys(gra)

    vin_bnd_keys = []
    for bnd_key in ste_bnd_keys:
        if bnd_key & vin_atm_keys:
            vin_bnd_keys.append(bnd_key)
    vin_bnd_keys = tuple(vin_bnd_keys)
    return vin_bnd_keys


def has_resonance_bond_stereo(gra):
    """does this graph have stereo at a resonance bond?

    :param gra: the molecular graph
    :rtype: bool
    """
    return bool(resonance_bond_stereo_keys(gra))


def has_vinyl_bond_stereo(gra):
    """does this graph have stereo at a vinyl bond?

    :param gra: the molecular graph
    :rtype: bool
    """
    return bool(vinyl_bond_stereo_keys(gra))


def has_nonkekule_bond_stereo(gra):
    """does this graph have stereo at a resonance bond that cannot be
    represented as a double bond in a kekule structure?

    :param gra: the molecular graph
    :rtype: bool
    """
    ste_bnd_keys = bond_stereo_keys(gra)
    bnd_ord_dcts = kekules_bond_orders(gra)

    nsingles_lst = []
    for bnd_ord_dct in bnd_ord_dcts:
        nsingles = len([k for k in ste_bnd_keys if bnd_ord_dct[k] == 1])
        nsingles_lst.append(nsingles)

    return min(nsingles_lst) > 0


def has_noninchi_stereo(gra):
    """does this graph have stereo that cannot be captured by an InChI string?

    :param gra: the molecular graph
    :rtype: bool
    """
    return (
        has_nonkekule_bond_stereo(gra)
        or has_vinyl_bond_stereo(gra)
        or has_atom_stereo(gra, symb="N")
    )


def radical_groups(gra):
    """returns a list of lists of groups attached each radical"""
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"

    groups = []
    rads = radical_atom_keys(gra, sing_res=True)
    for rad in rads:
        groups.append(branches(gra, rad))
    return groups


def radical_group_dct(gra):
    """return a dictionary of lists of groups attached each radical"""
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"

    groups = {}
    rads = list(radical_atom_keys(gra, sing_res=True))
    atms = atoms(gra)
    for rad in rads:
        key = atms[rad][0]
        if key in groups:
            groups[atms[rad][0]] += branches(gra, rad)
        else:
            groups[atms[rad][0]] = branches(gra, rad)

    return groups


def rigid_planar_bond_keys(gra):
    """Get keys to stereo candidate bonds, which will be stereogenic if their
    groups are distinct on either side

    Bonds will be considered stereo candidates if they are rigid and planar,
    which happens when:
        a. They are doubly bonded in at least one low-spin resonance structure
        b. AND both atoms are sp2, excluding cumulenes which are linear rather
        than planar and have sp1 atoms

    For TS graphs, bond keys in which both atoms are stereo candidates in the
    TS are excluded, to prevent inconsistency in, for example, elimination
    reactions.

    Question: Is it possible for a bond to be rigid and planar for the TS and
    the products, but not the reactants, without being captured by atom
    stereochemistry? I don't think so, and this implementation assumes not.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: The bond keys
    :rtype: frozenset[frozenset[int]]
    """
    ts_ = is_ts_graph(gra)
    if ts_:
        gras = [
            ts_reagents_graph_without_stereo(gra, prod=False),
            ts_reagents_graph_without_stereo(gra, prod=True),
        ]
    else:
        gras = [gra]

    rp_bnd_keys = set()
    for gra_ in gras:
        gra_ = without_pi_bonds(gra_)
        double_bnd_keys = dict_.keys_by_value(
            kekules_bond_orders_collated(gra_), lambda x: 2 in x
        )

        # make sure both ends are sp^2 (excludes cumulenes)
        atm_hyb_dct = atom_hybridizations(gra_)
        sp2_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 2)
        rp_bnd_keys |= {k for k in double_bnd_keys if k <= sp2_atm_keys}

    if ts_:
        tet_atm_keys = tetrahedral_atom_keys(gra)
        rp_bnd_keys = frozenset({k for k in rp_bnd_keys if k != k & tet_atm_keys})

    return rp_bnd_keys


def atom_centered_cumulene_keys(gra):
    """resonance dominant keys for atom-centered cumulenes

    the bond-centered cumulenes are described by
        (frozenset({end_atm_key1, end_atm_key2}), cent_atm_key)
    where the first pair contains the sp2 atoms at the cumulene ends and
    `cent_atm_key` is the key of the central atom
    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    cum_chains = _cumulene_chains(gra)
    cum_keys = set()
    for cum_chain in cum_chains:
        size = len(cum_chain)
        if size % 2 == 1:
            cum_keys.add(
                (frozenset({cum_chain[0], cum_chain[-1]}), cum_chain[size // 2])
            )
    cum_keys = frozenset(cum_keys)
    return cum_keys


def bond_centered_cumulene_keys(gra):
    """resonance dominant keys for bond-centered cumulenes

    the bond-centered cumulenes are described by
        (frozenset({end_atm_key1, end_atm_key2}),
         frozenset({cent_atm_key1, cent_atm_key2}))
    where the first pair contains the sp2 atoms at the cumulene ends and the
    second pair is the bond key for the central bond
    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    cum_chains = _cumulene_chains(gra)
    cum_keys = set()
    for cum_chain in cum_chains:
        size = len(cum_chain)
        if size % 2 == 0:
            cum_keys.add(
                (
                    frozenset({cum_chain[0], cum_chain[-1]}),
                    frozenset({cum_chain[size // 2 - 1], cum_chain[size // 2]}),
                )
            )
    cum_keys = frozenset(cum_keys)
    return cum_keys


# # helpers
def pi_system_atom_keys(gra):
    """Extract keys for each closed, connected pi-system of a molecule

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: keys for each closed, connected pi system
    """
    atm_unsat_dct = atom_unpaired_electrons(gra, bond_order=False)
    all_pi_keys = dict_.keys_by_value(atm_unsat_dct)
    pi_keys_lst = tuple(
        ks
        for ks in connected_components_atom_keys(subgraph(gra, all_pi_keys))
        if len(ks) > 1
    )
    return pi_keys_lst


def pi_system_kekules_bond_orders_brute_force(gra, pi_keys, log=False):
    """Determine kekules for a closed, connected pi-system

    A general brute-force algorithm

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pi_keys: keys of an closed, connected pi system
    :type pi_keys: frozenset[int]
    :param log: print a debugging log of the number of recursive calls?
    :type log: bool
    """
    pi_sy = subgraph(gra, pi_keys)
    atm_keys = list(atom_keys(pi_sy))
    bnd_keys = list(bond_keys(pi_sy))

    aus_dct = dict_.by_key(atom_unpaired_electrons(gra, bond_order=False), atm_keys)
    bus_dct = dict_.by_key(bond_unpaired_electrons(gra, bond_order=False), bnd_keys)

    spin_max = spin_min = sum(aus_dct.values())

    # Allows no more than 2 copies of a given key in the pool, to prohibit
    # quadruple bonding for [C][C]
    bnd_pool = list(
        itertools.chain(*(itertools.repeat(k, min(bus_dct[k], 2)) for k in bnd_keys))
    )

    bnd_ord_dcts = []
    niter = 0
    for inc in range(spin_max // 2, 0, -1):
        spin = spin_max - inc * 2

        if spin > spin_min:
            break

        for bkeys in itertools.combinations(bnd_pool, inc):
            niter += 1

            akeys = list(itertools.chain(*bkeys))
            excess = next((k for k in set(akeys) if aus_dct[k] < akeys.count(k)), None)
            if excess is None:
                spin_min = spin
                bnd_ord_dct = {k: 1 + bkeys.count(k) for k in set(bkeys)}
                if bnd_ord_dct not in bnd_ord_dcts:
                    bnd_ord_dcts.append(bnd_ord_dct)

    if log:
        print(f"niter: {niter}")

    return bnd_ord_dcts


def _cumulene_chains(gra):
    atm_hyb_dct = atom_hybridizations(gra)
    sp1_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 1)
    sp2_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 2)

    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    def _cumulene_chain(chain):
        ret = None
        atm_key = chain[-1]
        next_atm_keys = atm_ngb_keys_dct[atm_key] - {chain[-2]}
        if next_atm_keys:
            assert len(next_atm_keys) == 1
            (next_atm_key,) = next_atm_keys
            if next_atm_key in sp1_atm_keys:
                chain.append(next_atm_key)
                ret = _cumulene_chain(chain)
            elif next_atm_key in sp2_atm_keys:
                chain.append(next_atm_key)
                ret = chain
        return ret

    cum_chains = []
    for atm_key in sp2_atm_keys:
        sp1_atm_ngb_keys = atm_ngb_keys_dct[atm_key] & sp1_atm_keys
        chains = [[atm_key, atm_ngb_key] for atm_ngb_key in sp1_atm_ngb_keys]
        for chain in chains:
            cum_chain = _cumulene_chain(chain)
            if cum_chain is not None:
                cum_chains.append(cum_chain)

    cum_chains = tuple(map(tuple, cum_chains))
    return cum_chains


if __name__ == "__main__":
    # # C=C[CH2]
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})

    # # C=CC=C[CH2]
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 2, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({1, 2}): (1, None)})

    # # C=C=C=C
    # GRA = ({0: ('C', 2, None), 1: ('C', 0, None), 2: ('C', 0, None),
    #         3: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None)})

    # # C=CC=CC=C
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 2, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None)})

    # # C=C-C(-[CH2])=C
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 2, None), 6: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({2, 6}): (1, None)})

    # # C1=CC=CC=C1 (benzene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({0, 5}): (1, None)})

    # # C12=CC=C1C=C2  _  _
    # #              ||_||_||
    # GRA = ({0: ('C', 0, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({0, 3}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({0, 5}): (1, None)})

    # # C1=CC=C2C=CC=CC2=C1 (naphthalene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({6, 7}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({8, 3}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({0, 9}): (1, None)})

    # # C1=CC=C2C=C3C=CC=CC3=CC2=C1 (anthracene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 0, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 1, None),
    #         9: ('C', 1, None), 10: ('C', 0, None), 11: ('C', 1, None),
    #         12: ('C', 0, None), 13: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({0, 13}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({11, 12}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
    #         frozenset({8, 7}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({3, 12}): (1, None), frozenset({10, 5}): (1, None),
    #         frozenset({12, 13}): (1, None), frozenset({10, 11}): (1, None)})

    # # C1=CC=C2C(=C1)C=CC3=CC=CC=C32 (phenanthrene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None), 10: ('C', 1, None), 11: ('C', 1, None),
    #         12: ('C', 1, None), 13: ('C', 0, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({4, 6}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({11, 12}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
    #         frozenset({3, 13}): (1, None), frozenset({0, 5}): (1, None),
    #         frozenset({8, 7}): (1, None), frozenset({8, 13}): (1, None),
    #         frozenset({12, 13}): (1, None), frozenset({10, 11}): (1, None)})

    # # C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2 (pyrene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None), 10: ('C', 1, None), 11: ('C', 1, None),
    #         12: ('C', 0, None), 13: ('C', 0, None), 14: ('C', 1, None),
    #         15: ('C', 1, None)},
    #        {frozenset({4, 6}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({11, 12}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 5}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({14, 15}): (1, None), frozenset({10, 11}): (1, None),
    #         frozenset({3, 4}): (1, None), frozenset({3, 13}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({6, 7}): (1, None),
    #         frozenset({8, 13}): (1, None), frozenset({12, 14}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({2, 15}): (1, None),
    #         frozenset({12, 13}): (1, None)})

    # # C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7 (coronene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 0, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 0, None), 10: ('C', 0, None), 11: ('C', 0, None),
    #         12: ('C', 1, None), 13: ('C', 1, None), 14: ('C', 1, None),
    #         15: ('C', 1, None), 16: ('C', 0, None), 17: ('C', 0, None),
    #         18: ('C', 0, None), 19: ('C', 0, None), 20: ('C', 1, None),
    #         21: ('C', 1, None), 22: ('C', 1, None), 23: ('C', 1, None)},
    #        {frozenset({17, 10}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({11, 12}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 5}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({9, 4}): (1, None), frozenset({14, 15}): (1, None),
    #         frozenset({10, 11}): (1, None), frozenset({3, 4}): (1, None),
    #         frozenset({16, 17}): (1, None), frozenset({22, 23}): (1, None),
    #         frozenset({16, 15}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({18, 19}): (1, None), frozenset({18, 3}): (1, None),
    #         frozenset({11, 14}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({2, 21}): (1, None),
    #         frozenset({17, 18}): (1, None), frozenset({8, 13}): (1, None),
    #         frozenset({20, 21}): (1, None), frozenset({19, 20}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({16, 23}): (1, None),
    #         frozenset({19, 22}): (1, None), frozenset({12, 13}): (1, None)})

    # # [C][C]
    # GRA = ({0: ('C', 0, None), 1: ('C', 0, None)},
    #        {frozenset({0, 1}): (1, None)})

    # [CH]=CCC#CC1C=CC=CC=1
    GRA = (
        {
            0: ("C", 1, None),
            1: ("C", 1, None),
            2: ("C", 2, None),
            3: ("C", 1, None),
            4: ("C", 0, None),
            5: ("C", 1, None),
            6: ("C", 1, None),
            7: ("C", 0, None),
            8: ("C", 1, None),
            9: ("C", 1, None),
            10: ("C", 0, None),
        },
        {
            frozenset({9, 6}): (1, None),
            frozenset({9, 10}): (1, None),
            frozenset({10, 7}): (1, None),
            frozenset({1, 2}): (1, None),
            frozenset({0, 1}): (1, True),
            frozenset({3, 6}): (1, None),
            frozenset({8, 10}): (1, None),
            frozenset({2, 4}): (1, None),
            frozenset({3, 5}): (1, None),
            frozenset({8, 5}): (1, None),
            frozenset({4, 7}): (1, None),
        },
    )

    # # C#CC=C
    # GRA = ({0: ('C', 1, None), 1: ('C', 0, None), 2: ('C', 1, None),
    #         3: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None)})

    # # C#CC1C=CC=CC=1
    # GRA = ({3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 0, None), 8: ('C', 1, None),
    #         9: ('C', 1, None), 10: ('C', 0, None)},
    #        {frozenset({9, 6}): (1, None), frozenset({9, 10}): (1, None),
    #         frozenset({10, 7}): (1, None), frozenset({3, 6}): (1, None),
    #         frozenset({8, 10}): (1, None), frozenset({3, 5}): (1, None),
    #         frozenset({8, 5}): (1, None), frozenset({4, 7}): (1, None)})

    print(len(kekules(GRA)))
