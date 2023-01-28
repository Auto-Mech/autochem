""" graph functions associated with kekule (resonance) structures

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
import numpy
from automol.util import dict_
from automol.graph.base._core import subgraph
from automol.graph.base._core import implicit
from automol.graph.base._core import atoms
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bond_keys
from automol.graph.base._core import bond_orders
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import set_bond_orders
from automol.graph.base._core import atom_unsaturations
from automol.graph.base._core import bond_unsaturations
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atoms_bond_keys
from automol.graph.base._core import atom_bond_valences
from automol.graph.base._core import atom_lone_pair_counts
from automol.graph.base._core import dummy_atoms_neighbor_atom_key
from automol.graph.base._core import without_bond_orders
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import from_ts_graph
from automol.graph.base._algo import atom_groups
from automol.graph.base._algo import connected_components
from automol.graph.base._algo import connected_components_atom_keys
# from automol.graph.base._algo import weighted_maximal_matching


# # core functions
def kekule(gra, max_stereo_overlap=True):
    """ One low-spin kekule graph, ignoring current bond orders

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
    """ All possible low-spin kekule graphs, ignoring current bond orders

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
    """ Bond orders for one low-spin kekule graph

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
    """ Bond orders for all possible low-spin kekule graphs

        Low-spin kekule graphs have double and triple bonds assigned to
        minimize the number of unpaired electrons.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: bond orders for all possible low-spin kekule graphs
        :rtype: tuple[dict]
    """
    orig_bnd_ord_dct = bond_orders(gra)
    gra = implicit(from_ts_graph(gra))
    gra = without_bond_orders(gra)

    # identify all of the independent pi systems and assign kekules to each
    pi_keys_lst = pi_system_atom_keys(gra)
    pi_bord_dcts_lst = [
        pi_system_kekules_bond_orders_brute_force(gra, pi_keys)
        for pi_keys in pi_keys_lst]

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
    """ Bond orders for all possible low-spin kekule graphs, collated into a
        single dictionary

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: bond orders for all possible low-spin kekule graphs
        :rtype: tuple[dict]
    """
    bnd_keys = list(bond_keys(gra))
    bnd_ords_lst = list(dict_.values_by_key(d, bnd_keys) for d in
                        kekules_bond_orders(gra))
    bnd_ords_dct = dict(zip(bnd_keys, zip(*bnd_ords_lst)))
    return bnd_ords_dct


def kekules_bond_orders_averaged(gra):
    """ Bond orders for all possible low-spin kekule graphs, averaged

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: bond orders for all possible low-spin kekule graphs
        :rtype: tuple[dict]
    """
    bnd_ords_dct = kekules_bond_orders_collated(gra)
    avg_bnd_ord_dct = {k: sum(v)/len(v) for k, v in bnd_ords_dct.items()}
    return avg_bnd_ord_dct


# # derived properties
def linear_atom_keys(gra, dummy=True):
    """ atoms forming linear bonds, based on their hybridization

    :param gra: the graph
    :param dummy: whether or not to consider atoms connected to dummy atoms as
        linear, if different from what would be predicted based on their
        hybridization
    :returns: the linear atom keys
    :rtype: tuple[int]
    """
    gra = from_ts_graph(gra)
    atm_hyb_dct = atom_hybridizations_from_kekule(implicit(kekule(gra)))
    lin_atm_keys = set(dict_.keys_by_value(atm_hyb_dct, lambda x: x == 1))

    if dummy:
        dum_ngb_key_dct = dummy_atoms_neighbor_atom_key(gra)
        lin_atm_keys |= set(dum_ngb_key_dct.values())

    lin_atm_keys = tuple(sorted(lin_atm_keys))
    return lin_atm_keys


def linear_segments_atom_keys(gra, lin_keys=None):
    """ atom keys for linear segments in the graph
    """
    ngb_keys_dct = atoms_neighbor_atom_keys(without_dummy_atoms(gra))

    lin_keys = (linear_atom_keys(gra, dummy=True)
                if lin_keys is None else lin_keys)

    lin_keys = [k for k in lin_keys if len(ngb_keys_dct[k]) <= 2]

    lin_segs = connected_components(subgraph(gra, lin_keys))

    lin_keys_lst = []
    for lin_seg in lin_segs:
        lin_seg_keys = atom_keys(lin_seg)
        if len(lin_seg_keys) == 1:
            key, = lin_seg_keys
            lin_keys_lst.append([key])
        else:
            end_key1, end_key2 = sorted([
                key for key, ngb_keys in
                atoms_neighbor_atom_keys(lin_seg).items()
                if len(ngb_keys) == 1])
            ngb_keys_dct = atoms_neighbor_atom_keys(lin_seg)

            key = None
            keys = [end_key1]
            while key != end_key2:
                key, = ngb_keys_dct[keys[-1]] - set(keys)
                keys.append(key)
            lin_keys_lst.append(keys)

    lin_keys_lst = tuple(map(tuple, lin_keys_lst))
    return lin_keys_lst


def atom_hybridizations(gra):
    """ resonance-dominant atom hybridizations, by atom
    """
    gra = from_ts_graph(gra)
    atm_keys = list(atom_keys(gra))
    atm_hybs_by_res = [
        dict_.values_by_key(atom_hybridizations_from_kekule(g), atm_keys)
        for g in kekules(gra)]
    atm_hybs = [min(hybs) for hybs in zip(*atm_hybs_by_res)]
    atm_hyb_dct = dict(zip(atm_keys, atm_hybs))
    return atm_hyb_dct


def atom_hybridizations_from_kekule(gra):
    """ atom hybridizations, by atom
    """
    gra = from_ts_graph(gra)
    atm_keys = list(atom_keys(gra))
    atm_unsat_dct = atom_unsaturations(gra, bond_order=True)
    atm_bnd_vlc_dct = atom_bond_valences(gra, bond_order=False)     # note!!
    atm_unsats = numpy.array(
        dict_.values_by_key(atm_unsat_dct, atm_keys))
    atm_bnd_vlcs = numpy.array(dict_.values_by_key(atm_bnd_vlc_dct, atm_keys))
    atm_lpcs = numpy.array(
        dict_.values_by_key(atom_lone_pair_counts(gra), atm_keys))
    atm_hybs = atm_unsats + atm_bnd_vlcs + atm_lpcs - 1
    atm_hyb_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_hybs)), int)
    return atm_hyb_dct


def radical_atom_keys(gra, sing_res=False, min_valence=1.):
    """ Radical atom keys for this molecular graph

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
    gra = from_ts_graph(gra)
    atm_keys = list(atom_keys(gra))

    if sing_res:
        atm_rad_vlcs = dict_.values_by_key(
            atom_unsaturations(kekule(gra)), atm_keys)
    else:
        atm_rad_vlcs_by_res = [
            dict_.values_by_key(atom_unsaturations(dom_gra), atm_keys)
            for dom_gra in kekules(gra)]
        atm_rad_vlcs = [
            max(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]

    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs)
                             if atm_rad_vlc >= min_valence)
    return atm_rad_keys


def radical_atom_keys_from_kekule(gra, min_valence=1.):
    """ Radical atom keys for a particular kekule graph

    Assumes the graph already has assigned bond orders

    :param gra: a resonance-structure molecular graph
    :param min_valence: optionally, specify that only sites with at least a
        certain number of radical electrons be included
    :type min_valence: int
    :returns: the radical atom keys
    :rtype: frozenset[int]

    """
    gra = from_ts_graph(gra)
    atm_keys = list(atom_keys(gra))

    atm_rad_vlcs = dict_.values_by_key(
        atom_unsaturations(gra), atm_keys)

    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs)
                             if atm_rad_vlc >= min_valence)
    return atm_rad_keys


def nonresonant_radical_atom_keys(gra):
    """ keys for radical atoms that are not in resonance
    """
    gra = from_ts_graph(gra)
    atm_keys = list(atom_keys(gra))
    atm_rad_vlcs_by_res = [
        dict_.values_by_key(atom_unsaturations(g), atm_keys)
        for g in kekules(gra)]
    atm_rad_vlcs = [min(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]
    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs) if atm_rad_vlc)
    return atm_rad_keys


def vinyl_radical_atom_keys(gra):
    """ Vinyl radical atom keys for this molecular graph

    :param gra: the molecular graph
    :returns: the vinyl radical atom keys
    :rtype: frozenset[int]
    """
    gra = from_ts_graph(gra)
    atm_rad_keys = nonresonant_radical_atom_keys(gra)
    bnd_ords_dct = kekules_bond_orders_collated(gra)
    atm_bnd_keys_dct = atoms_bond_keys(gra)
    atm_vin_keys = []
    for atm_key in atm_rad_keys:
        for bnd_key in atm_bnd_keys_dct[atm_key]:
            if 2 in bnd_ords_dct[bnd_key]:
                atm_vin_keys.append(atm_key)
                break
    atm_vin_keys = frozenset(atm_vin_keys)
    return atm_vin_keys


def sigma_radical_atom_keys(gra):
    """ keys for sigma radical atoms

    :param gra: the molecular graph
    :returns: the sigma radical atom keys
    :rtype: frozenset[int]
    """
    gra = from_ts_graph(gra)
    atm_rad_keys = nonresonant_radical_atom_keys(gra)
    bnd_ords_dct = kekules_bond_orders_collated(gra)
    atm_bnd_keys_dct = atoms_bond_keys(gra)
    atm_sig_keys = []
    for atm_key in atm_rad_keys:
        for bnd_key in atm_bnd_keys_dct[atm_key]:
            if 3 in bnd_ords_dct[bnd_key]:
                atm_sig_keys.append(atm_key)
                break
    atm_sig_keys = frozenset(atm_sig_keys)
    return atm_sig_keys


def has_separated_radical_sites(gra):
    """ does this radical have two or more separated radical sites?

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
    """ does this graph have stereo at a resonance bond?

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
    """ does this graph have stereo at a resonance bond?

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
    """ does this graph have stereo at a resonance bond?

        :param gra: the molecular graph
        :rtype: bool
    """
    return bool(resonance_bond_stereo_keys(gra))


def has_vinyl_bond_stereo(gra):
    """ does this graph have stereo at a vinyl bond?

        :param gra: the molecular graph
        :rtype: bool
    """
    return bool(vinyl_bond_stereo_keys(gra))


def has_nonkekule_bond_stereo(gra):
    """ does this graph have stereo at a resonance bond that cannot be
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


def radical_groups(gra):
    """ returns a list of lists of groups attached each radical
    """
    gra = from_ts_graph(gra)

    groups = []
    rads = radical_atom_keys(gra, sing_res=True)
    for rad in rads:
        groups.append(atom_groups(gra, rad))
    return groups


def radical_group_dct(gra):
    """ return a dictionary of lists of groups attached each radical
    """
    gra = from_ts_graph(gra)

    groups = {}
    rads = list(radical_atom_keys(gra, sing_res=True))
    atms = atoms(gra)
    for rad in rads:
        key = atms[rad][0]
        if key in groups:
            groups[atms[rad][0]] += atom_groups(gra, rad)
        else:
            groups[atms[rad][0]] = atom_groups(gra, rad)

    return groups


def rigid_planar_bond_keys(gra):
    """ determine the sp2 bonds in this graph
    """
    gra = without_bond_orders(gra)
    bnd_keys = dict_.keys_by_value(
        kekules_bond_orders_collated(gra), lambda x: 2 in x)

    # make sure both ends are sp^2 (excludes cumulenes)
    atm_hyb_dct = atom_hybridizations(gra)
    sp2_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 2)
    bnd_keys = frozenset({bnd_key for bnd_key in bnd_keys
                          if bnd_key <= sp2_atm_keys})
    return bnd_keys


def atom_centered_cumulene_keys(gra):
    """ resonance dominant keys for atom-centered cumulenes

    the bond-centered cumulenes are described by
        (frozenset({end_atm_key1, end_atm_key2}), cent_atm_key)
    where the first pair contains the sp2 atoms at the cumulene ends and
    `cent_atm_key` is the key of the central atom
    """
    gra = from_ts_graph(gra)
    cum_chains = _cumulene_chains(gra)
    cum_keys = set()
    for cum_chain in cum_chains:
        size = len(cum_chain)
        if size % 2 == 1:
            cum_keys.add(
                (frozenset({cum_chain[0], cum_chain[-1]}),
                 cum_chain[size // 2])
            )
    cum_keys = frozenset(cum_keys)
    return cum_keys


def bond_centered_cumulene_keys(gra):
    """ resonance dominant keys for bond-centered cumulenes

    the bond-centered cumulenes are described by
        (frozenset({end_atm_key1, end_atm_key2}),
         frozenset({cent_atm_key1, cent_atm_key2}))
    where the first pair contains the sp2 atoms at the cumulene ends and the
    second pair is the bond key for the central bond
    """
    gra = from_ts_graph(gra)
    cum_chains = _cumulene_chains(gra)
    cum_keys = set()
    for cum_chain in cum_chains:
        size = len(cum_chain)
        if size % 2 == 0:
            cum_keys.add(
                (frozenset({cum_chain[0], cum_chain[-1]}),
                 frozenset({cum_chain[size // 2 - 1], cum_chain[size // 2]}))
            )
    cum_keys = frozenset(cum_keys)
    return cum_keys


# # helpers
def pi_system_atom_keys(gra):
    """ Extract keys for each closed, connected pi-system of a molecule

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: keys for each closed, connected pi system
    """
    atm_unsat_dct = atom_unsaturations(gra, bond_order=False)
    all_pi_keys = dict_.keys_by_value(atm_unsat_dct)
    pi_keys_lst = tuple(
        ks for ks in
        connected_components_atom_keys(subgraph(gra, all_pi_keys))
        if len(ks) > 1)
    return pi_keys_lst


def pi_system_kekules_bond_orders_brute_force(gra, pi_keys, log=False):
    """ Determine kekules for a closed, connected pi-system

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

    aus_dct = dict_.by_key(
        atom_unsaturations(gra, bond_order=False), atm_keys)
    bus_dct = dict_.by_key(
        bond_unsaturations(gra, bond_order=False), bnd_keys)

    spin_max = spin_min = sum(aus_dct.values())

    bnd_pool = list(itertools.chain(*(
        itertools.repeat(k, min(bus_dct[k], 2)) for k in bnd_keys)))

    inc_bkeys_iter = ((inc, bkeys) for inc in range(spin_max // 2, 0, -1)
                      for bkeys in itertools.combinations(bnd_pool, inc))

    bnd_ord_dcts = []
    niter = 0
    for inc, bkeys in inc_bkeys_iter:
        niter += 1
        spin = spin_max - inc * 2

        if spin > spin_min:
            break

        akeys = list(itertools.chain(*bkeys))
        excess = next(
            (k for k in set(akeys) if aus_dct[k] < akeys.count(k)), None)
        if excess is None:
            spin_min = spin
            bnd_ord_dct = {k: 1 + bkeys.count(k) for k in set(bkeys)}
            if bnd_ord_dct not in bnd_ord_dcts:
                bnd_ord_dcts.append(bnd_ord_dct)
        else:
            continue

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
            next_atm_key, = next_atm_keys
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


if __name__ == '__main__':
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
    GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 2, None),
            3: ('C', 1, None), 4: ('C', 0, None), 5: ('C', 1, None),
            6: ('C', 1, None), 7: ('C', 0, None), 8: ('C', 1, None),
            9: ('C', 1, None), 10: ('C', 0, None)},
           {frozenset({9, 6}): (1, None), frozenset({9, 10}): (1, None),
            frozenset({10, 7}): (1, None), frozenset({1, 2}): (1, None),
            frozenset({0, 1}): (1, True), frozenset({3, 6}): (1, None),
            frozenset({8, 10}): (1, None), frozenset({2, 4}): (1, None),
            frozenset({3, 5}): (1, None), frozenset({8, 5}): (1, None),
            frozenset({4, 7}): (1, None)})

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
