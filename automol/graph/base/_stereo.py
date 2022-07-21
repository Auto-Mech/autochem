""" graph functions that depend on stereo assignments

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools
import numpy
from automol import util
import automol.geom.base    # !!!!
from automol.graph.base._core import atom_keys
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import set_atom_stereo_parities
from automol.graph.base._core import set_bond_stereo_parities
from automol.graph.base._core import set_stereo_parities
from automol.graph.base._core import has_stereo
from automol.graph.base._core import from_ts_graph
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._algo import rings_atom_keys
from automol.graph.base._algo import branch
from automol.graph.base._algo import is_connected
from automol.graph.base._canon import stereogenic_atom_keys
from automol.graph.base._canon import stereogenic_bond_keys
from automol.graph.base._canon import stereogenic_keys
from automol.graph.base._canon import to_local_stereo
from automol.graph.base._canon import reflect_local_stereo
from automol.graph.base._canon import refine_priorities
from automol.graph.base._canon import parity_evaluator_from_geometry_
from automol.graph.base._canon import canonical_assignment_representation
from automol.graph.base._amchi import connected_amchi_with_indices


# # core functions
def expand_stereo(gra, enant=True, symeq=False):
    """ Obtain all possible stereoisomers of a graph, ignoring its assignments

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param enant: Include all enantiomers, or only canonical ones?
        :type enant: bool
        :param symeq: Include symmetrically equivalent stereoisomers?
        :type symeq: bool
        :returns: a series of molecular graphs for the stereoisomers
    """
    gpcs = expand_stereo_with_priorities_and_amchis(gra)
    gras = []
    seen_chis = []
    for gra_, _, chi in gpcs:
        is_inv = automol.amchi.base.is_inverted_enantiomer(chi)
        if enant or is_inv is not True:
            if symeq or chi not in seen_chis:
                gras.append(gra_)
                seen_chis.append(chi)

    return tuple(gras)


def expand_stereo_with_priorities_and_amchis(gra):
    """ Obtain all possible stereoisomers, along with their AMChI strings and
        canonical priorities

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: a sequence of graphs paired with AMChI strings
    """
    assert is_connected(gra), (
        "Currently requires a connected graph, but could be generalized.")

    bools = (False, True)
    gra = from_ts_graph(gra)
    gra = without_stereo_parities(gra)

    gps0 = None
    gps = [(gra, None)]

    # 1. Expand all possible stereoisomers, along with their priority mappings
    while gps0 != gps:
        gps0 = gps
        gps = []

        for gra_, pri_dct in gps0:
            pri_dct = refine_priorities(gra_, pri_dct=pri_dct)

            keys = stereogenic_keys(gra_, pri_dct=pri_dct)
            gras = [set_stereo_parities(gra_, dict(zip(keys, pars)))
                    for pars in itertools.product(bools, repeat=len(keys))]

            gps.extend([(g, pri_dct) for g in gras])

    # 2. Group enantiomers into pairs, sorting them so that the canonical
    # enantiomer comes first.
    loc_gras = []
    diast_gps = []
    enant_gps = []
    for ugra, upri_dct in gps:
        uloc_gra = to_local_stereo(ugra, pri_dct=upri_dct)
        rloc_gra = reflect_local_stereo(uloc_gra)

        ridx = next((i for i, g in enumerate(loc_gras) if g == rloc_gra), None)
        if ridx is None:
            diast_gps.append((ugra, upri_dct))
            loc_gras.append(uloc_gra)
        else:
            loc_gras.pop(ridx)
            rgra, rpri_dct = diast_gps.pop(ridx)

            urep = canonical_assignment_representation(ugra, upri_dct)
            rrep = canonical_assignment_representation(rgra, rpri_dct)

            if urep == rrep:
                diast_gps.extend([(ugra, upri_dct), (rgra, rpri_dct)])
                loc_gras.extend([uloc_gra, rloc_gra])
            else:
                if urep > rrep:
                    ugra, rgra = rgra, ugra
                    upri_dct, rpri_dct = rpri_dct, upri_dct
                enant_gps.append([(ugra, upri_dct), (rgra, rpri_dct)])

    gpcs = []
    for gra_, pri_dct in diast_gps:
        chi, _ = connected_amchi_with_indices(gra_, pri_dct=pri_dct,
                                              is_refl=None)
        gpcs.append((gra_, pri_dct, chi))

    for (ugra, upri_dct), (rgra, rpri_dct) in enant_gps:
        uchi, _ = connected_amchi_with_indices(ugra, pri_dct=upri_dct,
                                               is_refl=False)
        rchi, _ = connected_amchi_with_indices(ugra, pri_dct=upri_dct,
                                               is_refl=True)

        gpcs.append((ugra, upri_dct, uchi))
        gpcs.append((rgra, rpri_dct, rchi))

    # Sort them by AMChI string
    gpcs = tuple(sorted(gpcs, key=lambda x: x[-1]))
    return gpcs


# # stereo evaluation
def local_atom_stereo_parity_from_geometry(gra, atm_key, geo,
                                           geo_idx_dct=None):
    """ Determine the local stereo parity of an atom from a geometry

        Local stereo parities are given relative to the indices of neighboring
        atoms.

        :param gra: molecular graph with stereo parities
        :type gra: automol graph data structure
        :param atm_key: the atom whose parity is to be determined
        :type atm_key: int
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
    """
    atm_keys = atom_keys(gra)
    pri_dct = dict(zip(atm_keys, atm_keys))
    par_eval_ = parity_evaluator_from_geometry_(
        gra, geo, geo_idx_dct=geo_idx_dct)
    par = par_eval_(pri_dct)(atm_key)
    return par


def local_bond_stereo_parity_from_geometry(gra, bnd_key, geo,
                                           geo_idx_dct=None):
    """ Determine the local stereo parity of a bond from a geometry

        Local stereo parities are given relative to the indices of neighboring
        atoms.

        :param gra: molecular graph with stereo parities
        :type gra: automol graph data structure
        :param bnd_key: the bond whose parity is to be determined
        :type bnd_key: int
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
    """
    atm_keys = atom_keys(gra)
    pri_dct = dict(zip(atm_keys, atm_keys))
    par_eval_ = parity_evaluator_from_geometry_(
        gra, geo, geo_idx_dct=geo_idx_dct)
    par = par_eval_(pri_dct)(bnd_key)
    return par


# # stereo correction
def stereo_corrected_geometry(gra, geo, geo_idx_dct=None, local_stereo=False):
    """ Obtain a geometry corrected for stereo parities based on a graph

        :param gra: molecular graph with stereo parities
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
        :param local_stereo: is this graph using local instead of canonical
            stereo?
        :type local_stereo: bool
        :returns: a molecular geometry with corrected stereo
    """
    sgr = gra if local_stereo else to_local_stereo(gra)
    gra = without_stereo_parities(gra)

    if has_stereo(sgr):
        full_atm_par_dct = atom_stereo_parities(sgr)
        full_bnd_par_dct = bond_stereo_parities(sgr)

        atm_keys = set()
        bnd_keys = set()

        last_gra = None

        while last_gra != gra:
            last_gra = gra

            atm_keys.update(stereogenic_atom_keys(gra))
            bnd_keys.update(stereogenic_bond_keys(gra))

            atm_par_dct = {k: full_atm_par_dct[k] for k in atm_keys}
            bnd_par_dct = {k: full_bnd_par_dct[k] for k in bnd_keys}
            geo, gra = _local_atom_stereo_corrected_geometry(
                gra, atm_par_dct, geo, geo_idx_dct)
            geo, gra = _local_bond_stereo_corrected_geometry(
                gra, bnd_par_dct, geo, geo_idx_dct)

    return geo


def _local_atom_stereo_corrected_geometry(gra, atm_par_dct, geo,
                                          geo_idx_dct=None):
    """ Correct a geometry to match local atom stereo assignments.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param atm_par_dct: local atom parities (local means relative to the
            neighboring atom keys)
        :type atm_par_dct: dict
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
    """
    atm_keys = atom_keys(gra)
    ring_atm_keys = set(itertools.chain(*rings_atom_keys(gra)))
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(atm_keys)})

    # Create a parity evaluator
    pri_dct = dict(zip(atm_keys, atm_keys))
    par_eval_ = parity_evaluator_from_geometry_(
        gra, geo, geo_idx_dct=geo_idx_dct)
    par_ = par_eval_(pri_dct)

    ste_atm_keys = list(atm_par_dct.keys())
    for atm_key in ste_atm_keys:
        par = atm_par_dct[atm_key]
        curr_par = par_(atm_key)

        if curr_par != par:
            atm_ngb_keys = atm_ngb_keys_dct[atm_key]
            # for now, we simply exclude rings from the pivot keys
            # (will not work for stereo atom at the intersection of two rings)
            atm_piv_keys = list(atm_ngb_keys - ring_atm_keys)[:2]
            if len(atm_piv_keys) == 2:
                atm3_key, atm4_key = atm_piv_keys

                # get coordinates
                xyzs = automol.geom.base.coordinates(geo)
                atm_xyz = xyzs[geo_idx_dct[atm_key]]
                atm3_xyz = xyzs[geo_idx_dct[atm3_key]]
                atm4_xyz = xyzs[geo_idx_dct[atm4_key]]

                # do the rotation
                rot_axis = util.vec.unit_bisector(
                    atm3_xyz, atm4_xyz, orig_xyz=atm_xyz)

                rot_atm_keys = (
                    atom_keys(branch(gra, atm_key, {atm_key, atm3_key})) |
                    atom_keys(branch(gra, atm_key, {atm_key, atm4_key})))
                rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

                geo = automol.geom.rotate(
                    geo, rot_axis, numpy.pi, orig_xyz=atm_xyz, idxs=rot_idxs)
            else:
                # Handles nitrogens in rings
                assert len(atm_piv_keys) == 1
                assert len(atm_ngb_keys) == 3
                atm3_key, = atm_piv_keys
                atm1_key, atm2_key = atm_ngb_keys - {atm3_key}

                # get coordinates
                xyzs = automol.geom.base.coordinates(geo)
                atm_xyz = xyzs[geo_idx_dct[atm_key]]
                atm1_xyz = xyzs[geo_idx_dct[atm1_key]]
                atm2_xyz = xyzs[geo_idx_dct[atm2_key]]

                # do the rotation
                rot_axis = util.vec.unit_bisector(
                    atm1_xyz, atm2_xyz, orig_xyz=atm_xyz)

                rot_atm_keys = atom_keys(
                    branch(gra, atm_key, {atm_key, atm3_key}))
                rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

                geo = automol.geom.rotate(
                    geo, rot_axis, numpy.pi, orig_xyz=atm_xyz, idxs=rot_idxs)

        gra = set_atom_stereo_parities(gra, {atm_key: par})

    return geo, gra


def _local_bond_stereo_corrected_geometry(gra, bnd_par_dct, geo,
                                          geo_idx_dct=None):
    """ Correct a geometry to match local bond stereo assignments.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param bnd_par_dct: local bond parities (local means relative to the
            neighboring atom keys)
        :type bnd_par_dct: dict
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
    """
    atm_keys = atom_keys(gra)
    bnd_keys = list(bnd_par_dct.keys())

    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(atm_keys)})

    # Create a local parity evaluator
    pri_dct = dict(zip(atm_keys, atm_keys))
    par_eval_ = parity_evaluator_from_geometry_(
        gra, geo, geo_idx_dct=geo_idx_dct)
    par_ = par_eval_(pri_dct)

    for bnd_key in bnd_keys:
        par = bnd_par_dct[bnd_key]
        curr_par = par_(bnd_key)

        if curr_par != par:
            xyzs = automol.geom.base.coordinates(geo)

            atm1_key, atm2_key = bnd_key
            atm1_xyz = xyzs[geo_idx_dct[atm1_key]]
            atm2_xyz = xyzs[geo_idx_dct[atm2_key]]

            rot_axis = numpy.subtract(atm2_xyz, atm1_xyz)

            rot_atm_keys = atom_keys(
                branch(gra, atm1_key, {atm1_key, atm2_key}))

            rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

            geo = automol.geom.rotate(
                geo, rot_axis, numpy.pi, orig_xyz=atm1_xyz, idxs=rot_idxs)

        gra = set_bond_stereo_parities(gra, {bnd_key: par})

    return geo, gra


if __name__ == '__main__':
    import automol
    GRA = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
            3: ('C', 1, None), 4: ('C', 1, None), 5: ('O', 1, None),
            6: ('O', 1, None), 7: ('O', 0, None), 8: ('O', 0, None),
            9: ('O', 0, None), 10: ('O', 0, None)},
           {frozenset({10, 4}): (1, None), frozenset({8, 2}): (1, None),
            frozenset({3, 4}): (1, None), frozenset({9, 6}): (1, None),
            frozenset({9, 3}): (1, None), frozenset({10, 7}): (1, None),
            frozenset({0, 2}): (1, None), frozenset({2, 4}): (1, None),
            frozenset({8, 5}): (1, None), frozenset({1, 3}): (1, None)})
    GPCS = expand_stereo_with_priorities_and_amchis(GRA)
    for GRA, _, CHI in GPCS:
        print(CHI)
        CHI_ = automol.graph.amchi(GRA)
        assert CHI == CHI_

    print(len(expand_stereo(GRA, enant=True, symeq=True)))
    print(len(expand_stereo(GRA, enant=True, symeq=False)))
    print(len(expand_stereo(GRA, enant=False, symeq=True)))
    print(len(expand_stereo(GRA, enant=False, symeq=False)))
