""" graph functions that depend on stereo assignments

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools

import numpy

import automol.amchi.base  # !!!!
import automol.geom.base  # !!!!
from automol.graph.base._algo import (
    connected_components,
)
from automol.graph.base._amchi import connected_amchi_with_indices
from automol.graph.base._canon import (
    is_canonical_enantiomer,
    parity_evaluator_from_geometry_,
    refine_priorities,
    reflect_local_stereo,
    stereo_assignment_representation,
    stereogenic_atom_keys,
    stereogenic_bond_keys,
    stereogenic_keys,
    stereogenic_keys_from_priorities,
    to_local_stereo,
)
from automol.graph.base._core import (
    atom_keys,
    atom_stereo_keys,
    atom_stereo_parities,
    bond_stereo_keys,
    bond_stereo_parities,
    frozen,
    has_stereo,
    implicit,
    is_ts_graph,
    relabel,
    set_atom_stereo_parities,
    set_bond_stereo_parities,
    set_stereo_parities,
    stereo_parities,
    without_stereo,
)
from automol.graph.base._geom import (
    geometry_atom_parity,
    geometry_bond_parity,
    geometry_pseudorotate_atom,
    geometry_rotate_bond,
)


# # core functions
def expand_stereo(gra, enant=True, symeq=False):
    """Obtain all possible stereoisomers of a graph, ignoring its assignments

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param enant: Include all enantiomers, or only canonical ones?
    :type enant: bool
    :param symeq: Include symmetrically equivalent stereoisomers?
    :type symeq: bool
    :returns: a series of molecular graphs for the stereoisomers
    """

    bools = (False, True)

    gra0 = without_stereo(gra)
    gps0 = None
    gps = [(gra0, None)]

    # 1. Expand all possible stereoisomers, along with their priority mappings
    while gps0 != gps:
        gps0 = gps
        gps = []
        seen_reps = []

        for gra1, pri_dct in gps0:
            # a. Refine priorities based on current assignments
            pri_dct = refine_priorities(gra1, pri_dct=pri_dct)

            # b. Generate a representation of the current assignments
            rep = stereo_assignment_representation(gra1, pri_dct)
            #   i. If the representation has been seen, continue (skip)
            if rep in seen_reps and not symeq:
                continue
            #  ii. If not, add it to the list of seen representations
            seen_reps.append(rep)

            # c. Find stereogenic atoms and bonds based on current priorities
            keys = stereogenic_keys_from_priorities(gra1, pri_dct)

            # d. Assign True/False parities in all possible ways
            for pars in itertools.product(bools, repeat=len(keys)):
                gra2 = set_stereo_parities(gra1, dict(zip(keys, pars)))
                gps.append((gra2, pri_dct))

    # 2. If requested, filter out non-canonical enantiomers
    if not enant:
        # Save copies of graphs, priority dictionaries, and local stereo graphs
        # that have been seen. The last will be used to check for enantiomers.
        seen_gpls = []
        for ugra, upri_dct in gps.copy():
            # a. Convert to local stereo assignments
            uloc_gra = to_local_stereo(ugra, pri_dct=upri_dct)

            # b. Reflect local stereo assignments
            rloc_gra = reflect_local_stereo(uloc_gra)

            # c. Check if the reflection matches anything we've seen
            rgra, rpri_dct = next(
                ((g, p) for g, p, l in seen_gpls if l == rloc_gra), (None, None)
            )

            # d. Add the current stereoisomer to the list of seen ones
            seen_gpls.append((ugra, upri_dct, uloc_gra))

            # e. If it does, we have a pair of enantiomers. Identify which one
            # is non-canonical and remove it.
            if rgra is not None:
                is_can = is_canonical_enantiomer(ugra, upri_dct, rgra, rpri_dct)

                if is_can is True:
                    gps.remove((rgra, rpri_dct))
                elif is_can is False:
                    gps.remove((ugra, upri_dct))

    sgras = [sgra for sgra, _ in gps]
    sgras = tuple(sorted(sgras, key=frozen))
    return sgras


# # stereo correction
def stereo_corrected_geometry(gra, geo, geo_idx_dct=None, local_stereo=False):
    """Obtain a geometry corrected for stereo parities based on a graph

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
    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(atm_keys)} if geo_idx_dct is None else geo_idx_dct
    )
    gra = relabel(gra, geo_idx_dct)

    par_dct = stereo_parities(sgr)
    bnd_keys = bond_stereo_keys(sgr)
    atm_keys = atom_stereo_keys(sgr)

    for bnd_key in bnd_keys:
        curr_par = geometry_bond_parity(gra, geo, bnd_key)
        if curr_par != par_dct[bnd_key]:
            geo = geometry_rotate_bond(gra, geo, bnd_key, numpy.pi)

    for atm_key in atm_keys:
        curr_par = geometry_atom_parity(gra, geo, atm_key)
        if curr_par != par_dct[atm_key]:
            geo = geometry_pseudorotate_atom(gra, geo, atm_key)

            assert geo is not None, (
                f"Failed to correct the following geometry:"
                f"\ngeo:\n{geo}\ngra:\n{gra}"
            )

    return geo
