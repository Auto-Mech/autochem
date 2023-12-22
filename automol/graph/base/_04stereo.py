"""low-level stereochemistry functions
"""
from typing import Dict, Tuple, Union

from automol.graph.base._00core import (
    AtomKeys,
    CenterKey,
    is_ts_graph,
    tetrahedral_atom_keys,
    tetrahedral_atoms,
)
from automol.graph.base._03kekule import rigid_planar_bond_keys, rigid_planar_bonds
from automol.util import dict_

CenterNeighborDict = Dict[CenterKey, Union[AtomKeys, Tuple[AtomKeys, AtomKeys]]]


def stereocenter_candidates(
    gra, atom: bool = True, bond: bool = True
) -> CenterNeighborDict:
    """Get keys to stereocenter candidates in the graph

    Stereocenter candidates are atoms and bonds which are potentially stereogenic.
    The only thing left to check are the canonical priorities (atom symmetry classes)
    of their neighbors.

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param atom: Include atom stereocenter candidates? defaults to True
    :type atom: bool, optional
    :param bond: Include bond stereocenter candidates? defaults to True
    :type bond: bool, optional
    :returns: The keys of the stereocenter candidates
    :rtype: Dict[Union[int, frozenset[int]], tuple[tuple, tuple]]
    """
    cand_cent_dct = {}

    # 1. Atom stereocenter candidates: tetrahedral atoms
    # These need to be calculated no matter what
    cand_atm_dct = tetrahedral_atoms(gra, min_ncount=3)
    if atom:
        cand_cent_dct.update(cand_atm_dct)

    if bond:
        # 2. Bond stereocenter candidates: (a.) rigid, planar bonds, (b.) at least one
        # neighbor on each side, (c.) not locked in ring with <8 atoms
        cand_bnd_dct = rigid_planar_bonds(gra, min_ncount=1, min_ring_size=8)

        # 3. For TS graphs, remove redundant bond stereocenter candidates, where both
        # atoms are already atom stereocenter candidates
        cand_atm_keys = set(cand_atm_dct.keys())
        cand_bnd_dct = dict_.filter_by_key(
            cand_bnd_dct, lambda bk: not bk <= cand_atm_keys
        )

        cand_cent_dct.update(cand_bnd_dct)

    return cand_cent_dct


def stereocenter_candidate_keys(
    gra, atom: bool = True, bond: bool = True
) -> frozenset[Union[int, frozenset[int]]]:
    """Get keys to stereocenter candidates in the graph

    Stereocenter candidates are atoms and bonds which are potentially stereogenic.
    The only thing left to check are the canonical priorities (atom symmetry classes)
    of their neighbors.

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param atom: Include atom stereocenter candidates? defaults to True
    :type atom: bool, optional
    :param bond: Include bond stereocenter candidates? defaults to True
    :type bond: bool, optional
    :returns: The keys of the stereocenter candidates
    :rtype: frozenset[Union[int, frozenset[int]]]
    """
    keys = frozenset()

    # 1. Atom stereocenter candidates: tetrahedral atoms
    atm_keys = tetrahedral_atom_keys(gra)  # These need to be calculated no matter what
    if atom:
        keys |= atm_keys

    if bond:
        # 2. Bond stereocenter candidates: (a.) rigid, planar bonds, (b.) at least one
        # neighbor on each side, (c.) not locked in ring with <8 atoms
        bnd_keys = rigid_planar_bond_keys(gra, min_ncount=1, min_ring_size=8)

        # 3. For TS graphs, remove redundant bond stereocenter candidates, where both
        # atoms are already atom stereocenter candidates
        if is_ts_graph(gra):
            bnd_keys = frozenset({bk for bk in bnd_keys if not bk <= atm_keys})

        keys |= bnd_keys

    return keys
