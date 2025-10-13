"""low-level stereochemistry functions
"""

import itertools
import numbers
from typing import Any, Callable, Dict, Optional, Tuple, Union

import numpy

from ... import util
from ...geom import base as geom_base
from ...util import dict_
from ._00core import (
    AtomKey,
    AtomKeys,
    BondKey,
    CenterKey,
    CenterKeys,
    align_with_geometry,
    explicit,
    local_stereo_priorities,
    negate_hydrogen_stereo_priorities,
    stereo_keys,
    stereo_parities,
    tetrahedral_atoms,
    without_dummy_atoms,
)
from ._02algo import simple_paths_between_atoms
from ._03kekule import (
    rigid_planar_bonds,
    rigid_planar_bonds_with_ring_constraints,
)
from ._04class import (
    atom_transfers,
    insertions,
    substitutions,
    vinyl_addition_candidates,
)

AtomNeighborDict = Dict[AtomKey, AtomKeys]
BondNeighborDict = Dict[BondKey, Tuple[AtomKeys, AtomKeys]]
CenterNeighborDict = Dict[CenterKey, Union[AtomKeys, Tuple[AtomKeys, AtomKeys]]]
ParityEvaluator = Callable[
    [Any, CenterKeys, Dict[AtomKey, int], CenterNeighborDict, bool],
    Dict[CenterKey, int],
]


def stereocenter_candidates(
    gra, atom: bool = True, bond: bool = True, strict: bool = True
) -> CenterNeighborDict:
    """Get the stereocenter candidates in the graph, along with their stereo-determining
    neighbors, as a dictionary

    Stereocenter candidates are atoms and bonds which are potentially stereogenic.
    The only thing left to check are the canonical priorities (atom symmetry classes)
    of their neighbors.

    The neighbor keys are sorted by local priority

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param atom: Include atom stereocenter candidates? defaults to True
    :type atom: bool, optional
    :param bond: Include bond stereocenter candidates? defaults to True
    :type bond: bool, optional
    :param strict: Only include bonds that are guaranteed to be rigid?
    :type strict: bool, optional
    :returns: A mapping of candidates onto their stereo-determining neighbors
    :rtype: CenterNeighborDict
    """
    gra = without_dummy_atoms(gra)  # remove dummy atoms
    cand_dct = {}

    # 1. Atom stereocenter candidates: tetrahedral atoms
    # These need to be calculated no matter what
    cand_atm_dct = tetrahedral_atoms(gra, min_ncount=3)
    if atom:
        cand_dct.update(cand_atm_dct)

    if bond:
        # 2. Bond stereocenter candidates: (a.) rigid, planar bonds, (b.) at least one
        # neighbor on each side, (c.) not locked in ring with <8 atoms
        cand_bnd_dct = rigid_planar_bonds(
            gra, min_ncount=1, min_ring_size=8, strict=strict
        )

        # 3. For TS graphs, remove redundant bond stereocenter candidates, where both
        # atoms are already atom stereocenter candidates
        cand_atm_keys = set(cand_atm_dct.keys())
        cand_bnd_dct = dict_.filter_by_key(
            cand_bnd_dct, lambda bk: not bk <= cand_atm_keys
        )

        cand_dct.update(cand_bnd_dct)

    return cand_dct


def atom_stereodetermining_neighbor_keys(gra, key: AtomKey) -> AtomKeys:
    """Get the stereo-determining neighbor keys of an atom, sorted by local priority

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param key: The atom key
    :type key: AtomKey
    :return: The stereo-determining neighbor keys
    :rtype: AtomKeys
    """
    return stereocenter_candidates(gra, bond=False)[key]


def bond_stereodetermining_neighbor_keys(
    gra, key1: AtomKey, key2: AtomKey
) -> Tuple[AtomKeys, AtomKeys]:
    """Get the stereo-determining neighbor keys of a bond, sorted by local priority

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param key1: The first atom key
    :type key1: AtomKey
    :param key2: The second atom key
    :type key2: AtomKey
    :return: The stereo-determining neighbor keys
    :rtype: AtomKeys
    """
    keys = sorted([key1, key2])
    bkey = frozenset(keys)
    bnkeys = stereocenter_candidates(gra, strict=False)[bkey]
    nkey1s, nkey2s = (bnkeys[keys.index(k)] for k in (key1, key2))
    return nkey1s, nkey2s


def unassigned_stereocenter_keys_from_candidates(
    gra,
    cand_dct: CenterNeighborDict,
    pri_dct: Optional[Dict[AtomKey, int]] = None,
) -> CenterKeys:
    """Identify true, unassigned stereocenters from candidates, given a priority mapping

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param cand_dct: A mapping of stereocenter candidates onto their neighbor keys
    :type cand_dct: CenterNeighborDict
    :param pri_dct: An atom priority mapping
    :type pri_dct: Optional[Dict[AtomKey, int]]
    :return: The atom and bond stereocenter candidates as separate dictionaries
    :rtype: CenterKeys
    """
    pri_dct = negate_hydrogen_stereo_priorities(gra, pri_dct, with_none=True)
    ste_atm_dct, ste_bnd_dct = stereocenter_candidates_grouped(
        cand_dct, pri_dct=pri_dct, drop_nonste=True
    )
    ste_keys = frozenset(ste_atm_dct) | frozenset(ste_bnd_dct)
    ste_keys -= stereo_keys(gra)
    return ste_keys


def stereocenter_candidates_grouped(
    cand_dct: CenterNeighborDict,
    pri_dct: Optional[Dict[AtomKey, int]] = None,
    drop_nonste: bool = False,
) -> Tuple[AtomNeighborDict, BondNeighborDict]:
    """Group stereocenter candidates by type

    :param cand_dct: A mapping of stereocenter candidates onto their neighbor keys
    :type cand_dct: CenterNeighborDict
    :param pri_dct: An optional atom priority mapping for sorting the neighbor keys
    :type pri_dct: Optional[Dict[AtomKey, int]]
    :param drop_nonste: Drop non-stereogenic candidates, based on priorities?
    :type drop_nonste: bool, optional
    :return: The atom and bond stereocenter candidates as separate dictionaries
    :rtype: Tuple[AtomNeighborDict, BondNeighborDict]
    """
    assert not (
        drop_nonste and pri_dct is None
    ), "Droping non-stereogenics requires priorities"

    def _atom_is_stereogenic(nkeys):
        pris = list(map(pri_dct.get, nkeys))
        return len(set(pris)) == len(pris)

    def _bond_is_stereogenic(nkeys):
        pris = [list(map(pri_dct.get, nks)) for nks in nkeys]
        return all(len(ps) == 1 or len(set(ps)) == len(ps) for ps in pris)

    pri_ = (
        None if pri_dct is None else dict_.sort_value_(pri_dct, missing_val=-numpy.inf)
    )

    cand_atm_dct = {}
    cand_bnd_dct = {}

    for key, nkeys in cand_dct.items():
        # Atom
        if isinstance(key, numbers.Number):
            if pri_ is not None:
                nkeys = tuple(sorted(nkeys, key=pri_))
            if not drop_nonste or _atom_is_stereogenic(nkeys):
                cand_atm_dct[key] = nkeys
        # Bond
        else:
            if pri_ is not None:
                nkeys = tuple(tuple(sorted(nks, key=pri_)) for nks in nkeys)
            if not drop_nonste or _bond_is_stereogenic(nkeys):
                cand_bnd_dct[key] = nkeys

    return cand_atm_dct, cand_bnd_dct


def stereoatom_bridgehead_pairs(
    gra, cand_dct: CenterNeighborDict | None = None, migration_only: bool = False
) -> dict[tuple[AtomKey, AtomKey], tuple[AtomKeys, AtomKeys]]:
    r"""Identify pairs of interdependent bridgehead stereoatoms, if any.

    Bridgehead stereoatoms sharing the same three bridges are interdependent -- the
    configuration of one dictates that of the other (they must be opposite).

             .....
            /     \
        H--C--...--C--H
            \     /
             .....

    :param gra: A molecular graph
    :param cand_dct: A mapping of stereocenter candidates onto their neighbor keys
        (To avoid recalculating)
    :param migration_only: Only detect bridgehead pairs for atom migrations?
    :return: A dictionary mapping pairs of bridgehead atoms onto the connected neighbors
        of each, respectively
    """
    cand_dct = (
        stereocenter_candidates(gra, atom=True, bond=False)
        if cand_dct is None
        else cand_dct
    )
    cand_atm_dct, _ = stereocenter_candidates_grouped(cand_dct)

    bhp_dct = {}
    for key1, key2 in itertools.combinations(cand_atm_dct.keys(), r=2):
        paths = simple_paths_between_atoms(gra, key1, key2, overlap=False)
        conn_dct = {p[1]: p[-2] for p in paths}
        if len(set(conn_dct.keys())) == len(set(conn_dct.values())) == 3:
            conn_nkeys1, conn_nkeys2 = zip(*sorted(conn_dct.items()), strict=True)
            bhp_dct[(key1, key2)] = (conn_nkeys1, conn_nkeys2)

    if migration_only:
        mig_pairs = list(map(sorted, atom_transfers(gra).values()))
        bhp_dct = dict_.filter_by_key(bhp_dct, lambda k: sorted(k) in mig_pairs)

    return bhp_dct


# parity evaluation helpers
def geometry_atom_parity(gra, geo, key, nkeys=None, geo_idx_dct=None):
    r""" Calculate an atom parity directly from a geometry

    Neighboring atom keys (`nkeys`) must be passed in as a priority-sorted
    list. If `None`, a local parity calculation will occur based on the
    atom keys in the molecular graph.

    Atom parity is defined as follows:

    The four keys passed in are apices of a tetrahedron. Looking at 2, 3,
    and 4 from 1, they will either ascend in clockwise or counterclockwise
    order.

    If ascending in counterclockwise order, the parity is False ('-').
    If ascending in clockwise order, the parity is True ('+').

            2                   2
           /1\                 /1\
          3---4               4---3

        counterclockwise    clockwise
        False               True
        '-'                 '+'

    (Viewed looking down from 1)

    If only three keys are passed in, they will be treated as keys 2, 3,
    and 4 above and it will be assumed that there is a lone pair at 1.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param atm_key: the atom key whose parity is being evaluated
    :type atm_key: int
    :param nkeys: the neighboring atom keys, pre-sorted by priority
    :type nkeys: list[int]
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    """
    assert gra == explicit(gra), f"Implicit hydrogens are not a allowed here: {gra}"
    nkeys = atom_stereodetermining_neighbor_keys(gra, key) if nkeys is None else nkeys
    gra, geo, (key, nkeys), *_ = align_with_geometry(
        gra, geo, (key, nkeys), geo_idx_dct
    )

    # If there are only three groups, use the stereo atom itself as
    # the top apex of the tetrahedron.
    nkeys = [key, *nkeys] if len(nkeys) < 4 else nkeys

    assert len(nkeys) == 4, f"nkeys: {nkeys}, key: {key}, gra:\n{gra}"

    xyzs = geom_base.coordinates(geo, idxs=nkeys)
    det_mat = numpy.ones((4, 4))
    det_mat[:, 1:] = xyzs
    det_val = numpy.linalg.det(det_mat)
    assert det_val != 0.0  # for now, assume no four-atom planes
    par = bool(det_val > 0.0)
    return par


def geometry_bond_parity(gra, geo, bnd_key, bnd_nkeys=None, geo_idx_dct=None):
    r""" Calculate a bond parity directly from a geometry

    Neighboring bond keys (`bnd_nkeys`) must be passed in as a pair of
    priority-sorted lists corresponding to the first and second atoms in
    `bnd_key`. Note that the latter must be an *ordered list* in this case!
    If `None`, a local parity calculation will occur based on the atom keys
    in the molecular graph.

    Bond parity is defined as follows:

    For each atom in the double bond, find the heavy-atom neighbor with the
    higher canonical number. Although hydrogen atoms have higher canonical
    numbers, they are always given lowest priority.

    If the neighbors are cis to each other, the parity is False ('-').
    If the neighbors are trans to each other, the parity is True ('+').

        max     max      max     min
           \   /            \   /
            A=B              A=B
           /   \            /   \
        min     min      min     max

        cis              trans
        False            True
        '-'              '+'

    If one side only has a single neighbor, then it is compared with the
    maximum neighbor on the other side.

        max     nei      max
           \   /            \
            A=B              A=B
           /                /   \
        min              min     nei

        cis              trans
        False            True
        '-'              '+'

    If both sides have only single neighbors, then they are compared to
    each other.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param bnd_key: the bond key. If using `bnd_nkeys`, this must be an
        ordered list!
    :type bnd_key: list[int]
    :param bnd_nkeys: a pair of lists of neighboring keys for the first and
        second atoms in `bnd_key`, respectively.
    :type bnd_nkeys: list[list[int]]
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    """
    assert gra == explicit(gra), f"Implicit hydrogens are not a allowed here: {gra}"
    bnd_nkeys = (
        bond_stereodetermining_neighbor_keys(gra, *bnd_key)
        if bnd_nkeys is None
        else bnd_nkeys
    )
    gra, geo, (bnd_key, bnd_nkeys), *_ = align_with_geometry(
        gra, geo, (bnd_key, bnd_nkeys), geo_idx_dct
    )

    bnd_key = sorted(bnd_key)

    key1, key2 = bnd_key
    nkey1, nkey2 = (nks[-1] for nks in bnd_nkeys)

    xyz1, xyz2, nxyz1, nxyz2 = geom_base.coordinates(
        geo, idxs=(key1, key2, nkey1, nkey2)
    )

    # Project out the central bond direction
    bvec = numpy.subtract(xyz2, xyz1)
    nvec1 = numpy.subtract(nxyz1, xyz1)
    nvec2 = numpy.subtract(nxyz2, xyz2)
    pvec1 = util.vector.orthogonalize(bvec, nvec1)
    pvec2 = util.vector.orthogonalize(bvec, nvec2)

    dot_val = numpy.vdot(pvec1, pvec2)
    assert dot_val != 0.0  # for now, assume not collinear
    par = bool(dot_val < 0.0)
    return par


def local_flipped_parities(
    gra,
    par_dct: Dict[CenterKey, bool],
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
):
    """Get local parity assignments from canonical ones, or vice versa

    :param gra: A graph
    :type gra: automol graph data structure
    :param par_dct: The original parity assignments
    :type par_dct: Dict[CenterKey, bool]
    :param pri_dct: The canonical priority mapping
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
    :type nkeys_dct: CenterNeighborDict
    """
    # Read in the parities
    par_dct = par_dct.copy()
    keys = dict_.keys_by_value(par_dct, lambda x: x is not None)

    nkeys_dct = dict_.by_key(nkeys_dct, keys)
    loc_pri_dct = local_stereo_priorities(gra)

    # Group atoms and bonds and sort their neighbors, canonically and locally
    can_nkeys_dct, can_bnkeys_dct = stereocenter_candidates_grouped(
        nkeys_dct, pri_dct=pri_dct
    )
    loc_nkeys_dct, loc_bnkeys_dct = stereocenter_candidates_grouped(
        nkeys_dct, pri_dct=loc_pri_dct
    )

    # Determine atom parity flips
    for key, can_nkeys in can_nkeys_dct.items():
        loc_nkeys = loc_nkeys_dct[key]

        par_dct[key] ^= util.is_odd_permutation(loc_nkeys, can_nkeys)

    # Determine bond parity flips
    for bkey, (can_nkeys1, can_nkeys2) in can_bnkeys_dct.items():
        (loc_nkeys1, loc_nkeys2) = loc_bnkeys_dct[bkey]

        # Get the maximum-priority neighbor for each atom
        can_nmax1, can_nmax2, loc_nmax1, loc_nmax2 = (
            nks[-1] for nks in (can_nkeys1, can_nkeys2, loc_nkeys1, loc_nkeys2)
        )

        par_dct[bkey] ^= (loc_nmax1 != can_nmax1) ^ (loc_nmax2 != can_nmax2)

    return par_dct


def substitution_reversal_parities(
    tsg,
    par_dct: Dict[CenterKey, bool],
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
) -> Dict[AtomKey, bool]:
    r"""Update parity assignments for a substitution reversal


        Case 1: Reversal causes local parity flip (reversal value: `True`):

                  3                     3
                  |                     |
            1-----C  +  2   <=>   1  +  C-----2
                 / \                   / \
                5   4                 5   4

                clockwise             counterclockwise
                ('+')                 ('-')

        Case 2: Reversal leaves local parity unchanged (reversal value: `False`):

                  2                     2
                  |                     |
            1-----C  +  5   <=>   1  +  C-----5
                 / \                   / \
                4   3                 4   3

                clockwise             clockwise
                ('+')                 ('+')

        Given a reversal value, r_flip, the local parity upon reversing TS direction is
        given by the following:

            p_rev = p_forw ^ r_flip

    :param tsg: A TS graph
    :type tsg: automol graph data structure
    :param par_dct: The original parity assignments
    :type par_dct: Dict[CenterKey, bool]
    :param pri_dct: The priority mapping
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
    :type nkeys_dct: CenterNeighborDict
    :return: A mapping identifying which stereoatoms flip upon reversal
    :rtype: Dict[AtomKey, bool]
    """
    # Read in the parities
    par_dct = par_dct.copy()
    keys = dict_.keys_by_value(par_dct, lambda x: x is not None)

    pri_ = dict_.sort_value_(pri_dct, missing_val=-numpy.inf)

    subst_dct = substitutions(tsg)

    for key in keys:
        if key in subst_dct:
            don_key, acc_key = subst_dct[key]
            nkeys0 = sorted(nkeys_dct[key], key=pri_)
            nkeys0[nkeys0.index(don_key)] = acc_key
            nkeys = sorted(nkeys0, key=pri_)

            par_dct[key] ^= util.is_even_permutation(nkeys0, nkeys)
    return par_dct


def vinyl_addition_reactant_parities(
    tsg,
    par_dct: Dict[CenterKey, bool],
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
) -> Dict[BondKey, bool]:
    """Parity flips for the reactant parities of a vinyl addition TS graph

    :param tsg: A TS graph
    :type tsg: automol graph data structure
    :param par_dct: The original parity assignments
    :type par_dct: Dict[CenterKey, bool]
    :param pri_dct: The local priority mapping
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the *reactant* graph
    :type nkeys_dct: CenterNeighborDict
    :return: The updated parities
    :rtype: Dict[AtomKey, bool]
    """
    # Read in the parities
    par_dct = par_dct.copy()
    keys = dict_.keys_by_value(par_dct, lambda x: x is not None)

    pri_ = dict_.sort_value_(pri_dct, missing_val=-numpy.inf)

    vin_add_dct = vinyl_addition_candidates(tsg)

    for bkey in keys:
        if bkey in vin_add_dct:
            key, add_key = vin_add_dct[bkey]
            idx = sorted(bkey).index(key)

            nkeys0 = list(nkeys_dct[bkey][idx])
            nkeys = nkeys0 + [add_key]

            nmax0 = max(nkeys0, key=pri_)
            nmax = max(nkeys, key=pri_)

            par_dct[bkey] ^= nmax0 != nmax

    return par_dct


def insertion_reactant_parities(
    tsg,
    par_dct: Dict[CenterKey, bool],
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
) -> Dict[BondKey, bool]:
    """Parity flips for the reactant parities of an insertion TS graph

    Determines the insertion *bond parity* as a flip of its constituent atom parities:

        p_bond = p_atom1 ^ p_atom2 ^ ins_flip

    :param tsg: A TS graph
    :type tsg: automol graph data structure
    :param par_dct: The original parity assignments
    :type par_dct: Dict[CenterKey, bool]
    :param pri_dct: The local priority mapping
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the *reactant* graph
    :type nkeys_dct: CenterNeighborDict
    :return: The updated parities
    :rtype: Dict[AtomKey, bool]
    """
    # Read in the parities
    all_par_dct = stereo_parities(tsg)
    par_dct = par_dct.copy()
    keys = list(par_dct.keys())

    pri_ = dict_.sort_value_(pri_dct, missing_val=-numpy.inf)

    ins_dct = insertions(tsg)

    for bkey in keys:
        if bkey in ins_dct:
            key1, key2 = sorted(bkey)
            lea_key1, lea_key2, *_ = ins_dct[bkey]

            nkeys1b, nkeys2b = (sorted(nks, key=pri_) for nks in nkeys_dct[bkey])
            nkeys1 = [key2, lea_key1] + nkeys1b
            nkeys2 = [key1, lea_key2] + nkeys2b
            nkeys1_ = sorted(nkeys1, key=pri_)
            nkeys2_ = sorted(nkeys2, key=pri_)

            sgn1 = util.is_odd_permutation(nkeys1, nkeys1_)
            sgn2 = util.is_odd_permutation(nkeys2, nkeys2_)

            par1, par2 = map(all_par_dct.get, (key1, key2))
            if par1 is not None and par2 is not None:
                par_dct[bkey] = not par1 ^ par2 ^ sgn1 ^ sgn2

    return par_dct


def ring_bond_reactant_parities(
    tsg,
    par_dct: Dict[CenterKey, bool],
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
) -> Dict[BondKey, bool]:
    """Ring-constrained bond parities for the reactants of a TS graph

    :param tsg: A TS graph
    :type tsg: automol graph data structure
    :param par_dct: The original parity assignments
    :type par_dct: Dict[CenterKey, bool]
    :param pri_dct: The local priority mapping
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the *reactant* graph
    :type nkeys_dct: CenterNeighborDict
    :return: The updated parities
    :rtype: Dict[AtomKey, bool]
    """
    # Read in the parities
    par_dct = par_dct.copy()
    keys = list(par_dct.keys())

    pri_ = dict_.sort_value_(pri_dct, missing_val=-numpy.inf)

    vin_add_dct = vinyl_addition_candidates(tsg)

    # Insert vinyl addition neighbors back in, so they are captured by the ring
    # constraint
    nkeys_dct = nkeys_dct.copy()
    for bkey in list(nkeys_dct.keys()):
        if bkey in vin_add_dct:
            vin_key, vin_nkey = vin_add_dct[bkey]
            bnkeys = list(map(list, nkeys_dct[bkey]))
            bnkeys[sorted(bkey).index(vin_key)].append(vin_nkey)
            nkeys_dct[bkey] = tuple(map(tuple, bnkeys))

    rng_const_dct = rigid_planar_bonds_with_ring_constraints(tsg, nkeys_dct)

    for bkey in keys:
        if bkey in rng_const_dct:
            nkeys1, nkeys2 = nkeys_dct[bkey]
            nmax1 = max(nkeys1, key=pri_)
            nmax2 = max(nkeys2, key=pri_)

            rkey1, rkey2 = rng_const_dct[bkey]

            par_dct[bkey] = (nmax1 == rkey1) ^ (nmax2 == rkey2)

    return par_dct


# parity evaluators
def parity_evaluator_measure_from_geometry_(
    geo,
    local_stereo: bool = False,
    geo_idx_dct: Optional[Dict[AtomKey, AtomKey]] = None,
) -> ParityEvaluator:
    """Creates a parity evaluator that measures parities from a geometry

    :param geo: A geometry
    :type geo: automol geom data structure
    :param local_stereo: Return local stereo assignments? defaults to False
    :type local_stereo: bool, optional
    :param geo_idx_dct: A mapping of graph keys onto geometry indices
    :type geo_idx_dct: Optional[Dict[AtomKey, AtomKey]]
    :returns: A parity evaluator
    :rtype: ParityEvaluator
    """

    def parity_evaluator(
        gra,
        keys: CenterKeys,
        pri_dct: Dict[AtomKey, int],
        nkeys_dct: CenterNeighborDict,
        is_rev_ts: bool = False,
    ):
        """A parity evaluator that measures parities from a geometry

        :param gra: A graph
        :type gra: automol graph data structure
        :param keys: The stereocenter keys at which to evaluate the parity
        :type keys: CenterKeys
        :param pri_dct: A canonical atom priority mapping
        :type pri_dct: Dict[AtomKey, int]
        :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
        :type nkeys_dct: CenterNeighborDict
        :param is_rev_ts: Is this a reversed TS graph?, defaults to False
        :type is_rev_ts: bool, optional
        :return: A parity assignment mapping for the given stereocenter keys
        :rtype: Dict[CenterKey, bool]
        """
        pri_dct = negate_hydrogen_stereo_priorities(gra, pri_dct, with_none=True)

        # 0. Determine local priories, if requesting local stereo
        pri_dct = local_stereo_priorities(gra) if local_stereo else pri_dct

        # 1. Get sorted neighbor keys for stereocenters
        nkeys_dct = dict_.by_key(nkeys_dct, keys)
        atm_nkeys_dct, bnd_nkeys_dct = stereocenter_candidates_grouped(
            nkeys_dct, pri_dct=pri_dct
        )

        # 2. Measure parities for atoms and bonds
        atm_par_dct = {
            k: geometry_atom_parity(gra, geo, k, ns, geo_idx_dct=geo_idx_dct)
            for k, ns in atm_nkeys_dct.items()
        }
        bnd_par_dct = {
            k: geometry_bond_parity(gra, geo, sorted(k), ns, geo_idx_dct=geo_idx_dct)
            for k, ns in bnd_nkeys_dct.items()
        }

        # 3. Combine parity mappings
        par_dct = {**atm_par_dct, **bnd_par_dct}

        # 4. Apply substitution reversal flips, if this is a reverse TS graph
        if local_stereo and is_rev_ts:
            par_dct = substitution_reversal_parities(gra, par_dct, pri_dct, nkeys_dct)

        return par_dct

    return parity_evaluator


def parity_evaluator_read_from_graph(
    gra,
    keys: CenterKeys,
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
    is_rev_ts: bool = False,
) -> Dict[CenterKey, bool]:
    """A parity evaluator that reads parities from a graph

    :param gra: A graph
    :type gra: automol graph data structure
    :param keys: The stereocenter keys at which to evaluate the parity
    :type keys: CenterKeys
    :param pri_dct: Canonical priorities for the atoms
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
    :type nkeys_dct: CenterNeighborDict
    :param is_rev_ts: Is this a reversed TS graph?, defaults to False
    :type is_rev_ts: bool, optional
    :return: A mapping of these stereocenter keys onto their parity assignments
    :rtype: Dict[CenterKey, bool]
    """
    assert nkeys_dct or not nkeys_dct
    assert pri_dct or not pri_dct
    assert is_rev_ts or not is_rev_ts

    return dict_.by_key(stereo_parities(gra), keys)


def parity_evaluator_flip_from_graph(
    gra,
    keys: CenterKeys,
    pri_dct: Dict[AtomKey, int],
    nkeys_dct: CenterNeighborDict,
    is_rev_ts: bool = False,
) -> Dict[CenterKey, bool]:
    """A parity evaluator that reads parities from a graph and flips them to convert
    canonical to local parities (or vice versa)

    :param gra: A graph
    :type gra: automol graph data structure
    :param keys: The stereocenter keys at which the parity is being evaluated
    :type keys: CenterKeys
    :param pri_dct: Canonical priorities for the atoms
    :type pri_dct: Dict[AtomKey, int]
    :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
    :type nkeys_dct: CenterNeighborDict
    :param is_rev_ts: Is this a reversed TS graph?, defaults to False
    :type is_rev_ts: bool, optional
    :return: A mapping of these stereocenter keys onto their parity assignments
    :rtype: Dict[CenterKey, bool]
    """
    pri_dct = negate_hydrogen_stereo_priorities(gra, pri_dct, with_none=True)

    # 0. Read in the parities
    par_dct = dict_.by_key(stereo_parities(gra), keys)

    # 1. Apply local parity flips to convert to/from canonical stereo
    par_dct = local_flipped_parities(gra, par_dct, pri_dct, nkeys_dct)

    # 2. Apply substitution reversal flips, if this is a reverse TS graph
    if is_rev_ts:
        par_dct = substitution_reversal_parities(gra, par_dct, pri_dct, nkeys_dct)

    return par_dct


def parity_evaluator_reactants_from_local_ts_graph_(
    loc_tsg,
    local_stereo: bool = False,
) -> ParityEvaluator:
    """Creates a parity evaluator that evaluates reactant parities from a local TS graph

    :param loc_tsg: A TS graph with local parities
    :type loc_tsg: automol graph data structure
    :param local_stereo: Return local stereo assignments? defaults to False
    :type local_stereo: bool, optional
    :returns: A parity evaluator
    :rtype: ParityEvaluator
    """
    ts_par_dct = stereo_parities(loc_tsg)

    # Determine local TS stereo parities
    loc_pri_dct = local_stereo_priorities(loc_tsg)

    def parity_evaluator(
        gra,
        keys: CenterKeys,
        pri_dct: Dict[AtomKey, int],
        nkeys_dct: CenterNeighborDict,
        is_rev_ts: bool = False,
    ):
        """A parity evaluator that evaluates reactant parities from a local TS graph

        :param gra: A graph
        :type gra: automol graph data structure
        :param keys: The stereocenter keys at which to evaluate the parity
        :type keys: CenterKeys
        :param pri_dct: A canonical atom priority mapping
        :type pri_dct: Dict[AtomKey, int]
        :param nkeys_dct: Stereo-determining neighbors for the stereocenter keys
        :type nkeys_dct: CenterNeighborDict
        :param is_rev_ts: Is this a reversed TS graph?, defaults to False
        :type is_rev_ts: bool, optional
        :return: A parity assignment mapping for the given stereocenter keys
        :rtype: Dict[CenterKey, bool]
        """
        pri_dct = negate_hydrogen_stereo_priorities(gra, pri_dct, with_none=True)

        assert gra or not gra
        assert is_rev_ts or not is_rev_ts

        # 0. Read in local TS parities for these keys
        par_dct = dict_.by_key(ts_par_dct, keys)

        # 1. Determine constrained-ring bond parities
        par_dct = ring_bond_reactant_parities(loc_tsg, par_dct, loc_pri_dct, nkeys_dct)

        # 2. Correct vinyl addition parity flips
        par_dct = vinyl_addition_reactant_parities(
            loc_tsg, par_dct, loc_pri_dct, nkeys_dct
        )

        # 3. Determine insertion bonds from atom parities with insertion parity flips
        par_dct = insertion_reactant_parities(loc_tsg, par_dct, loc_pri_dct, nkeys_dct)

        # 4. Apply local parity flips to convert to canonical stereo, if requested
        if not local_stereo:
            par_dct = local_flipped_parities(gra, par_dct, pri_dct, nkeys_dct)

        return par_dct

    return parity_evaluator
