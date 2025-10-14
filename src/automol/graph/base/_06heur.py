"""functions for working with rotational bonds and groups

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools
from typing import List, Optional

import more_itertools as mit
import numpy
from phydat import phycon

from ...util import heuristic
from ._00core import (
    atom_implicit_hydrogens,
    atom_neighbor_atom_keys,
    atom_symbols,
    atom_zmat_neighbor_atom_key,
    atoms_neighbor_atom_keys,
    bond_keys,
    explicit,
    implicit,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reactants_graph_without_stereo,
    ts_reacting_atom_keys,
    ts_transferring_atoms,
    without_dummy_atoms,
    without_reacting_bonds,
)
from ._02algo import branch_atom_keys, rings_bond_keys
from ._03kekule import (
    atom_hybridizations,
    kekules_bond_orders_collated,
    linear_segments_atom_keys,
    rigid_planar_bonds,
    sigma_radical_atom_bond_keys,
    vinyl_radical_atom_bond_keys,
)

# bond angles
TET_ANG = 109.4712  # degrees
TRI_ANG = 120.0  # degrees
LIN_ANG = 180.0  # degrees


# heuristic coordinate values
def heuristic_bond_distance(
    gra,
    key1: int,
    key2: int,
    fdist_factor: float = 1.1,
    bdist_factor: float = 0.9,
    angstrom: bool = True,
    check: bool = False,
) -> float:
    """The heuristic bond distance between two bonded atoms

    For non-reacting atoms, returns whichever is smaller of (a.) the sum of covalent
    radii, and (b.) the average vdw radius.

    For reacting atoms, returns a multiple of whichever is larger of the same distances;
    The multiple for forming bonds is set by `fdist_factor` and that for breaking bonds
    is set by `bdist_factor`

    :param gra: Molecular graph
    :type gra: automol graph data structure
    :param key1: The first atom key
    :type key1: int
    :param key2: The second atom key
    :type key2: int
    :param fdist_factor: Set the forming bond distance to this times the average
        van der Waals radius, defaults to 1.1
    :type fdist_factor: float, optional
    :param bdist_factor: Set the breaking bond distance to this times the average
        van der Waals radius, defaults to 0.9
    :param angstrom: Return in angstroms intead of bohr?, defaults to True
    :type angstrom: bool, optional
    :param check: Check that these atoms are in fact bonded, defaults to False
    :type check: bool, optional
    :return: The heuristic distance
    :rtype: float
    """
    if check:
        assert key1 in atoms_neighbor_atom_keys(gra)[key2]

    symb_dct = atom_symbols(gra)
    symb1, symb2 = map(symb_dct.__getitem__, [key1, key2])

    bkey = frozenset({key1, key2})
    if bkey in ts_forming_bond_keys(gra):
        dist = fdist_factor * heuristic.bond_distance_limit(
            symb1, symb2, angstrom=angstrom
        )
    elif bkey in ts_breaking_bond_keys(gra):
        dist = bdist_factor * heuristic.bond_distance_limit(
            symb1, symb2, angstrom=angstrom
        )
    else:
        dist = heuristic.bond_distance(symb1, symb2, angstrom=angstrom)

    return dist


def heuristic_bond_distance_limit(
    gra,
    key1: int,
    key2: int,
    dist_factor: float = None,
    angstrom: bool = True,
) -> float:
    """The heuristic bond distance between two bonded atoms

    Returns `dist_factor` times whichever is larger of of (a.) the sum of covalent
    radii, and (b.) the average vdw radius.

    :param gra: Molecular graph
    :type gra: automol graph data structure
    :param key1: The first atom key
    :type key1: int
    :param key2: The second atom key
    :type key2: int
    :param dist_factor: The multiplier on the distance limit, defaults to None
    :type dist_factor: float, optional
    :param angstrom: Return in angstroms intead of bohr?, defaults to True
    :type angstrom: bool, optional
    :return: The heuristic bond distance limit
    :rtype: float
    """
    symb_dct = atom_symbols(gra)
    symb1, symb2 = map(symb_dct.__getitem__, [key1, key2])
    return heuristic.bond_distance_limit(symb1, symb2, dist_factor, angstrom=angstrom)


def heuristic_bond_angle(
    gra, key1, key2, key3, degree=False, check=False, hyb_dct=None
):
    """heuristic bond angle

    If being reused multiple times, you can speed this up by passing in the
    hybridizations, so they don't need to be recalculated
    """
    if check:
        assert {key1, key3} <= set(atoms_neighbor_atom_keys(gra)[key2])

    if hyb_dct is None:
        hyb_dct = atom_hybridizations(gra)

    hyb2 = hyb_dct[key2]
    if hyb2 == 3:
        ang = TET_ANG
    elif hyb2 == 2:
        ang = TRI_ANG
    else:
        assert hyb2 == 1
        ang = LIN_ANG

    ang *= 1 if degree else phycon.DEG2RAD

    return ang


# heuristic structural properties
def rotational_bond_keys(
    gra,
    lin_keys: Optional[List[int]] = None,
    with_h_rotors: bool = True,
    with_ch_rotors: bool = True,
    with_rings_rotors: bool = False,
):
    """Get all rotational bonds for a graph

    For TS graphs, this will include only bonds which are not pi bonds for
    *either* reactants or products.

    :param gra: the graph
    :param lin_keys: keys to linear atoms in the graph
    :type lin_keys: list[int]
    :param with_h_rotors: Include rotors neighbored by only hydrogens on one side?
    :param with_ch_rotors: Include rotors with C neighbored by only hydrogens?
    :param with_rings_rotors: Include rotors in rings?
    :returns: The rotational bond keys
    :rtype: frozenset[frozenset[{int, int}]]
    """
    rot_skeys_lst = rotational_segment_keys(
        gra,
        lin_keys=lin_keys,
        with_h_rotors=with_h_rotors,
        with_ch_rotors=with_ch_rotors,
        with_rings_rotors=with_rings_rotors,
        extend_lin_seg=True,
    )
    rot_bkeys = [frozenset(ks[-2:]) for ks in rot_skeys_lst]
    rot_bkeys = frozenset(sorted(rot_bkeys, key=sorted))
    return rot_bkeys


def rotational_segment_keys(
    gra,
    lin_keys: Optional[List[int]] = None,
    with_h_rotors: bool = True,
    with_ch_rotors: bool = True,
    with_rings_rotors: bool = False,
    extend_lin_seg: bool = True,
):
    """Get the keys for all rotational segments (bonds or linear segments)

    For TS graphs, this will include only bonds which are not pi bonds for
    *either* reactants or products.

    :param gra: the graph
    :param lin_keys: keys to linear atoms in the graph
    :type lin_keys: list[int]
    :param with_h_rotors: Include rotors neighbored by only hydrogens on one side?
    :param with_ch_rotors: Include rotors with C neighbored by only hydrogens?
    :param with_rings_rotors: Include rotors in rings?
    :returns: The rotational bond keys
    :rtype: frozenset[frozenset[{int, int}]]
    """
    gra = explicit(gra)
    sym_dct = atom_symbols(gra)
    nkeys_dct = atoms_neighbor_atom_keys(gra)
    bord_dct = kekules_bond_orders_collated(gra)
    rng_bkeys = list(itertools.chain(*rings_bond_keys(gra)))

    def _is_rotational_bond(bkey):
        """Not guaranteed to have out-of-line neighbors

        This is taken care of below by subtracting bonds in linear segments
        """
        ngb_keys_lst = [nkeys_dct[k] - bkey for k in bkey]

        is_single = max(bord_dct[bkey]) <= 1
        has_neighbors = all(ngb_keys_lst)
        not_in_ring = bkey not in rng_bkeys

        is_h_rotor = any(
            set(map(sym_dct.__getitem__, ks)) == {"H"} for ks in ngb_keys_lst
        )
        is_chx_rotor = is_h_rotor and any(sym_dct[k] == "C" for k in bkey)

        return (
            is_single
            and has_neighbors
            and (not_in_ring or with_rings_rotors)
            and (not is_h_rotor or with_h_rotors)
            and (not is_chx_rotor or with_ch_rotors)
        )

    # 1. Find the rotational bonds
    rot_bkeys = frozenset(filter(_is_rotational_bond, bond_keys(gra)))

    # 2. Find the linear segments, extended to include in-line neighbors
    lin_seg_keys_lst = linear_segments_atom_keys(
        gra, lin_keys=lin_keys, extend=extend_lin_seg
    )

    # 3. Start the rotational segment key list with linear segments that can rotate
    keys_lst = []
    for seg_keys in lin_seg_keys_lst:
        seg_bkeys = list(map(frozenset, mit.pairwise(seg_keys)))
        end_bkeys = {seg_bkeys[0], seg_bkeys[-1]}
        if end_bkeys <= rot_bkeys:
            keys_lst.append(seg_keys)

        # (Remove bonds on either end from list of rotational bonds)
        rot_bkeys -= set(seg_bkeys)

    # 4. Add the remaining rotational bonds to the list
    keys_lst.extend(map(sorted, rot_bkeys))

    keys_lst = frozenset(map(tuple, sorted(keys_lst, key=sorted)))
    return keys_lst


def rotational_coordinates(
    gra,
    segment: bool = True,
    lin_keys: Optional[List[int]] = None,
    with_h_rotors: bool = True,
    with_ch_rotors: bool = True,
    extend_lin_seg: bool = True,
):
    """Get torsion coordinates for rotational segments

    For rotational linear segments, the coordinate is either based on the ends of the
    segment (segment=True), or based on the final bond in the segment, in which case
    there must be a dummy atom for this to work.

    Note that only the latter case will be a valid z-matrix coordinate.

    :param gra: the graph
    :type gra: automol graph data structure
    :param segment: Stretch the coordinates across linear segments, instead of using the
        final bond with a dummy atom? defaults to True
    :type segment: bool, optional
    :param lin_keys: keys to linear atoms in the graph
    :type lin_keys: list[int]
    :param with_h_rotors: Include H rotors?
    :type with_h_rotors: bool
    :param with_ch_rotors: Include CH rotors?
    :type with_ch_rotors: bool
    :returns: The rotational bond keys
    :rtype: frozenset[frozenset[{int, int}]]
    """
    skeys_lst = rotational_segment_keys(
        gra,
        lin_keys=lin_keys,
        with_h_rotors=with_h_rotors,
        with_ch_rotors=with_ch_rotors,
        extend_lin_seg=extend_lin_seg,
    )

    coo_keys = []
    for skeys in skeys_lst:
        end_key1 = skeys[0] if segment else skeys[-2]
        end_key2 = skeys[-1]

        end_nkey1 = atom_zmat_neighbor_atom_key(gra, end_key1, excl_keys=skeys)

        excl_keys = list(skeys) + [end_nkey1]
        end_nkey2 = atom_zmat_neighbor_atom_key(gra, end_key2, excl_keys=excl_keys)

        assert end_nkey1 is not None, f"Missing dummy atom for key {end_key1}?\n{gra}"
        assert end_nkey2 is not None, f"Missing dummy atom for key {end_key2}?\n{gra}"

        coo_keys.append((end_nkey1, end_key1, end_key2, end_nkey2))

    return frozenset(coo_keys)


def rotational_groups(gra, key1, key2, dummy=False):
    """get the rotational groups for a given rotational axis

    :param gra: the graph
    :param key1: the first atom key
    :param key2: the second atom key
    """

    if not dummy:
        gra = without_dummy_atoms(gra)

    grp1 = branch_atom_keys(gra, key2, key1) - {key1}
    grp2 = branch_atom_keys(gra, key1, key2) - {key2}
    grp1 = tuple(sorted(grp1))
    grp2 = tuple(sorted(grp2))
    return grp1, grp2


def rotational_symmetry_number(gra, key1, key2, lin_keys=None):
    """Get the rotational symmetry number along a given rotational axis.

    :param gra: the graph
    :param key1: the first atom key
    :param key2: the second atom key
    """
    # Gather extended linear segments, which include "caps" (in-line neighbors), but
    # drop the ones that have their rotational symmetry broken by an out-of-line
    # breaking/forming bond. If both atoms are in one of these extended linear segments,
    # its caps should be used for symmetry number determination.
    lin_keys_lst = linear_segments_atom_keys(gra, lin_keys=lin_keys, extend=True)
    ts_bkeys = list(map(set, ts_forming_bond_keys(gra) | ts_breaking_bond_keys(gra)))
    lin_keys_lst = [
        ks
        for ks in lin_keys_lst
        if all(len(set.difference(bk, ks)) in (0, 2) for bk in ts_bkeys)
    ]

    # Replace the keys with linear segment caps if both are in the segment
    axis_keys = {key1, key2}
    for lin_keys_ in lin_keys_lst:
        if key1 in lin_keys_ and key2 in lin_keys_:
            key1, *_, key2 = lin_keys_
            axis_keys |= set(lin_keys_)
            break

    # Determine the symmetry number
    nkeys_dct = atoms_neighbor_atom_keys(without_dummy_atoms(gra))
    imp_hyd_dct = atom_implicit_hydrogens(implicit(gra))
    sym_num = 1
    for key in (key1, key2):
        if key in imp_hyd_dct:
            nkeys = nkeys_dct[key] - axis_keys
            if len(nkeys) == imp_hyd_dct[key] == 3:
                sym_num = 3
                break
    return sym_num


def ts_reacting_atom_plane_keys(tsg, key: int, include_self: bool = True):
    """Keys used to define a plane for forming the TS geometry

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param key: The key of a bond-forming atom
    :type key: int
    :param include_self: Whether to include the key itself; defaults to `True`
    :type include_self: bool, optional
    """
    nrbs_gra = without_reacting_bonds(tsg)
    rcts_gra = ts_reactants_graph_without_stereo(tsg)

    nkeys_rct = atom_neighbor_atom_keys(rcts_gra, key)
    nkeys_nrb = atom_neighbor_atom_keys(nrbs_gra, key)

    pkeys = {key} if include_self else set()
    pkeys |= nkeys_nrb if len(nkeys_rct) > 3 else nkeys_rct

    rp_dct = rigid_planar_bonds(rcts_gra, min_ncount=0, min_ring_size=0, strict=False)
    rp_bkey = next((bk for bk in rp_dct if key in bk), None)
    if rp_bkey is not None:
        (key_,) = rp_bkey - {key}
        nkeys_ = rp_dct[rp_bkey][sorted(rp_bkey).index(key_)]
        pkeys |= {key_}
        pkeys |= set(nkeys_)

    return frozenset(pkeys)


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
    rcts_gra = ts_reactants_graph_without_stereo(tsg)
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
