""" geometry embedding using the distance geometry algorithm

This algorithm generates approximate geometries by a randomized guess at the
distance matrix, within heuristic distance bounds, which is then converted to
an approximate geometry that satifies these distances.

Blaney, J. M.; Dixon, J. S. “Distance Geometry in Molecular Modeling”.
Reviews in Computational Chemistry; VCH: New York, 1994.

The steps in the algorithm are as follows:

    1. Generate distance bounds matrix B. (Here, B is replaced with L and U).
    2. Narrow the bounds in B by triangle smoothing.
    3. Generate a distance matrix D by uniform sampling within the bounds.
    4. Generate the metric matrix G (matrix of position vector dot products).
    5. Diagonalize G and determine the principal components.
    6. The three largest eigenvectors and eigenvalues of G can be used to
    generate x, y, z coordinates for the molecule which approximately
    correspond to the distance matrix D.
    7. Do error refinement to clean up the structure and enforce correct
    chirality.

Step 1 is performed by the function distance_bounds_matrices() below.
Steps 2-6 are performed by the function automol.embed.sample_raw_distance_coordinates()
Step 7 is performed by the function automol.embed.cleaned_up_coordinates()
"""

import itertools
from typing import Any, Optional, Tuple

import numpy

from phydat import phycon

from .. import embed, geom, zmat
from ..geom import base as geom_base
from ..util import dict_
from .base import (
    add_bonds,
    atom_hybridizations,
    atom_keys,
    atom_shortest_paths,
    atom_stereo_keys,
    atom_stereo_parities,
    atom_symbols,
    atoms_neighbor_atom_keys,
    bond_keys,
    bond_orders,
    bond_stereo_keys,
    heuristic_bond_angle,
    heuristic_bond_distance,
    heuristic_bond_distance_limit,
    is_ts_graph,
    linear_atom_keys,
    rigid_planar_bond_keys,
    rings_atom_keys,
    rotational_bond_keys,
    set_stereo_from_geometry,
    stereo_parities,
    stereocenter_candidates,
    string,
    to_local_stereo,
    ts,
    without_bonds_by_orders,
    without_dummy_atoms,
    without_stereo,
)


# # geometry embedding functions
def clean_geometry(
    gra: Any,
    geo: Any,
    stereo: bool = True,
    local_stereo: bool = False,
    none_if_failed: bool = True,
    geos: Tuple[Any, ...] = (),
    geos_keys: Tuple[Tuple, ...] = (),
    relax_angles: bool = False,
    log: bool = False,
) -> Any:
    """Clean up a geometry based on this graph, removing any bonds that
    aren't supposed to be there

    :param gra: A molecular graph
    :param geo: A molecular geometry
    :param stereo: Take stereochemistry into consideration? defaults to True
    :param local_stereo: Does the graph have local stereo assignments? defaults to False
    :param relax_angles: Allow angles to relax?, defaults to False
    :param none_if_failed: Return `None` if the geometry doesn't match? defaults to True
    :param geos: Geometries for one or more subgraphs of `gra`
    :param geos_keys: Graph keys for the geometries in `geos`
    :param relax_angles: Relax the angles in `geos`?
    :returns: The cleaned-up geometry
    """
    if geo is None:
        return None

    hard_geo = hardcoded_geometry(gra)
    if hard_geo is not None:
        return hard_geo

    gra = gra if local_stereo else to_local_stereo(gra)
    rgra = ts.reactants_graph(gra) if is_ts_graph(gra) else gra

    xmat = geom_base.coordinates(geo, angstrom=True)

    # Only enforce planarity based on the reactants
    pla_dct = planarity_constraint_bounds(rgra)
    chi_dct = chirality_constraint_bounds(gra) if stereo else {}

    # If no auxiliary geometries were passed in, use the reference geometry
    if not geos:
        geos = [geo]
        geos_keys = [list(range(geom.count(geo)))]

    lmat, umat = distance_bounds_matrices(
        gra=gra,
        geos=geos,
        geos_keys=geos_keys,
        relax_angles=relax_angles,
        relax_torsions=True,
        angstrom=True,
    )

    xmat, conv = embed.cleaned_up_coordinates(
        xmat, lmat, umat, chi_dct=chi_dct, pla_dct=pla_dct
    )

    if log:
        print("Converged!" if conv else "Did not converge.")

    syms = geom_base.symbols(geo)
    xyzs = xmat[:, :3]
    geo = geom_base.from_data(syms, xyzs, angstrom=True)

    # If the clean-up failed, return `None`
    if none_if_failed and not geometry_matches(
        gra, geo, stereo=stereo, local_stereo=True, log=log
    ):
        return None

    return geo


def geometry_matches(
    gra,
    geo,
    stereo: bool = True,
    local_stereo: bool = False,
    check_ts_bonds: bool = True,
    log: bool = False,
) -> bool:
    """Check whether a geometry matches the graph

    :param gra: molecular graph with stereo parities
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: Take stereochemistry into consideration? defaults to True
    :param local_stereo: Does the graph have local stereo assignments? defaults to False
    :param check_ts_bonds: Check reacting bonds for TS graphs? If `True`, this will
        check that the atoms in reacting bonds are the closest atoms to each other
    :param log: Log information to the screen? defaults to False
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    if geo is None:
        return False

    cgra = without_stereo(gra)
    cgra0 = geom.graph_without_stereo(geo, fix_hyper=False)

    cgra_ = cgra0
    if is_ts_graph(cgra):
        r_bkeys = ts.reacting_bond_keys(cgra)
        r_ord_dct = dict_.by_key(bond_orders(cgra), r_bkeys)
        cgra_ = add_bonds(cgra_, keys=r_bkeys, ord_dct=r_ord_dct, check=False)

    # 1. Check that the connectivities match
    # Note: Do not remove dummy atoms up front! The graph and geometry must be
    # consistent for the stereo assessment
    matches = without_dummy_atoms(cgra) == without_dummy_atoms(cgra_)
    if log:
        print("Connectivity matches?", matches)
        print(f"{cgra}\n == ? \n{cgra_}")

    # 2. Check that the TS reacting bonds look reasonable
    if is_ts_graph(cgra) and check_ts_bonds:
        brk_bkeys = ts.breaking_bond_keys(cgra)
        frm_bkeys = ts.forming_bond_keys(cgra)

        # 1. Check that the breaking bond pairs are bonded in the geometry
        bkeys = bond_keys(cgra0)
        matches &= all(bk in bkeys for bk in brk_bkeys)

        # 2. Check that the forming bond pairs form closest non-bonded pairs
        no_fb_cgra_ = without_bonds_by_orders(cgra0, [0.1])
        bmatches = [
            geom.could_be_forming_bond(geo, *b, gra=no_fb_cgra_) for b in frm_bkeys
        ]
        matches &= all(bmatches)

        if log:
            print("TS forming bonds match?", dict(zip(frm_bkeys, bmatches)))

    # 3. Check that the stereo parities match
    if stereo:
        gra = gra if local_stereo else to_local_stereo(gra)
        gra_ = set_stereo_from_geometry(cgra_, geo, local_stereo=True)

        ste_akeys = atom_stereo_keys(gra)
        ste_bkeys = bond_stereo_keys(gra)

        # Exclude bonds that are likely to be near-linear
        lin_keys = linear_atom_keys(gra)
        ste_bkeys = {bk for bk in ste_bkeys if not bk & lin_keys}

        # Exclude bonds in small rings -- if marked as stereogenic rather than excluded,
        # these are ring-opening TS graphs and the bond has ambiguous stereochemistry
        rng_keys = [set(ks) for ks in rings_atom_keys(gra) if len(ks) < 8]
        ste_bkeys = {
            bk for bk in ste_bkeys if not any(bk & ks == bk for ks in rng_keys)
        }

        ste_keys = sorted(ste_akeys) + sorted(ste_bkeys)
        pars = dict_.values_by_key(stereo_parities(gra), ste_keys)
        pars_ = dict_.values_by_key(stereo_parities(gra_), ste_keys)
        matches &= pars == pars_

        if log:
            print(f"Stereo parities match at {ste_keys}?\n{pars} ==? {pars_}\n")

    return matches


def zmatrix_matches(
    gra,
    zma,
    stereo: bool = True,
    local_stereo: bool = False,
    check_ts_bonds: bool = True,
    log: bool = False,
) -> bool:
    """Check whether a z-matrix matches the graph

    :param gra: molecular graph with stereo parities
    :type gra: automol graph data structure
    :param zma: molecular z-matrix
    :type zma: automol z-matrix data structure
    :param stereo: Take stereochemistry into consideration? defaults to True
    :type stereo: bool, optional
    :param local_stereo: Does the graph have local stereo assignments? defaults to False
    :type local_stereo: bool, optional
    :param check_ts_bonds: Check reacting bonds for TS graphs? If `True`, this will
        check that the atoms in reacting bonds are the closest atoms to each other
    :type check_ts_bonds: bool, optional
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    geo = zmat.geometry(zma, dummy=True)
    return geometry_matches(
        gra,
        geo,
        stereo=stereo,
        local_stereo=local_stereo,
        check_ts_bonds=check_ts_bonds,
        log=log,
    )


# # bounds matrices
def distance_bounds_matrices(
    gra,
    geos: Tuple[Any, ...] = (),
    geos_keys: Tuple[Tuple, ...] = (),
    relax_angles: bool = False,
    relax_torsions: bool = False,
    angstrom: bool = True,
):
    """generates initial distance bounds matrices for various different
    scenarios, allowing the geometry to be manipulated in different ways

    :param gra: A molecular graph
    :param geos: Geometries for one or more subgraphs of `gra`
    :param geos_keys: Graph keys for the geometries in `geos`
    :param relax_angles: Relax the angles in `geos`?
    :param relax_torsions: Relax the torsions in `geos`?
    :param angstrom: Return units of angstroms?
    """
    keys = sorted(atom_keys(gra))
    rgra = ts.reactants_graph(gra) if is_ts_graph(gra) else gra
    sp_dct = atom_shortest_paths(rgra)

    natms = len(keys)

    lmat, umat = _distance_bounds_matrices(rgra, sp_dct=sp_dct)

    # save the current values so that we can overwrite the fixed torsions below
    lmat_old = numpy.copy(lmat)
    umat_old = numpy.copy(umat)

    # 1. Set lower and upper bonds based on geometries, if requested
    for geo, geo_keys in zip(geos, geos_keys):
        xmat = geom_base.coordinates(geo, angstrom=angstrom)
        dmat = embed.distance_matrix_from_coordinates(xmat)

        lmat[numpy.ix_(geo_keys, geo_keys)] = dmat
        umat[numpy.ix_(geo_keys, geo_keys)] = dmat

    # 2. reopen bounds on the angles from the reactant
    # (also triggers torsional re-opening, for sufficient flexibility)
    if relax_angles:
        relax_torsions = True

        ang_ijs = [
            [i, j]
            for i, j in itertools.combinations(range(natms), 2)
            if j in sp_dct[i] and len(sp_dct[i][j]) == 3
        ]
        ang_ijs += list(map(list, map(reversed, ang_ijs)))

        ang_idxs = tuple(map(list, zip(*ang_ijs)))

        if ang_idxs:
            lmat[ang_idxs] = lmat_old[ang_idxs]
            umat[ang_idxs] = umat_old[ang_idxs]

    # 3. reopen bounds on the torsions from the reactant
    # (also opens >5-atom distances)
    if relax_torsions:
        rot_bkeys = rotational_bond_keys(gra)

        # open up torsions at rotational bonds
        tors_ijs = [
            [i, j]
            for i, j in itertools.combinations(range(natms), 2)
            if j in sp_dct[i]
            and len(sp_dct[i][j]) == 4
            and frozenset(sp_dct[i][j][1:3]) in rot_bkeys
        ]
        # open up >5-atom distances
        tors_ijs.extend(
            [i, j]
            for i, j in itertools.combinations(range(natms), 2)
            if j in sp_dct[i] and len(sp_dct[i][j]) > 4
        )

        tors_ijs += list(map(list, map(reversed, tors_ijs)))

        tors_idxs = tuple(map(list, zip(*tors_ijs)))

        if tors_idxs:
            lmat[tors_idxs] = lmat_old[tors_idxs]
            umat[tors_idxs] = umat_old[tors_idxs]

    # 4. set distance bounds for the forming bonds
    for bkey in ts.forming_bond_keys(gra):
        key1, key2 = sorted(bkey)
        fdist = ts.heuristic_bond_distance(gra, key1, key2, angstrom=True)
        lmat[(key1, key2)] = lmat[(key2, key1)] = fdist
        umat[(key1, key2)] = umat[(key2, key1)] = fdist

    return lmat, umat


def _distance_bounds_matrices(gra, sp_dct=None):
    """initial distance bounds matrices

    :param gra: molecular graph
    :param keys: atom keys specifying the order of indices in the matrix
    :param sp_dct: a 2d dictionary giving the shortest path between any pair of
        atoms in the graph
    """
    keys = sorted(atom_keys(gra))
    sp_dct = atom_shortest_paths(gra) if sp_dct is None else sp_dct

    bounds_ = path_distance_bounds_(gra)

    natms = len(keys)
    umat = numpy.zeros((natms, natms))
    lmat = numpy.zeros((natms, natms))
    for (idx1, key1), (idx2, key2) in itertools.combinations(enumerate(keys), 2):
        if key2 in sp_dct[key1]:
            path = sp_dct[key1][key2]
            ldist, udist = bounds_(path)
            lmat[idx1, idx2] = lmat[idx2, idx1] = ldist
            umat[idx1, idx2] = umat[idx2, idx1] = udist
        else:
            # they are disconnected
            lmat[idx1, idx2] = lmat[idx2, idx1] = closest_approach(gra, key1, key2)
            umat[idx1, idx2] = umat[idx2, idx1] = 999

        assert lmat[idx1, idx2] <= umat[idx1, idx2], (
            "Lower bound exceeds upper bound. This is a bug!\n"
            f"{string(gra)}\npath: {str(path)}\n"
        )

    return lmat, umat


# # constraint dictionaries
def chirality_constraint_bounds(loc_gra):
    """bounds for enforcing chirality restrictions"""
    keys = sorted(atom_keys(loc_gra))
    ste_keys = set(atom_stereo_keys(loc_gra)) & set(keys)
    par_dct = atom_stereo_parities(loc_gra)
    nkeys_dct = stereocenter_candidates(loc_gra, atom=True, bond=False)

    def _chirality_constraint(key):
        nkeys = nkeys_dct[key]
        # Add the key itself as apex, if we only have 3 neighbors
        if len(nkeys) < 4:
            nkeys = (key,) + nkeys
        idxs = tuple(map(keys.index, nkeys))
        vol_range = (-999.0, -1.0) if par_dct[key] else (+1.0, +999.0)
        return idxs, vol_range

    chi_dct = dict(map(_chirality_constraint, ste_keys))
    return chi_dct


def planarity_constraint_bounds(gra):
    """bounds for enforcing planarity restrictions"""
    keys = sorted(atom_keys(gra))
    nkeys_dct = atoms_neighbor_atom_keys(gra)
    bkeys = rigid_planar_bond_keys(gra, min_ncount=0, min_ring_size=0)

    def _planarity_constraints(bkey):
        key1, key2 = sorted(bkey)
        nkey1s = sorted(nkeys_dct[key1] - {key2})
        nkey2s = sorted(nkeys_dct[key2] - {key1})

        lst = []

        # I don't think the order of the keys matters, but I tried to be
        # roughly consistent with Figure 8 in the Blaney Dixon paper
        if len(nkey1s) == 2 and len(nkey2s) == 2:
            lst.append(tuple(map(keys.index, nkey1s + nkey2s)))
        if len(nkey1s) == 2:
            lst.append(tuple(map(keys.index, [key1, key2] + nkey1s)))
        if len(nkey2s) == 2:
            lst.append(tuple(map(keys.index, [key1, key2] + nkey2s)))
        if (len(nkey1s) == 2 and len(nkey2s) == 1) or (
            len(nkey1s) == 1 and len(nkey2s) == 2
        ):
            lst.append(tuple(map(keys.index, [key1] + nkey1s + nkey2s)))
            lst.append(tuple(map(keys.index, [key2] + nkey1s + nkey2s)))
        if len(nkey1s) == 1 and len(nkey2s) == 1:
            lst.append(tuple(map(keys.index, [key1, key2] + nkey1s + nkey2s)))

        return tuple(lst)

    const_dct = {
        idxs: (-0.5, +0.5)
        for idxs in itertools.chain(*map(_planarity_constraints, bkeys))
    }

    return const_dct


# # heuristic coordinate values
def heuristic_bond_angle_distance(
    gra,
    key1,
    key2,
    key3,
    d12=None,
    d23=None,
    angstrom=True,
    a123=None,
    degree=True,
    hyb_dct=None,
):
    """heuristic distance between atoms at two ends of a bond angle

    :param angstrom: whether or not to return the distance in angstroms
    :type angstrom: bool
    :param a123: (optional) specify the value of the angle
    :type a123: float
    :param degree: units for the angle, if specified
    :type degree: bool

    uses the law of cosines:

        d13 = sqrt(d12^2 + d23^2 - 2*d12*d23*cos(a123))
    """
    if a123 is None:
        a123 = heuristic_bond_angle(
            gra, key1, key2, key3, degree=False, hyb_dct=hyb_dct
        )
    else:
        a123 = a123 * phycon.DEG2RAD if degree else a123

    d12 = heuristic_bond_distance(gra, key1, key2, angstrom=angstrom)
    d23 = heuristic_bond_distance(gra, key2, key3, angstrom=angstrom)
    d13 = numpy.sqrt(d12**2 + d23**2 - 2 * d12 * d23 * numpy.cos(a123))
    return d13


def closest_approach(gra, key1, key2):
    """closest approach between atoms, based on their van der Waals radii

    Warning: The scaling factor on the van der waals radii was arbitrarily
    chosen based on limited tests and may need to be lowered
    """
    return 1.5 * heuristic_bond_distance_limit(gra, key1, key2, angstrom=True)


def path_distance_bounds_(gra):
    """upper distance bound between two ends of a path

    :param gra: molecular graph
    :param path: the shortest path between two atoms
    :type path: list or tuple
    """

    hyb_dct = atom_hybridizations(gra)
    rng_keys_lst = rings_atom_keys(gra)

    def _distance_bounds(path):
        # if the path is 0, the atoms are disconnected and could be arbitrarily
        # far apart
        rsz = shared_ring_size(path, rng_keys_lst)
        if len(path) == 1:
            ldist = udist = 0
        elif len(path) == 2:
            ldist = udist = heuristic_bond_distance(gra, *path)
        elif len(path) == 3:
            if rsz == 0:
                ldist = udist = heuristic_bond_angle_distance(
                    gra, *path, hyb_dct=hyb_dct
                )
            else:
                a123 = (rsz - 2.0) * 180.0 / rsz
                la123 = a123 - 10.0
                ua123 = a123 + 10.0
                lrdist = heuristic_bond_angle_distance(gra, *path, a123=la123)
                urdist = heuristic_bond_angle_distance(gra, *path, a123=ua123)
                odist = heuristic_bond_angle_distance(gra, *path, hyb_dct=hyb_dct)
                ldist = min(lrdist, odist)
                udist = max(urdist, odist)
        # otherwise, just do the sum of the distances between atoms along the
        # path
        else:
            # we can't handle disconnected points, because in that case the
            # path is [] and there is no way to recover the keys
            assert len(path) > 2
            ldist = closest_approach(gra, path[0], path[-1])
            udist = 999.0

        return ldist, udist

    return _distance_bounds


# # helpers
def shared_ring_size(keys, rng_keys_lst):
    """determine whether these keys share a ring and, if so, determine the
    size
    """
    rng_keys = next(
        (rng_keys for rng_keys in rng_keys_lst if set(keys) <= set(rng_keys)), ()
    )
    natms = len(rng_keys)
    return natms


def hardcoded_geometry(gra: Any) -> Optional[Any]:
    """Generate a hardcoded geometry if this is a monatomic or diatomic molecule

    :param gra: A molecular graph
    :return: A geometry, if the graph is simple, otherwise `None`
    """
    symb_dct = atom_symbols(gra)

    # Build monatomics and diatomics directly
    if len(symb_dct) == 1:
        symbs = list(symb_dct.values())
        xyzs = [[0.0, 0.0, 0.0]]
        return geom_base.from_data(symbs, xyzs, angstrom=True)

    if len(symb_dct) == 2:
        bkey = frozenset(symb_dct.keys())
        symbs = list(symb_dct.values())
        key1, key2 = bkey
        bdist = ts.heuristic_bond_distance(gra, key1, key2, angstrom=True)
        xyzs = [[0.0, 0.0, 0.0], [bdist, 0.0, 0.0]]
        return geom_base.from_data(symbs, xyzs, angstrom=True)

    return None
