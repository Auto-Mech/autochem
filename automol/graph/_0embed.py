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
Step 7 is performed byt the function automol.embed.cleaned_up_coordinates()
"""
import itertools

import numpy
from automol import embed, error, geom, zmat
from automol.geom import base as geom_base
from automol.graph.base import (
    add_bonds,
    atom_hybridizations,
    atom_keys,
    atom_shortest_paths,
    atom_stereo_keys,
    atom_stereo_parities,
    atom_symbols,
    atoms_neighbor_atom_keys,
    bond_keys,
    bond_neighborhoods,
    bond_orders,
    bond_stereo_keys,
    bond_stereo_parities,
    explicit,
    geometry_atom_parity,
    geometry_bond_parity,
    heuristic_bond_angle,
    heuristic_bond_distance,
    heuristic_bond_distance_limit,
    is_ts_graph,
    rigid_planar_bond_keys,
    rings_atom_keys,
    rotational_bond_keys,
    set_stereo_from_geometry,
    stereo_keys,
    stereo_parities,
    stereocenter_candidates,
    string,
    subgraph,
    to_local_stereo,
    ts,
    without_bonds_by_orders,
    without_dummy_atoms,
    without_stereo,
)
from automol.util import dict_, heuristic
from phydat import phycon


# # geometry embedding functions
def embed_geometry(gra, keys=None, ntries=5, max_dist_err=0.2):
    """sample a qualitatively-correct stereo geometry

    :param gra: the graph, which may or may not have stereo
    :param keys: graph keys, in the order in which they should appear in the
        geometry
    :param ntries: number of tries for finding a valid geometry
    :param max_dist_err: maximum distance error convergence threshold

    Qualitatively-correct means it has the right connectivity and the right
    stero parities, but its bond lengths and bond angles may not be
    quantitatively realistic
    """
    assert gra == explicit(gra), (
        "Graph => geometry conversion requires explicit hydrogens!\n"
        "Use automol.graph.explicit() to convert to an explicit graph."
    )

    symbs = atom_symbols(gra)
    if len(symbs) == 1:
        symb = list(symbs.values())[0]
        geo = ((symb, (0.00, 0.00, 0.00)),)
        return geo

    # For simplicity, convert to local stereo and use local stereo
    # throughout
    loc_gra = to_local_stereo(gra)

    # 0. Get keys and symbols
    symb_dct = atom_symbols(loc_gra)

    keys = sorted(atom_keys(loc_gra)) if keys is None else keys
    symbs = tuple(map(symb_dct.__getitem__, keys))

    # 1. Generate bounds matrices
    lmat, umat = distance_bounds_matrices(loc_gra, keys)
    chi_dct = chirality_constraint_bounds(loc_gra, keys)
    pla_dct = planarity_constraint_bounds(loc_gra, keys)
    conv1_ = qualitative_convergence_checker_(loc_gra, keys)
    conv2_ = embed.distance_convergence_checker_(lmat, umat, max_dist_err)

    def conv_(xmat, err, grad):
        return conv1_(xmat, err, grad) & conv2_(xmat, err, grad)

    # 2. Generate coordinates with correct stereo, trying a few times
    for _ in range(ntries):
        xmat = embed.sample_raw_distance_coordinates(lmat, umat, dim4=True)
        xmat, conv = embed.cleaned_up_coordinates(
            xmat, lmat, umat, pla_dct=pla_dct, chi_dct=chi_dct, conv_=conv_
        )
        if conv:
            break

    if not conv:
        raise error.FailedGeometryGenerationError(f"Bad gra {string(loc_gra)}")

    # 3. Generate a geometry data structure from the coordinates
    xyzs = xmat[:, :3]
    geo = geom_base.from_data(symbs, xyzs, angstrom=True)

    return geo


def clean_geometry(
    gra,
    geo,
    rct_geos=None,
    dist_range_dct=None,
    max_dist_err=0.2,
    stereo=True,
    local_stereo=False,
    relax_angles=False,
    log=False,
    none_if_failed=True,
):
    """Clean up a geometry based on this graph, removing any bonds that
    aren't supposed to be there

    :param gra: molecular graph with stereo parities
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param rct_geos: Use reactant geometries to determine distance bounds for a TS
    :type rct_geos: List[automol geometry data structure], optional
    :param dist_range_dct: Override a subsets of distance ranges, defaults to None
    :type dist_range_dct: dict, optional
    :param stereo: Take stereochemistry into consideration? defaults to True
    :type stereo: bool, optional
    :param local_stereo: Does the graph have local stereo assignments? defaults to False
    :type local_stereo: bool, optional
    :param relax_angles: Allow angles to relax?, defaults to False
    :type relax_angles: bool, optional
    :param none_if_failed: Return `None` if the geometry doesn't match? defaults to True
    type none_if_failed: bool, optional
    """
    gra = gra if local_stereo else to_local_stereo(gra)

    symb_dct = atom_symbols(gra)

    # Build monatomics and diatomics directly
    if len(symb_dct) == 1:
        symbs = list(symb_dct.values())
        xyzs = [[0.0, 0.0, 0.0]]
        return geom_base.from_data(symbs, xyzs, angstrom=True)

    if len(symb_dct) == 2:
        bkey = frozenset(symb_dct.keys())
        symbs = list(symb_dct.values())
        if dist_range_dct and bkey in dist_range_dct:
            bdist = sum(dist_range_dct[bkey]) / 2.0
        else:
            key1, key2 = bkey
            bdist = heuristic_bond_distance(gra, key1, key2, angstrom=True)
        xyzs = [[0.0, 0.0, 0.0], [bdist, 0.0, 0.0]]
        return geom_base.from_data(symbs, xyzs, angstrom=True)

    rct_geos = [geo] if rct_geos is None else rct_geos

    dist_gra = ts.reactants_graph(gra) if is_ts_graph(gra) else gra

    keys = sorted(atom_keys(gra))
    xmat = geom_base.coordinates(geo, angstrom=True)
    lmat, umat = distance_bounds_matrices(
        dist_gra,
        keys,
        rct_geos=rct_geos,
        dist_range_dct=dist_range_dct,
        relax_angles=relax_angles,
        relax_torsions=True,
    )

    pla_dct = planarity_constraint_bounds(gra, keys)
    if stereo:
        chi_dct = chirality_constraint_bounds(gra, keys)
    else:
        chi_dct = {}

    xmat, conv = embed.cleaned_up_coordinates(
        xmat,
        lmat,
        umat,
        chi_dct=chi_dct,
        pla_dct=pla_dct,
        max_dist_err=max_dist_err,
        log=log,
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
    cgra = without_stereo(gra)
    cgra0 = geom.graph_without_stereo(geo)

    cgra = without_dummy_atoms(cgra)
    cgra0 = without_dummy_atoms(cgra0)

    cgra_ = cgra0
    if is_ts_graph(cgra):
        r_bkeys = ts.reacting_bond_keys(cgra)
        r_ord_dct = dict_.by_key(bond_orders(cgra), r_bkeys)
        cgra_ = add_bonds(cgra_, keys=r_bkeys, ord_dct=r_ord_dct, check=False)

    # 1. Check that the connectivities match
    matches = cgra == cgra_
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

        ste_keys = stereo_keys(gra)
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


# # convergence checking
def qualitative_convergence_checker_(
    loc_gra,
    keys,
    dist_factor=None,
    bond_nobond_diff=0.3,
):
    """a convergence checker for error minimization, checking that the
    geometry is qualitatively correct (correct connectivity and stereo)
    """
    symb_dct = atom_symbols(loc_gra)
    pairs = set(map(frozenset, itertools.combinations(keys, 2)))

    bnd_keys = pairs & bond_keys(loc_gra)
    nob_keys = pairs - bond_keys(loc_gra)

    nob_symbs = tuple(tuple(map(symb_dct.__getitem__, nob_key)) for nob_key in nob_keys)
    bnd_symbs = tuple(tuple(map(symb_dct.__getitem__, bnd_key)) for bnd_key in bnd_keys)
    nob_idxs = tuple(tuple(map(keys.index, nob_key)) for nob_key in nob_keys)
    bnd_idxs = tuple(tuple(map(keys.index, bnd_key)) for bnd_key in bnd_keys)

    bnd_udists = tuple(
        heuristic.bond_distance_limit(s1, s2, dist_factor=dist_factor, angstrom=True)
        for s1, s2 in bnd_symbs
    )

    diff = bond_nobond_diff
    nob_ldists = tuple(
        diff
        + heuristic.bond_distance_limit(s1, s2, dist_factor=dist_factor, angstrom=True)
        for s1, s2 in nob_symbs
    )

    bnd_idxs += tuple(map(tuple, map(reversed, bnd_idxs)))
    bnd_idx_vecs = tuple(map(list, zip(*bnd_idxs)))
    bnd_udists *= 2

    nob_idxs += tuple(map(tuple, map(reversed, nob_idxs)))
    nob_idx_vecs = tuple(map(list, zip(*nob_idxs)))
    nob_ldists *= 2

    symbs = tuple(map(symb_dct.__getitem__, keys))
    geo_idx_dct = dict(map(reversed, enumerate(keys)))
    atm_ste_keys = atom_stereo_keys(loc_gra) & set(keys)
    bnd_ste_keys = bond_stereo_keys(loc_gra) & bnd_keys
    atm_ste_par_dct = atom_stereo_parities(loc_gra)
    bnd_ste_par_dct = bond_stereo_parities(loc_gra)

    def _is_converged(xmat, err, grad):
        assert err or not err
        assert grad or not grad
        xyzs = xmat[:, :3]
        dmat = embed.distance_matrix_from_coordinates(xyzs)

        # check for correct connectivity
        connectivity_check = (
            numpy.all(dmat[bnd_idx_vecs] < bnd_udists) if bnd_udists else True
        ) and (numpy.all(dmat[nob_idx_vecs] > nob_ldists) if nob_ldists else True)

        # check for correct stereo parities
        geo = geom_base.from_data(symbs, xyzs, angstrom=True)
        atom_stereo_check = all(
            (
                geometry_atom_parity(loc_gra, geo, k, geo_idx_dct=geo_idx_dct)
                == atm_ste_par_dct[k]
            )
            for k in atm_ste_keys
        )

        bond_stereo_check = all(
            (
                geometry_bond_parity(loc_gra, geo, k, geo_idx_dct=geo_idx_dct)
                == bnd_ste_par_dct[k]
            )
            for k in bnd_ste_keys
        )

        return connectivity_check and atom_stereo_check and bond_stereo_check

    return _is_converged


# # bounds matrices
def distance_bounds_matrices(
    gra,
    keys,
    dist_range_dct=(),
    rct_geos=None,
    relax_angles=False,
    relax_torsions=False,
    sp_dct=None,
    angstrom=True,
):
    """generates initial distance bounds matrices for various different
    scenarios, allowing the geometry to be manipulated in different ways

    :param gra: molecular graph:
    :param keys: atom keys specifying the order of indices in the matrix
    :param dist_range_dct: distance ranges for specific atoms in the graph
    :param rct_geos: Use reactant geometries to determine distance bounds for a TS
    :type rct_geos: List[automol geometry data structure], optional
    :param relax_angles: whether or not to allow angles to change from
        their value in the reactant geometries
    :param relax_torsions: whether or not to allow torsions to change from
        their value in the reactant geometries
    :param sp_dct: a 2d dictionary giving the shortest path between any pair of
        atoms in the graph
    """
    sp_dct = atom_shortest_paths(gra) if sp_dct is None else sp_dct

    natms = len(keys)

    lmat, umat = _distance_bounds_matrices(gra, keys, sp_dct=sp_dct)

    # save the current values so that we can overwrite the fixed torsions below
    lmat_old = numpy.copy(lmat)
    umat_old = numpy.copy(umat)

    # 1. set known geometric parameters
    if rct_geos:
        xmats = [geom_base.coordinates(geo, angstrom=angstrom) for geo in rct_geos]
        dmats = list(map(embed.distance_matrix_from_coordinates, xmats))

        start = 0
        for dmat in dmats:
            dim, _ = numpy.shape(dmat)
            end = start + dim

            lmat[start:end, start:end] = dmat
            umat[start:end, start:end] = dmat

            start = end

    # 2. reopen bounds on the torsions from the reactant
    if relax_torsions:
        rot_bnd_keys = rotational_bond_keys(gra)

        tors_ijs = [
            [i, j]
            for i, j in itertools.combinations(range(natms), 2)
            if j in sp_dct[i]
            and len(sp_dct[i][j]) >= 4
            and frozenset(sp_dct[i][j][1:3]) in rot_bnd_keys
        ]

        tors_ijs += list(map(list, map(reversed, tors_ijs)))

        tors_idxs = tuple(map(list, zip(*tors_ijs)))

        if tors_idxs:
            lmat[tors_idxs] = lmat_old[tors_idxs]
            umat[tors_idxs] = umat_old[tors_idxs]

    # 3. reopen bounds on the angles from the reactant
    if relax_angles:
        ang_ijs = [
            [i, j]
            for i, j in itertools.combinations(range(natms), 2)
            if j in sp_dct[i] and len(sp_dct[i][j]) >= 3
        ]
        ang_ijs += list(map(list, map(reversed, ang_ijs)))

        ang_idxs = tuple(map(list, zip(*ang_ijs)))

        if ang_idxs:
            lmat[ang_idxs] = lmat_old[ang_idxs]
            umat[ang_idxs] = umat_old[ang_idxs]

    # 4. set distance bounds for the forming bonds
    if dist_range_dct:
        for bnd, (ldist, udist) in dist_range_dct.items():
            idx1 = tuple(bnd)
            idx2 = tuple(reversed(idx1))
            lmat[idx1] = lmat[idx2] = ldist
            umat[idx1] = umat[idx2] = udist

    return lmat, umat


def _distance_bounds_matrices(gra, keys, sp_dct=None):
    """initial distance bounds matrices

    :param gra: molecular graph
    :param keys: atom keys specifying the order of indices in the matrix
    :param sp_dct: a 2d dictionary giving the shortest path between any pair of
        atoms in the graph
    """
    assert set(keys) <= set(atom_keys(gra))

    sub_gra = subgraph(gra, keys, stereo=True)
    sp_dct = atom_shortest_paths(sub_gra) if sp_dct is None else sp_dct

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
def chirality_constraint_bounds(loc_gra, keys):
    """bounds for enforcing chirality restrictions"""
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


def planarity_constraint_bounds(gra, keys):
    """bounds for enforcing planarity restrictions"""
    nkeys_dct = atoms_neighbor_atom_keys(gra)
    # rp_dct = rigid_planar_bond_keys(gra, min_ncount=0, min_ring_size=0)
    ngb_dct = bond_neighborhoods(gra)
    bkeys = [
        bkey
        for bkey in rigid_planar_bond_keys(gra, min_ncount=0, min_ring_size=0)
        if atom_keys(ngb_dct[bkey]) <= set(keys)
    ]

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


def heuristic_torsion_angle_distance(
    gra,
    key1,
    key2,
    key3,
    key4,
    angstrom=True,
    a123=None,
    a234=None,
    d1234=None,
    degree=True,
    hyb_dct=None,
):
    """heuristic max distance between atoms at two ends of a torsion angle

    Formula rearranged from eq. 2.15 in the following paper. Note that in the
    denominator of this formula there is a typo: d02 and d12 should be switched

    Havel, T. F.; Kunz, I.D.; Crippen, G. M.; "The Theory and Practice of
    Distance Geometry"; Bulletin of Mathematical Biology; Vol. 45, No. 5.,
    pp.  665-720, 1983.
    """
    if d1234 is None:
        d1234 = numpy.pi
    else:
        d1234 *= phycon.DEG2RAD if degree else 1.0

    d12 = heuristic_bond_distance(gra, key1, key2, angstrom=angstrom)
    d23 = heuristic_bond_distance(gra, key2, key3, angstrom=angstrom)
    d34 = heuristic_bond_distance(gra, key3, key4, angstrom=angstrom)
    d13 = heuristic_bond_angle_distance(
        gra, key1, key2, key3, angstrom=angstrom, a123=a123, hyb_dct=hyb_dct
    )
    d24 = heuristic_bond_angle_distance(
        gra, key2, key3, key4, angstrom=angstrom, a123=a234, hyb_dct=hyb_dct
    )

    term1 = (d12**2 + d23**2 - d13**2) * (d23**2 + d34**2 - d24**2)
    term2 = 2 * d23**2 * (d13**2 + d24**2 - d23**2)
    term3 = (4 * d12**2 * d23**2 - (d12**2 + d23**2 - d13**2) ** 2) * (
        4 * d23**2 * d34**2 - (d23**2 + d34**2 - d24**2) ** 2
    )
    if term3 < 0.0:
        assert numpy.allclose(term3, 0)
        term3 = 0.0

    denom = numpy.sqrt(term3)

    term4 = (term1 + term2 - numpy.cos(d1234) * denom) / (2 * d23**2)
    if term4 < 0.0:
        assert numpy.allclose(term4, 0)
        term4 = 0.0

    d14 = numpy.sqrt(term4)
    return d14


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
    atm_ngb_keys = atoms_neighbor_atom_keys(gra)

    ste_bnd_keys = bond_stereo_keys(gra)
    bnd_par_dct = bond_stereo_parities(gra)

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
        elif len(path) == 4:
            if rsz == 0:
                key2, key3 = path[1:3]
                bnd_key23 = frozenset({key2, key3})
                # handle bond stereo here
                if bnd_key23 in ste_bnd_keys:
                    key2_ngbs = sorted(atm_ngb_keys[key2] - {key3})
                    key3_ngbs = sorted(atm_ngb_keys[key3] - {key2})
                    pos2 = key2_ngbs.index(path[0])
                    pos3 = key3_ngbs.index(path[-1])
                    cis = bnd_par_dct[bnd_key23] != (pos2 != pos3)
                    dih = 0.0 if cis else 180.0
                    ldist = udist = heuristic_torsion_angle_distance(
                        gra, *path, d1234=dih, degree=True, hyb_dct=hyb_dct
                    )
                else:
                    ldist = heuristic_torsion_angle_distance(
                        gra, *path, d1234=0.0, degree=True, hyb_dct=hyb_dct
                    )
                    udist = heuristic_torsion_angle_distance(
                        gra, *path, d1234=180.0, degree=True, hyb_dct=hyb_dct
                    )
            else:
                ang = (rsz - 2.0) * 180.0 / rsz
                rdist = heuristic_torsion_angle_distance(
                    gra, *path, d1234=0.0, degree=True, a123=ang, a234=ang
                )
                cdist = heuristic_torsion_angle_distance(
                    gra, *path, d1234=0.0, degree=True, hyb_dct=hyb_dct
                )
                tdist = heuristic_torsion_angle_distance(
                    gra, *path, d1234=180.0, degree=True, hyb_dct=hyb_dct
                )
                ldist = min(rdist, cdist)
                udist = max(rdist, tdist)
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
