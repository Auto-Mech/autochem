""" objects that everything else depends on
"""
import itertools
import numpy
from phydat import phycon
from phydat import ptab
from automol import error
from automol import embed
from automol.graph._graph_dep import atoms
from automol.graph._graph_dep import bonds
from automol.graph._graph_dep import atom_keys
from automol.graph._graph_dep import bond_keys
from automol.graph._graph_dep import atom_symbols
from automol.graph._graph_dep import atom_stereo_parities
from automol.graph._graph_dep import bond_stereo_parities
from automol.graph._graph_dep import set_atom_stereo_parities
from automol.graph._graph_dep import set_bond_stereo_parities
from automol.graph._graph_dep import string
from automol.graph._graph_dep import relabel
# getters
from automol.graph._graph_dep import implicit
from automol.graph._graph_dep import explicit
from automol.graph._graph_dep import atom_neighborhoods
from automol.graph._graph_dep import bond_neighborhoods
from automol.graph._graph_dep import atoms_neighbor_atom_keys
from automol.graph._graph_dep import subgraph
from automol.graph._graph_dep import atom_stereo_sorted_neighbor_atom_keys
from automol.graph._graph_dep import resonance_dominant_atom_hybridizations
from automol.graph._graph_dep import sp2_bond_keys
# stereo
from automol.graph._graph_dep import atom_stereo_keys
from automol.graph._graph_dep import bond_stereo_keys
from automol.util import dict_
import automol.create as _create
from automol.graph import _networkx
from automol.graph import _igraph
# geometry tools
from automol.graph.geom import coordinates
from automol.graph.geom import geometry_join


# embed
def fake_stereo_geometry(gra, ntries=5, max_dist_err=0.5):
    """ generate a fake stereo geometry
    """
    # determine stereo "groups" with geometrically interdependent chirality
    atm_ngbs_dct = atom_neighborhoods(gra)
    bnd_ngbs_dct = bond_neighborhoods(gra)
    atm_ste_keys = atom_stereo_keys(gra)
    bnd_ste_keys = bond_stereo_keys(gra)
    atm_ste_groups = list(
        map(atom_keys, map(atm_ngbs_dct.__getitem__, atm_ste_keys)))
    bnd_ste_groups = list(
        map(atom_keys, map(bnd_ngbs_dct.__getitem__, bnd_ste_keys)))

    ste_groups = _aggregate_connected_groups(atm_ste_groups + bnd_ste_groups)

    ste_groups = list(map(sorted, ste_groups))

    natms = 0
    geo_idx_dct = {}
    geo = ()
    for group in ste_groups:
        group_geo = geometry(
            gra, keys=group, ntries=ntries,
            max_dist_err=max_dist_err)
        group_natms = len(group)

        idxs = list(range(natms, natms+group_natms))
        geo_idx_dct.update(dict(zip(group, idxs)))

        natms += group_natms
        geo = geometry_join(geo, group_geo)

    return geo, geo_idx_dct


def _aggregate_connected_groups(keys_lst):
    """ from a list of groups, join the ones with common atoms
    """
    groups = []
    for keys in keys_lst:
        idxs = list(idx for idx, group in enumerate(groups)
                    if group & keys)
        if not idxs:
            groups.append(keys)
        else:
            groups[idxs[0]] |= keys
            for idx in idxs[1:]:
                groups[idxs[0]] |= groups[idx]
                groups.pop(idx)

    groups = tuple(groups)

    return groups


def geometry(gra, keys=None, ntries=5, max_dist_err=0.2):
    """ sample a qualitatively-correct stereo geometry

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
        "Use automol.graph.explicit() to convert to an explicit graph.")

    # 0. Get keys and symbols
    symb_dct = atom_symbols(gra)

    keys = sorted(atom_keys(gra)) if keys is None else keys
    symbs = tuple(map(symb_dct.__getitem__, keys))

    # 1. Generate bounds matrices
    lmat, umat = distance_bounds_matrices(gra, keys)
    chi_dct = chirality_constraint_bounds(gra, keys)
    pla_dct = planarity_constraint_bounds(gra, keys)
    conv1_ = qualitative_convergence_checker_(gra, keys)
    conv2_ = embed.distance_convergence_checker_(lmat, umat, max_dist_err)

    def conv_(xmat, err, grad):
        return conv1_(xmat, err, grad) & conv2_(xmat, err, grad)

    # 2. Generate coordinates with correct stereo, trying a few times
    for _ in range(ntries):
        xmat = embed.sample_raw_distance_coordinates(lmat, umat, dim4=True)
        xmat, conv = embed.cleaned_up_coordinates(
            xmat, lmat, umat, pla_dct=pla_dct, chi_dct=chi_dct, conv_=conv_)
        if conv:
            break

    if not conv:
        raise error.FailedGeometryGenerationError

    # 3. Generate a geometry data structure from the coordinates
    xyzs = xmat[:, :3]
    geo = _create.geom.from_data(symbs, xyzs, angstrom=True)

    return geo


#
def distance_bounds_matrices(gra, keys, sp_dct=None):
    """ initial distance bounds matrices

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
    for (idx1, key1), (idx2, key2) in itertools.combinations(
            enumerate(keys), 2):
        if key2 in sp_dct[key1]:
            path = sp_dct[key1][key2]
            ldist, udist = bounds_(path)
            lmat[idx1, idx2] = lmat[idx2, idx1] = ldist
            umat[idx1, idx2] = umat[idx2, idx1] = udist
        else:
            # they are disconnected
            lmat[idx1, idx2] = lmat[idx2, idx1] = closest_approach(
                gra, key1, key2)
            umat[idx1, idx2] = umat[idx2, idx1] = 999

        assert lmat[idx1, idx2] <= umat[idx1, idx2], (
            "Lower bound exceeds upper bound. This is a bug!\n"
            "{}\npath: {}\n"
            .format(string(gra, one_indexed=False), str(path)))

    return lmat, umat


def chirality_constraint_bounds(gra, keys):
    """ bounds for enforcing chirality restrictions
    """
    ste_keys = set(atom_stereo_keys(gra)) & set(keys)
    par_dct = atom_stereo_parities(gra)
    ngb_key_dct = atoms_neighbor_atom_keys(gra)

    def _chirality_constraint(key):
        ngb_keys = ngb_key_dct[key]
        ngb_keys = atom_stereo_sorted_neighbor_atom_keys(gra, key, ngb_keys)
        idxs = tuple(map(keys.index, ngb_keys))
        vol_range = (-999., -7.) if par_dct[key] else (+7., +999.)
        return idxs, vol_range

    chi_dct = dict(map(_chirality_constraint, ste_keys))
    return chi_dct


def planarity_constraint_bounds(gra, keys):
    """ bounds for enforcing planarity restrictions
    """
    ngb_key_dct = atoms_neighbor_atom_keys(gra)
    ngb_dct = bond_neighborhoods(gra)
    bnd_keys = [bnd_key for bnd_key in sp2_bond_keys(gra)
                if atom_keys(ngb_dct[bnd_key]) <= set(keys)]

    def _planarity_constraints(bnd_key):
        key1, key2 = sorted(bnd_key)
        key1ab = sorted(ngb_key_dct[key1] - {key2})
        key2ab = sorted(ngb_key_dct[key2] - {key1})

        lst = []

        # I don't think the order of the keys matters, but I tried to be
        # roughly consistent with Figure 8 in the Blaney Dixon paper
        if len(key1ab) == 2 and len(key2ab) == 2:
            lst.append(tuple(map(keys.index, key1ab + key2ab)))
        if len(key1ab) == 2:
            lst.append(tuple(map(keys.index, [key1, key2] + key1ab)))
        if len(key2ab) == 2:
            lst.append(tuple(map(keys.index, [key1, key2] + key2ab)))
        if (len(key1ab) == 2 and len(key2ab) == 1) or (
                len(key1ab) == 1 and len(key2ab) == 2):
            lst.append(tuple(map(keys.index, [key1] + key1ab + key2ab)))
            lst.append(tuple(map(keys.index, [key2] + key1ab + key2ab)))
        if len(key1ab) == 1 and len(key2ab) == 1:
            lst.append(tuple(map(keys.index, [key1, key2] + key1ab + key2ab)))

        return tuple(lst)

    const_dct = {
        idxs: (-0.5, +0.5) for idxs in
        itertools.chain(*map(_planarity_constraints, bnd_keys))}

    return const_dct


def qualitative_convergence_checker_(gra, keys, rqq_bond_max=1.8,
                                     rqh_bond_max=1.3, rhh_bond_max=1.1,
                                     bond_nobond_diff=0.3):
    """ a convergence checker for error minimization, checking that the
    geometry is qualitatively correct (correct connectivity and stereo)
    """
    symb_dct = atom_symbols(gra)
    pairs = set(map(frozenset, itertools.combinations(keys, 2)))

    bnd_keys = pairs & bond_keys(gra)
    nob_keys = pairs - bond_keys(gra)

    nob_symbs = tuple(tuple(map(symb_dct.__getitem__, nob_key))
                      for nob_key in nob_keys)
    bnd_symbs = tuple(tuple(map(symb_dct.__getitem__, bnd_key))
                      for bnd_key in bnd_keys)
    nob_idxs = tuple(tuple(map(keys.index, nob_key)) for nob_key in nob_keys)
    bnd_idxs = tuple(tuple(map(keys.index, bnd_key)) for bnd_key in bnd_keys)

    bnd_udists = tuple((rqq_bond_max if 'H' not in symb else
                        rhh_bond_max if set(symb) == {'H'} else
                        rqh_bond_max) for symb in bnd_symbs)

    diff = bond_nobond_diff
    nob_ldists = tuple((rqq_bond_max+diff if 'H' not in symb else
                        rhh_bond_max+diff if set(symb) == {'H'} else
                        rqh_bond_max+diff) for symb in nob_symbs)

    bnd_idxs += tuple(map(tuple, map(reversed, bnd_idxs)))
    bnd_idx_vecs = tuple(map(list, zip(*bnd_idxs)))
    bnd_udists *= 2

    nob_idxs += tuple(map(tuple, map(reversed, nob_idxs)))
    nob_idx_vecs = tuple(map(list, zip(*nob_idxs)))
    nob_ldists *= 2

    symbs = tuple(map(symb_dct.__getitem__, keys))
    geo_idx_dct = dict(map(reversed, enumerate(keys)))
    atm_ste_keys = atom_stereo_keys(gra) & set(keys)
    bnd_ste_keys = bond_stereo_keys(gra) & bnd_keys
    atm_ste_par_dct = atom_stereo_parities(gra)
    bnd_ste_par_dct = bond_stereo_parities(gra)

    def _is_converged(xmat, err, grad):
        assert err and numpy.any(grad)
        xyzs = xmat[:, :3]
        dmat = embed.distance_matrix_from_coordinates(xyzs)

        # check for correct connectivity
        connectivity_check = (
            (numpy.all(dmat[bnd_idx_vecs] < bnd_udists)
             if bnd_udists else True) and
            (numpy.all(dmat[nob_idx_vecs] > nob_ldists)
             if nob_ldists else True))

        # check for correct stereo parities
        geo = _create.geom.from_data(symbs, xyzs, angstrom=True)
        atom_stereo_check = all(
            (_atom_stereo_parity_from_geometry(gra, atm_key, geo, geo_idx_dct)
             == atm_ste_par_dct[atm_key])
            for atm_key in atm_ste_keys)

        bond_stereo_check = all(
            (_bond_stereo_parity_from_geometry(gra, bnd_key, geo, geo_idx_dct)
             == bnd_ste_par_dct[bnd_key])
            for bnd_key in bnd_ste_keys)

        return connectivity_check and atom_stereo_check and bond_stereo_check

    return _is_converged


def atom_shortest_paths(gra):
    """ shortest paths between any two atoms in the graph

    :returns: a 2d dictionary keyed by pairs of atoms
    """
    nxg = _networkx.from_graph(gra)
    sp_dct = dict(_networkx.all_pairs_shortest_path(nxg))
    return sp_dct


def van_der_waals_radius(gra, key):
    """ van der waals radius for an atom in the graph (in angstroms)
    """

    symb_dct = atom_symbols(gra)
    symb = symb_dct[key]
    rad = ptab.van_der_waals_radius(symb) * phycon.BOHR2ANG

    return rad


def closest_approach(gra, key1, key2):
    """ closest approach between atoms, based on their van der Waals radii

    Warning: The scaling factor on the van der waals radii was arbitrarily
    chosen based on limited tests and may need to be lowered
    """
    vdw_scaling_factor = 0.75
    dist = (van_der_waals_radius(gra, key1) +
            van_der_waals_radius(gra, key2)) * vdw_scaling_factor
    return dist


def path_distance_bounds_(gra):
    """ upper distance bound between two ends of a path

    :param gra: molecular graph
    :param path: the shortest path between two atoms
    :type path: list or tuple
    """

    hyb_dct = resonance_dominant_atom_hybridizations(gra)
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
                    gra, *path, hyb_dct=hyb_dct)
            else:
                a123 = (rsz - 2.) * 180. / rsz
                la123 = a123 - 10.
                ua123 = a123 + 10.
                lrdist = heuristic_bond_angle_distance(gra, *path, a123=la123)
                urdist = heuristic_bond_angle_distance(gra, *path, a123=ua123)
                odist = heuristic_bond_angle_distance(
                    gra, *path, hyb_dct=hyb_dct)
                ldist = min(lrdist, odist)
                udist = max(urdist, odist)
        elif len(path) == 4:
            if rsz == 0:
                key2, key3 = path[1:3]
                bnd_key23 = frozenset({key2, key3})
                # handle bond stereo here
                if bnd_key23 in ste_bnd_keys:
                    key2_ngbs = atom_stereo_sorted_neighbor_atom_keys(
                        gra, key2, atm_ngb_keys[key2]-{key3})
                    key3_ngbs = atom_stereo_sorted_neighbor_atom_keys(
                        gra, key3, atm_ngb_keys[key3]-{key2})
                    pos2 = key2_ngbs.index(path[0])
                    pos3 = key3_ngbs.index(path[-1])
                    cis = bnd_par_dct[bnd_key23] != (pos2 != pos3)
                    dih = 0. if cis else 180.
                    ldist = udist = heuristic_torsion_angle_distance(
                        gra, *path, d1234=dih, degree=True, hyb_dct=hyb_dct)
                else:
                    ldist = heuristic_torsion_angle_distance(
                        gra, *path, d1234=0., degree=True, hyb_dct=hyb_dct)
                    udist = heuristic_torsion_angle_distance(
                        gra, *path, d1234=180., degree=True, hyb_dct=hyb_dct)
            else:
                ang = (rsz - 2.) * 180. / rsz
                rdist = heuristic_torsion_angle_distance(
                    gra, *path, d1234=0., degree=True, a123=ang, a234=ang)
                cdist = heuristic_torsion_angle_distance(
                    gra, *path, d1234=0., degree=True, hyb_dct=hyb_dct)
                tdist = heuristic_torsion_angle_distance(
                    gra, *path, d1234=180., degree=True, hyb_dct=hyb_dct)
                ldist = min(rdist, cdist)
                udist = max(rdist, tdist)
        # otherwise, just do the sum of the distances between atoms along the
        # path
        else:
            # we can't handle disconnected points, because in that case the
            # path is [] and there is no way to recover the keys
            assert len(path) > 2
            ldist = closest_approach(gra, path[0], path[-1])
            udist = 999.

        return ldist, udist

    return _distance_bounds


# bond distances
XY_DIST = 1.5       # angstroms
XH_DIST = 1.1       # angstroms

# bond angles
TET_ANG = 109.4712  # degrees
TRI_ANG = 120.      # degrees
LIN_ANG = 180.      # degrees


def heuristic_bond_distance(gra, key1, key2, angstrom=True, check=False):
    """ heuristic bond distance (in angstroms)
    """
    if check:
        assert key1 in atoms_neighbor_atom_keys(gra)[key2]

    symb_dct = atom_symbols(gra)
    symb1 = symb_dct[key1]
    symb2 = symb_dct[key2]

    if ptab.to_number(symb1) == 1 or ptab.to_number(symb2) == 1:
        dist = XH_DIST
    else:
        dist = XY_DIST

    dist *= 1 if angstrom else phycon.ANG2BOHR

    return dist


def heuristic_bond_angle(gra, key1, key2, key3, degree=False, check=False,
                         hyb_dct=None):
    """ heuristic bond angle

    If being reused multiple times, you can speed this up by passing in the
    hybridizations, so they don't need to be recalculated
    """
    if check:
        assert {key1, key3} <= set(atoms_neighbor_atom_keys(gra)[key2])

    if hyb_dct is None:
        hyb_dct = resonance_dominant_atom_hybridizations(gra)

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


def heuristic_bond_angle_distance(gra, key1, key2, key3,
                                  d12=None, d23=None, angstrom=True,
                                  a123=None, degree=True, hyb_dct=None):
    """ heuristic distance between atoms at two ends of a bond angle

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
            gra, key1, key2, key3, degree=False, hyb_dct=hyb_dct)
    else:
        a123 = a123 * phycon.DEG2RAD if degree else a123

    d12 = heuristic_bond_distance(gra, key1, key2, angstrom=angstrom)
    d23 = heuristic_bond_distance(gra, key2, key3, angstrom=angstrom)
    d13 = numpy.sqrt(d12**2 + d23**2 - 2*d12*d23*numpy.cos(a123))
    return d13


def heuristic_torsion_angle_distance(gra, key1, key2, key3, key4,
                                     angstrom=True,
                                     a123=None, a234=None, d1234=None,
                                     degree=True, hyb_dct=None):
    """ heuristic max distance between atoms at two ends of a torsion angle

    Formula rearranged from eq. 2.15 in the following paper. Note that in the
    denominator of this formula there is a typo: d02 and d12 should be switched

    Havel, T. F.; Kunz, I.D.; Crippen, G. M.; "The Theory and Practice of
    Distance Geometry"; Bulletin of Mathematical Biology; Vol. 45, No. 5.,
    pp.  665-720, 1983.
    """
    if d1234 is None:
        d1234 = numpy.pi
    else:
        d1234 *= phycon.DEG2RAD if degree else 1.

    d12 = heuristic_bond_distance(gra, key1, key2, angstrom=angstrom)
    d23 = heuristic_bond_distance(gra, key2, key3, angstrom=angstrom)
    d34 = heuristic_bond_distance(gra, key3, key4, angstrom=angstrom)
    d13 = heuristic_bond_angle_distance(
        gra, key1, key2, key3, angstrom=angstrom, a123=a123, hyb_dct=hyb_dct)
    d24 = heuristic_bond_angle_distance(
        gra, key2, key3, key4, angstrom=angstrom, a123=a234, hyb_dct=hyb_dct)

    term1 = (d12**2 + d23**2 - d13**2)*(d23**2 + d34**2 - d24**2)
    term2 = 2*d23**2 * (d13**2 + d24**2 - d23**2)
    denom = numpy.sqrt((4*d12**2 * d23**2 - (d12**2 + d23**2 - d13**2)**2) *
                       (4*d23**2 * d34**2 - (d23**2 + d34**2 - d24**2)**2))

    d14 = numpy.sqrt((term1 + term2 - numpy.cos(d1234) * denom) /
                     (2 * d23**2))
    return d14


def shared_ring_size(keys, rng_keys_lst):
    """ determine whether these keys share a ring and, if so, determine the
    size
    """
    rng_keys = next((rng_keys for rng_keys in rng_keys_lst
                     if set(keys) <= set(rng_keys)), ())
    natms = len(rng_keys)
    return natms


def rings_atom_keys(gra):
    """ atom keys for each ring in the graph sorted by connectivity (minimal basis)
    """
    rng_atm_keys_lst = frozenset(
        map(sorted_ring_atom_keys_from_bond_keys, rings_bond_keys(gra)))
    return rng_atm_keys_lst


def sorted_ring_atom_keys(rng):
    """ get a ring's atom keys, sorted in order of connectivity
    """
    return sorted_ring_atom_keys_from_bond_keys(bond_keys(rng))


def sorted_ring_atom_keys_from_bond_keys(rng_bnd_keys):
    """ get a ring's atom keys, sorted in order of connectivity, from its bond
    keys
    """
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


def _set_atom_stereo_from_geometry(gra, atm_keys, geo, geo_idx_dct):
    assert gra == explicit(gra)

    atm_pars = [
        _atom_stereo_parity_from_geometry(gra, atm_key, geo, geo_idx_dct)
        for atm_key in atm_keys]
    gra = set_atom_stereo_parities(gra, dict(zip(atm_keys, atm_pars)))
    return gra


def _set_bond_stereo_from_geometry(gra, bnd_keys, geo, geo_idx_dct):
    assert gra == explicit(gra)

    bnd_pars = [
        _bond_stereo_parity_from_geometry(gra, bnd_key, geo, geo_idx_dct)
        for bnd_key in bnd_keys]
    gra = set_bond_stereo_parities(gra, dict(zip(bnd_keys, bnd_pars)))
    return gra


def _bond_stereo_parity_from_geometry(gra, bnd_key, geo, geo_idx_dct):
    """ get the current stereo parity of a bond from its geometry
    """
    atm1_key, atm2_key = bnd_key
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    atm1_ngb_keys = atm_ngb_keys_dct[atm1_key] - {atm2_key}
    atm2_ngb_keys = atm_ngb_keys_dct[atm2_key] - {atm1_key}

    atm1_ngb_keys = atom_stereo_sorted_neighbor_atom_keys(
        gra, atm1_key, atm1_ngb_keys)
    atm2_ngb_keys = atom_stereo_sorted_neighbor_atom_keys(
        gra, atm2_key, atm2_ngb_keys)

    # get the top priority neighbor keys on each side
    atm1_ngb_key = atm1_ngb_keys[0]
    atm2_ngb_key = atm2_ngb_keys[0]

    # determine the parity based on the coordinates
    xyzs = coordinates(geo)
    atm1_xyz = xyzs[geo_idx_dct[atm1_key]]
    atm2_xyz = xyzs[geo_idx_dct[atm2_key]]
    atm1_ngb_xyz = xyzs[geo_idx_dct[atm1_ngb_key]]
    atm2_ngb_xyz = xyzs[geo_idx_dct[atm2_ngb_key]]
    atm1_bnd_vec = numpy.subtract(atm1_ngb_xyz, atm1_xyz)
    atm2_bnd_vec = numpy.subtract(atm2_ngb_xyz, atm2_xyz)
    dot_val = numpy.vdot(atm1_bnd_vec, atm2_bnd_vec)
    assert dot_val != 0.  # for now, assume no collinear
    par = dot_val > 0.
    return par


def connected_components(gra):
    """ connected components in the graph
    """
    cmp_gra_atm_keys_lst = connected_components_atom_keys(gra)
    cmp_gras = tuple(subgraph(gra, cmp_gra_atm_keys, stereo=True)
                     for cmp_gra_atm_keys in cmp_gra_atm_keys_lst)
    return cmp_gras


def connected_components_atom_keys(gra):
    """ atom keys for each connected component in the graph
    """
    nxg = _networkx.from_graph(gra)
    cmp_gra_atm_keys_lst = _networkx.connected_component_atom_keys(nxg)
    return cmp_gra_atm_keys_lst


def union(gra1, gra2, check=True):
    """ a union of two graphs
    """
    if check:
        assert not atom_keys(gra1) & atom_keys(gra2)
    atm_dct = {}
    atm_dct.update(atoms(gra1))
    atm_dct.update(atoms(gra2))

    bnd_dct = {}
    bnd_dct.update(bonds(gra1))
    bnd_dct.update(bonds(gra2))
    return _create.graph.from_atoms_and_bonds(atm_dct, bnd_dct)


def transform_keys(gra, atm_key_func):
    """ transform atom keys with a function
    """
    atm_keys = atom_keys(gra)
    atm_key_dct = dict(zip(atm_keys, map(atm_key_func, atm_keys)))
    return relabel(gra, atm_key_dct)


def _atom_stereo_parity_from_geometry(gra, atm_key, geo, geo_idx_dct):
    """ get the current stereo parity of an atom from its geometry
    """
    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    # sort the neighbor keys by stereo priority
    atm_ngb_keys = atom_stereo_sorted_neighbor_atom_keys(
        gra, atm_key, atm_ngb_keys)

    # determine the parity based on the coordinates
    xyzs = coordinates(geo)
    atm_ngb_idxs = dict_.values_by_key(geo_idx_dct, atm_ngb_keys)
    atm_ngb_xyzs = [xyzs[idx] for idx in atm_ngb_idxs]
    det_mat = numpy.ones((4, 4))
    det_mat[:, :3] = atm_ngb_xyzs
    det_val = numpy.linalg.det(det_mat)
    assert det_val != 0.  # for now, assume no four-atom planes
    par = det_val > 0.
    return par


def backbone_isomorphic(gra1, gra2):
    """ are these molecular graphs backbone isomorphic?
    """
    return backbone_isomorphism(gra1, gra2) is not None


# def backbone_isomorphism(gra1, gra2, igraph=True):
def backbone_isomorphism(gra1, gra2, igraph=False):
    """ graph backbone isomorphism

    for implicit graphs, this is the relabeling of `gra1` to produce `gra2`
    for other graphs, it gives the correspondences between backbone atoms
    """
    gra1 = implicit(gra1)
    gra2 = implicit(gra2)
    if igraph:
        igr1 = _igraph.from_graph(gra1)
        igr2 = _igraph.from_graph(gra2)
        iso_dcts = _igraph.isomorphisms(igr1, igr2)
        iso_dct = iso_dcts[0] if iso_dcts else None
    else:
        nxg1 = _networkx.from_graph(gra1)
        nxg2 = _networkx.from_graph(gra2)
        iso_dct = _networkx.isomorphism(nxg1, nxg2)
    return iso_dct
