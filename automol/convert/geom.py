""" geometry conversions
"""

import itertools
import functools
import numpy
from phydat import phycon, ptab
from automol import util
from automol import error
from automol.convert import _pyx2z
from automol.convert import _util
import automol.graph.vmat
from automol.convert.graph import inchi_with_sort_from_geometry
from automol.graph._graph_dep import atom_keys
from automol.graph._graph_dep import explicit
from automol.graph._graph_dep import atoms_neighbor_atom_keys
from automol.graph._graph_dep import atoms_sorted_neighbor_atom_keys
from automol.graph._graph_dep import dummy_atoms_neighbor_atom_key
from automol.graph._graph_dep import (
    without_dummy_atoms as _without_dummy_atoms_graph)
from automol.graph._graph_dep import implicit
from automol.graph._embed_dep import transform_keys
from automol.graph._embed_dep import union
from automol.graph._embed_dep import connected_components
from automol.graph._embed_dep import backbone_isomorphic
from automol.graph._stereo_geom import set_stereo_from_geometry
from automol.convert.inchi import object_to_hardcoded_inchi_by_key
from automol.convert.inchi import object_from_hardcoded_inchi_by_key
from automol.convert.inchi import equivalent
from automol.convert.inchi import has_stereo
from automol.convert.inchi import standard_form
from automol.convert.inchi import same_connectivity
from automol.convert.inchi import split as _split
from automol import create
from automol.convert import _rdkit
from automol.convert import _pybel
from automol.graph.geom import coordinates
from automol.graph.geom import count
from automol.graph.geom import is_atom
from automol.graph.geom import geometry_join
from automol.graph.geom import translate
from automol.graph.geom import symbols


def external_symmetry_factor(geo):
    """ Obtain the external symmetry factor for a geometry using x2z interface
        which determines the initial symmetry factor and then divides by the
        enantiomeric factor.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :rtype: float
    """

    if is_atom(geo):
        ext_sym_fac = 1.
    else:
        oriented_geom = _pyx2z.to_oriented_geometry(geo)
        ext_sym_fac = oriented_geom.sym_num()
        if oriented_geom.is_enantiomer():
            ext_sym_fac *= 0.5

    return ext_sym_fac


# geometry => graph
def connectivity_graph(geo,
                       rqq_bond_max=3.45, rqh_bond_max=2.6, rhh_bond_max=1.9):
    """ Generate a molecular graph from the molecular geometry that has information
        about bond connectivity.

        :param rqq_bond_max: maximum distance between heavy atoms
        :type rqq_bond_max: float
        :param rqh_bond_max: maximum distance between heavy atoms and hydrogens
        :type rqh_bond_max: float
        :param rhh_bond_max: maximum distance between hydrogens
        :type rhh_bond_max: float
        :rtype: automol molecular graph structure
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo)

    def _distance(idx_pair):
        xyz1, xyz2 = map(xyzs.__getitem__, idx_pair)
        dist = numpy.linalg.norm(numpy.subtract(xyz1, xyz2))
        return dist

    def _are_bonded(idx_pair):
        sym1, sym2 = map(symbs.__getitem__, idx_pair)
        dist = _distance(idx_pair)
        return (False if 'X' in (sym1, sym2) else
                (dist < rqh_bond_max) if 'H' in (sym1, sym2) else
                (dist < rhh_bond_max) if (sym1 == 'H' and sym2 == 'H') else
                (dist < rqq_bond_max))

    idxs = range(len(xyzs))
    atm_symb_dct = dict(enumerate(symbs))
    bnd_keys = tuple(
        map(frozenset, filter(_are_bonded, itertools.combinations(idxs, r=2))))

    bnd_ord_dct = {bnd_key: 1 for bnd_key in bnd_keys}

    gra = create.graph.from_data(atom_symbols=atm_symb_dct, bond_keys=bnd_keys,
                                 bond_orders=bnd_ord_dct)
    return gra


def graph(geo, stereo=True):
    """ Generate a molecular graph from the molecular geometry that has information
        about bond connectivity and if requested, stereochemistry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph data structure
    """
    gra = connectivity_graph(geo)
    if stereo:
        gra = set_stereo_from_geometry(gra, geo)

    return gra


# geometry => inchi
def inchi(geo, stereo=True):
    """ Generate an InChI string from a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: str
    """
    ich, _ = inchi_with_sort(geo, stereo=stereo)

    return ich


def inchi_with_sort(geo, stereo=True):
    """ Generate an InChI string from a molecular geometry. (Sort?)

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: str
    """
    ich = object_to_hardcoded_inchi_by_key(
        'geom', geo, comp=_compare)
    nums = None
    if ich is None:
        gra = connectivity_graph(geo)
        if not stereo:
            geo = None
            geo_idx_dct = None
        else:
            geo_idx_dct = dict(enumerate(range(count(geo))))
        ich, nums = inchi_with_sort_from_geometry(
            gra=gra, geo=geo, geo_idx_dct=geo_idx_dct)

    return ich, nums


def _compare(geo1, geo2):
    """ Check if the backbone atoms of two molecular geometries are similar.

        :param geo1: molecular geometry 1
        :type geo1: automol geometry data structure
        :param geo2: molecular geometry 2
        :type geo2: automol geometry data structure
        :rtype: bool
    """
    gra1 = _without_dummy_atoms_graph(connectivity_graph(geo1))
    gra2 = _without_dummy_atoms_graph(connectivity_graph(geo2))

    return backbone_isomorphic(gra1, gra2)


# geometry => formula
def formula(geo):
    """ Generate a stoichiometric formula dictionary from a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :type: dict[str: int]
    """

    symbs = symbols(geo)
    fml = _util.formula(symbs)

    return fml


def formula_string(geo):
    """ Generate a stoichiometric formula string from a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :type: dict[str: int]
    """

    symbs = symbols(geo)
    fml = _util.formula(symbs)
    fml_str = automol.formula.string(fml)

    return fml_str


def from_subset(geo, idxs):
    """ Generate a new molecular geometry from a subset of the atoms in an
        input geometry.

        (Rename this and put it under operations?)

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param idxs: indices representing the subset of atoms
        :type idxs: tuple(int)
        :rtype: automol moleculer geometry data structure
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo)

    symbs = list(map(symbs.__getitem__, idxs))
    xyzs = list(map(xyzs.__getitem__, idxs))

    return automol.create.geom.from_data(symbs, xyzs)


def insert_dummies_on_linear_atoms(geo, lin_idxs=None, gra=None, dist=1.,
                                   tol=5.):
    """ Insert dummy atoms over linear atoms in the geometry.

        :param geo: the geometry
        :type geo: automol molecular geometry data structure
        :param lin_idxs: the indices of the linear atoms; if None, indices are
            automatically determined from the geometry based on the graph
        :type lin_idxs: tuple(int)
        :param gra: the graph describing connectivity; if None, a connectivity
            graph will be generated using default distance thresholds
        :type gra: automol molecular graph data structure
        :param dist: distance of dummy atom from the linear atom, in angstroms
        :type dist: float
        :param tol: the tolerance threshold for linearity, in degrees
        :type tol: float
        :returns: geometry with dummy atoms inserted, along with a dictionary
            mapping the linear atoms onto their associated dummy atoms
        :rtype: automol molecular geometry data structure
    """

    lin_idxs = linear_atoms(geo) if lin_idxs is None else lin_idxs
    gra = connectivity_graph(geo) if gra is None else gra

    dummy_ngb_idxs = set(
        dummy_atoms_neighbor_atom_key(gra).values())
    assert not dummy_ngb_idxs & set(lin_idxs), (
        "Attempting to add dummy atoms on atoms that already have them: {}"
        .format(dummy_ngb_idxs & set(lin_idxs)))

    ngb_idxs_dct = atoms_sorted_neighbor_atom_keys(gra)

    xyzs = coordinates(geo, angstrom=True)

    def _perpendicular_direction(idxs):
        """ find a nice perpendicular direction for a series of linear atoms
        """
        triplets = []
        for idx in idxs:
            for n1idx in ngb_idxs_dct[idx]:
                for n2idx in ngb_idxs_dct[n1idx]:
                    if n2idx != idx:
                        ang = central_angle(
                            geo, idx, n1idx, n2idx, degree=True)
                        if numpy.abs(ang - 180.) > tol:
                            triplets.append((idx, n1idx, n2idx))

        if triplets:
            idx1, idx2, idx3 = min(triplets, key=lambda x: x[1:])
            xyz1, xyz2, xyz3 = map(xyzs.__getitem__, (idx1, idx2, idx3))
            r12 = util.vec.unit_direction(xyz1, xyz2)
            r23 = util.vec.unit_direction(xyz2, xyz3)
            direc = util.vec.orthogonalize(r12, r23, normalize=True)
        else:
            if len(idxs) > 1:
                idx1, idx2 = idxs[:2]
            else:
                idx1, = idxs
                idx2, = ngb_idxs_dct[idx1]

            xyz1, xyz2 = map(xyzs.__getitem__, (idx1, idx2))
            r12 = util.vec.unit_direction(xyz1, xyz2)
            for i in range(3):
                disp = numpy.zeros((3,))
                disp[i] = -1.
                alt = numpy.add(r12, disp)
                direc = util.vec.unit_perpendicular(r12, alt)
                if numpy.linalg.norm(direc) > 1e-2:
                    break

        return direc

    # partition the linear atoms into adjacent groups, to be handled together
    lin_idxs_lst = sorted(map(sorted, util.equivalence_partition(
        lin_idxs, lambda x, y: x in ngb_idxs_dct[y])))

    dummy_key_dct = {}

    for idxs in lin_idxs_lst:
        direc = _perpendicular_direction(idxs)
        for idx in idxs:
            xyz = numpy.add(xyzs[idx], numpy.multiply(dist, direc))
            dummy_key_dct[idx] = count(geo)

            geo = insert(geo, 'X', xyz, angstrom=True)

    return geo, dummy_key_dct


def linear_atoms(geo, gra=None, tol=5.):
    """ find linear atoms in a geometry (atoms with 180 degree bond angle)

        :param geo: the geometry
        :type geo: automol geometry data structure
        :param gra: the graph describing connectivity; if None, a connectivity
            graph will be generated using default distance thresholds
        :type gra: automol graph data structure
        :param tol: the tolerance threshold for linearity, in degrees
        :type tol: float
        :rtype: tuple(int)
    """

    gra = connectivity_graph(geo) if gra is None else gra
    ngb_idxs_dct = atoms_neighbor_atom_keys(gra)

    lin_idxs = []
    for idx in range(count(geo)):
        nidxs = ngb_idxs_dct[idx]
        if len(nidxs) >= 2:
            for nidx1, nidx2 in itertools.combinations(nidxs, 2):
                ang = central_angle(geo, nidx1, idx, nidx2, degree=True)
                if numpy.abs(ang - 180.) < tol:
                    lin_idxs.append(idx)

    lin_idxs = tuple(lin_idxs)

    return lin_idxs


def central_angle(geo, idx1, idx2, idx3, degree=False):
    """ Measure the angle inscribed by three atoms in a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param idx1: index of atom 1 in the triplet to be measured
        :type idx1: int
        :param idx2: index of atom 2 in the triplet to be measured
        :type idx2: int
        :param idx3: index of atom 3 in the triplet to be measured
        :type idx3: int
        :param degree: parameter to control conversion to degree
        :type degree: bool
        :rtype: float
    """

    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    xyz3 = xyzs[idx3]
    ang = util.vec.central_angle(xyz1, xyz2, xyz3)
    ang *= phycon.RAD2DEG if degree else 1
    return ang


def insert(geo, symb, xyz, idx=None, angstrom=False):
    """ Insert an atom into a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param symb: symbol of atom to add
        :type symb: str
        :param xyz: xyz coordinates of atom to add
        :type xyz: tuple(float)
        :param idx: index of geometry to place atom
        :type idx: int
        :rtype: automol geometry date structure
    """

    symbs = list(symbols(geo))
    xyzs = list(coordinates(geo, angstrom=angstrom))

    idx = idx if idx is not None else len(symbs)

    symbs.insert(idx, symb)
    xyzs.insert(idx, xyz)

    return automol.create.geom.from_data(symbs, xyzs, angstrom=angstrom)


def without_dummy_atoms(geo):
    """ Return a copy of the molecular geometry without dummy atoms.
    """

    symbs = symbols(geo)
    non_dummy_idxs = [idx for idx, symb in enumerate(symbs)
                      if ptab.to_number(symb)]

    return from_subset(geo, non_dummy_idxs)


def geometry(ich):
    """ Generate a molecular geometry from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: automol molecular geometry data structure
    """

    # rdkit fails for multi-component inchis, so we split it up and space out
    # the geometries
    ichs = _split(ich)
    geos = list(map(_connected_geometry, ichs))
    geos = [translate(geo, [50. * idx, 0., 0.])
            for idx, geo in enumerate(geos)]
    geo = functools.reduce(geometry_join, geos)
    return geo


def _connected_geometry(ich):
    """ Generate a molecular geometry from an InChI string where
        all atoms are connected by at least one bond.

        :param ich: InChI string
        :type ich: str
        :rtype: automol molecular geometry data structure
    """

    geo = object_from_hardcoded_inchi_by_key('geom', ich)
    if geo is None:
        ich = standard_form(ich)

        def _gen1(ich):
            rdm = _rdkit.from_inchi(ich)
            geo, = _rdkit.to_conformers(rdm, nconfs=1)
            return geo

        def _gen2(ich):
            pbm = _pybel.from_inchi(ich)
            geo = _pybel.to_geometry(pbm)
            return geo

        def _gen3(ich):
            if has_stereo(ich):
                raise ValueError

            gra = inchi_graph(ich, stereo=False)
            gra = explicit(gra)
            geo = automol.graph.embed.geometry(gra)
            return geo

        for gen_ in [_gen1, _gen1, _gen1, _gen2, _gen3]:
            success = False
            # geo = gen_(ich)
            # geo_ich = automol.convert.geom.inchi(geo)
            # same_conn = automol.inchi.same_connectivity(ich, geo_ich)
            # conn = automol.geom.connected(geo)
            # has_stereo = automol.inchi.has_stereo(ich)
            try:
                geo = gen_(ich)
                geo_ich = inchi(geo)
                # Check connectivity
                same_conn = same_connectivity(ich, geo_ich)
                conn = connected(geo)
                _has_stereo = has_stereo(ich)
                ich_equiv = equivalent(ich, geo_ich)
                if (same_conn and conn) and (not _has_stereo or ich_equiv):
                    success = True
                    break
            except (RuntimeError, TypeError, ValueError):
                continue

        if not success:
            raise error.FailedGeometryGenerationError

    return geo


def conformers(ich, nconfs):
    """ Generate a list of molecular geometries for various conformers
        of a species from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param nconfs: number of conformer structures to generate
        :type: int
        :rtype: automol molecular geometry data structure
    """

    geo = object_from_hardcoded_inchi_by_key('geom', ich)
    if geo is None:
        ich = standard_form(ich)

        def _gen1(ich):
            rdm = _rdkit.from_inchi(ich)
            geos = _rdkit.to_conformers(rdm, nconfs)
            return geos

        # def _gen2(ich):
        #     pbm = _pybel.from_inchi(ich)
        #     geos = _pybel.to_conformers(pbm)
        #     return geos

        for gen_ in [_gen1]:
            success = False
            try:
                geos = gen_(ich)
                for geo in geos:
                    geo_ich = inchi(geo)
                    if same_connectivity(ich, geo_ich) and (
                            not has_stereo(ich) or
                            equivalent(ich, geo_ich)):
                        success = True  # fix
                        break
            except (RuntimeError, TypeError, ValueError):
                continue

        if not success:
            raise error.FailedGeometryGenerationError

    return geos


def components_graph(geo, stereo=True):
    """ Generate a list of molecular graphs where each element is a graph that
        consists of fully connected (bonded) atoms. Stereochemistry is included
        if requested.
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph data structure
    """
    return connected_components(graph(geo, stereo=stereo))


def connected(geo, stereo=True):
    """ Determine if all atoms in geometry are completely connected.
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: bool
    """
    return len(components_graph(geo, stereo=stereo)) == 1


def dihedral_angle(geo, idx1, idx2, idx3, idx4, degree=False):
    """ Measure the angle inscribed by three atoms in a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param idx1: index of atom 1 in the quartet to be measured
        :type idx1: int
        :param idx2: index of atom 2 in the quartet to be measured
        :type idx2: int
        :param idx3: index of atom 3 in the quartet to be measured
        :type idx3: int
        :param idx4: index of atom 4 in the quartet to be measured
        :type idx4: int
        :param degree: parameter to control conversion to degree
        :type degree: bool
        :rtype: float
    """

    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    xyz3 = xyzs[idx3]
    xyz4 = xyzs[idx4]
    dih = util.vec.dihedral_angle(xyz1, xyz2, xyz3, xyz4)
    dih *= phycon.RAD2DEG if degree else 1

    return dih


def distance(geo, idx1, idx2, angstrom=False):
    """ Measure the distance between two atoms in a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param idx1: index of atom 1 in the pair to be measured
        :type idx1: int
        :param idx2: index of atom 2 in the pair to be measured
        :type idx2: int
        :param angstrom: parameter to control conversion to Angstrom
        :type angstrom: bool
        :rtype: float
    """

    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    dist = util.vec.distance(xyz1, xyz2)
    dist *= phycon.BOHR2ANG if angstrom else 1
    return dist


def inchi_graph(ich, stereo=True):
    """ Generate a molecular graph from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph
    """

    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = _split(ich)
    gras = [_inchi_connected_graph(ich, stereo=stereo) for ich in ichs]
    for idx, gra in enumerate(gras):
        if idx == 0:
            num = 0
        else:
            num = max(map(max, map(atom_keys, gras[:idx]))) + 1
        gras[idx] = transform_keys(gra, num.__add__)
    gra = functools.reduce(union, gras)
    return gra


def _inchi_connected_graph(ich, stereo=True):
    """ Generate a molecular graph from an InChI string where
        all all atoms are connected by at least one bond.

        :param ich: InChI string
  d      :type ich: str
        :param remove_stereo: parameter to include stereochemistry information
        :type remove_stereo: bool
        :rtype: automol molecular graph
    """

    gra = object_from_hardcoded_inchi_by_key('graph', ich)
    if gra is None:
        ich = standard_form(ich)
        if not stereo or not has_stereo(ich):
            rdm = _rdkit.from_inchi(ich)
            gra = _rdkit.to_connectivity_graph(rdm)
        else:
            geo = geometry(ich)
            gra = graph(geo, stereo=stereo)

    gra = implicit(gra)
    return gra
